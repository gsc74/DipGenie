#!/bin/bash
set -euo pipefail

# Usage: bash run_tool.sh
# Runs Paragraph + Beagle pipeline for all samples at all coverages

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)

THREADS=128
GRAPH_DIR="../Graphs"
READS_DIR="../Reads"
VCF_DIR="../VCFs"
REF_DIR="Reference"
OUT_BASE="Results"
GMAP="../MHC_CHM13.gmap"

NUMA="numactl --cpunodebind=0 --membind=0"
command -v numactl >/dev/null || { echo "[FATAL] numactl not found"; exit 127; }

mkdir -p logs "$OUT_BASE" "$REF_DIR"

# Prepare per-sample reference files (CHM13 renamed per sample for contig matching)
prepare_ref() {
  local S="$1"
  local REF_FA="${REF_DIR}/${S}.fa"
  if [[ ! -f "$REF_FA" ]]; then
    # Rename CHM13 contig to SAMPLE name for bcftools compatibility
    awk -v S="$S" '/^>/{print ">" S; next} {print}' ../hprc_haps/MHC-CHM13.0.fa > "$REF_FA"
    samtools faidx "$REF_FA"
  fi
}

run_one() {
  local SAMPLE="$1" TAG="$2" R1="$3" R2="$4"
  [[ -f "$R1" ]] || { echo "[WARN] Missing R1: $R1"; return; }
  [[ -f "$R2" ]] || { echo "[WARN] Missing R2: $R2"; return; }

  local REF="${REF_DIR}/${SAMPLE}.fa"
  local PANEL_VCF="${VCF_DIR}/${SAMPLE}.vcf"
  local CONTIG="${SAMPLE}"
  local OUT_DIR="${OUT_BASE}/${SAMPLE}/${SAMPLE}_${TAG}"
  mkdir -p "$OUT_DIR"

  [[ -f "$PANEL_VCF" ]] || { echo "[WARN] Missing panel VCF: $PANEL_VCF"; return; }

  echo "[RUN] Paragraph+Beagle: $SAMPLE $TAG"

  local WORK
  WORK=$(mktemp -d "${SAMPLE}_${TAG}_impute.XXXXXX")
  trap 'rm -rf "${WORK}"' RETURN
  mkdir -p "${WORK}"/{ref,align,panel,phased,logs}

  # Prepare reference
  cp -f "$REF" "${WORK}/ref/${SAMPLE}.fa"
  local REF_LOCAL="${WORK}/ref/${SAMPLE}.fa"
  samtools faidx "$REF_LOCAL"

  # BWA-MEM2 index
  if [[ ! -f "${REF}.bwt.2bit.64" ]]; then
    ${NUMA} bwa-mem2 index "$REF" 2> "${WORK}/logs/${SAMPLE}.bwa-index.log"
  fi

  # Align
  local BAM="${WORK}/align/${SAMPLE}.sorted.bam"
  ${NUMA} bwa-mem2 mem -t "$THREADS" \
    -R "@RG\\tID:${SAMPLE}_${TAG}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:lib1" \
    "$REF" "$R1" "$R2" 2> "${WORK}/logs/${SAMPLE}.bwa-mem2.log" \
  | ${NUMA} samtools sort -@ "$THREADS" -O bam -o "$BAM"
  ${NUMA} samtools index -@ "$THREADS" "$BAM"

  # Prepare panel
  local PANEL_GZ="${WORK}/panel/panel.phased.vcf.gz"
  local PANEL_SITES="${WORK}/panel/panel.sites.vcf.gz"

  ${NUMA} bcftools annotate --rename-chrs <(echo "CHM13#0#0 ${SAMPLE}") "$PANEL_VCF" -Ou \
  | ${NUMA} bcftools norm -m -any -Ou \
  | ${NUMA} bcftools view -m 2 -M 2 -Oz -o "$PANEL_GZ"
  ${NUMA} bcftools index -f "$PANEL_GZ"

  # Sites (SNPs only for mpileup)
  ${NUMA} bcftools view -G -v snps -Oz -o "$PANEL_SITES" "$PANEL_GZ"
  ${NUMA} bcftools index -f "$PANEL_SITES"

  # GL generation
  local SITES_TSV="${WORK}/panel/panel.sites.tsv.gz"
  ${NUMA} bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$PANEL_SITES" | bgzip -c > "$SITES_TSV"
  tabix -s1 -b2 -e2 "$SITES_TSV"

  local GL_VCF="${WORK}/align/${SAMPLE}.gl.vcf.gz"
  ${NUMA} bcftools mpileup -f "$REF_LOCAL" -I -E -a 'FORMAT/DP' -T "$PANEL_SITES" -r "$CONTIG" "$BAM" -Ou \
  | ${NUMA} bcftools call -Aim -a GQ -C alleles -T "$SITES_TSV" -Ou \
  | ${NUMA} bcftools +fixploidy -Ou -- -f 2 \
  | ${NUMA} bcftools filter -S . -e 'FMT/GQ<10' -Oz -o "$GL_VCF"
  ${NUMA} bcftools index -f "$GL_VCF"

  # Paragraph SV genotyping
  local PANEL_SVS="${WORK}/panel/panel.svs.vcf.gz"
  ${NUMA} bcftools view -i 'STRLEN(ALT)>=50 || STRLEN(REF)>=50' -Oz -o "$PANEL_SVS" "$PANEL_GZ"
  ${NUMA} bcftools index -f "$PANEL_SVS"

  local SV_COUNT
  SV_COUNT=$(bcftools view -H "$PANEL_SVS" | wc -l)
  local SV_VCF="${WORK}/align/${SAMPLE}.sv.vcf.gz"

  if [[ "$SV_COUNT" -gt 0 ]]; then
    local MEAN_DP
    MEAN_DP=$(samtools depth -a "$BAM" | awk '{s+=$3; n++} END{if(n>0) printf "%.1f", s/n; else print "0"}')
    local READ_LEN
    READ_LEN=$(samtools stats "$BAM" | awk -F'\t' '/^SN\tmaximum length:/{print $3}')

    local MANIFEST="${WORK}/align/manifest.txt"
    printf "id\tpath\tdepth\tread length\n" > "$MANIFEST"
    printf "%s\t%s\t%s\t%s\n" "$SAMPLE" "$(realpath "$BAM")" "$MEAN_DP" "$READ_LEN" >> "$MANIFEST"

    local PARAGRAPH_OUT="${WORK}/align/paragraph_out"
    (
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate paragraph
      ${NUMA} multigrmpy.py -i "$PANEL_SVS" -m "$MANIFEST" -r "$REF_LOCAL" \
        -o "$PARAGRAPH_OUT" --threads "$THREADS" 2> "${WORK}/logs/paragraph.stderr.log"
      conda deactivate
    )

    if [[ -f "${PARAGRAPH_OUT}/genotypes.vcf.gz" ]]; then
      echo "${SAMPLE}" > "${WORK}/align/sv_sample_name.txt"
      ${NUMA} bcftools reheader -s "${WORK}/align/sv_sample_name.txt" "${PARAGRAPH_OUT}/genotypes.vcf.gz" \
      | ${NUMA} bcftools +fixploidy -Ou -- -f 2 \
      | ${NUMA} bcftools view -m 2 -M 2 -Oz -o "$SV_VCF"
      ${NUMA} bcftools index -f "$SV_VCF"
    else
      SV_VCF=""
    fi
  else
    SV_VCF=""
  fi

  # Merge SNP + SV
  local MERGED_GL="${WORK}/align/${SAMPLE}.merged.gl.vcf.gz"
  if [[ -n "$SV_VCF" && -s "$SV_VCF" ]]; then
    ${NUMA} bcftools concat -a "$GL_VCF" "$SV_VCF" -Ou \
    | ${NUMA} bcftools sort -T "$WORK" -Oz -o "$MERGED_GL"
    ${NUMA} bcftools index -f "$MERGED_GL"
  else
    cp -f "$GL_VCF" "$MERGED_GL"
    bcftools index -f "$MERGED_GL"
  fi

  # Beagle imputation
  local BEAGLE_PREFIX="${WORK}/phased/${SAMPLE}.${CONTIG}.imputed"
  local GMAP_PLINK="${WORK}/panel/${SAMPLE}.beagle.map"
  awk -v chr="$SAMPLE" 'BEGIN{OFS="\t"} NR>1 {print chr, ".", $3, $1+1}' "$GMAP" > "$GMAP_PLINK"

  ${NUMA} beagle -Xmx16g \
    ref="$PANEL_GZ" gt="$MERGED_GL" out="$BEAGLE_PREFIX" \
    map="$GMAP_PLINK" nthreads="$THREADS" ne=5000 \
    > "${WORK}/logs/beagle.stdout.log" 2> "${WORK}/logs/beagle.stderr.log" \
  || { echo "[ERROR] Beagle failed for $SAMPLE $TAG"; cat "${WORK}/logs/beagle.stderr.log" >&2; return; }

  local FINAL_VCF="${BEAGLE_PREFIX}.vcf.gz"
  ${NUMA} bcftools index -f "$FINAL_VCF"

  cp -f "$FINAL_VCF" "${OUT_DIR}/MHC_${SAMPLE}_${TAG}.vcf.gz"
  bcftools index -t -f "${OUT_DIR}/MHC_${SAMPLE}_${TAG}.vcf.gz"

  echo "[DONE] $SAMPLE/$TAG"
}

for SAMPLE in "${SAMPLES[@]}"; do
  echo "[INFO] Processing $SAMPLE"
  prepare_ref "$SAMPLE"

  run_one "$SAMPLE" "1x"    "${READS_DIR}/1x/${SAMPLE}.mhc.1x.R1.fastq"       "${READS_DIR}/1x/${SAMPLE}.mhc.1x.R2.fastq"
  run_one "$SAMPLE" "2x"    "${READS_DIR}/2x/${SAMPLE}.mhc.2x.R1.fastq"       "${READS_DIR}/2x/${SAMPLE}.mhc.2x.R2.fastq"
  run_one "$SAMPLE" "Total" "${READS_DIR}/Total/${SAMPLE}.mhc.Total.R1.fastq"  "${READS_DIR}/Total/${SAMPLE}.mhc.Total.R2.fastq"

done

echo "All done."
