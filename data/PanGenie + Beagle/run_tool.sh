#!/bin/bash
set -euo pipefail

# Usage: bash run_tool.sh
# Runs PanGenie + Beagle pipeline for all samples at all coverages

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
JELLY_HASH=3000000000

NUMA="numactl --cpunodebind=0 --membind=0"
command -v numactl >/dev/null || { echo "[FATAL] numactl not found"; exit 127; }

PANGENIE="$(pwd)/pangenie/build/src/PanGenie"
PANGENIE_INDEX="$(pwd)/pangenie/build/src/PanGenie-index"

mkdir -p logs "$OUT_BASE" "$REF_DIR" PanGenie_index

# Prepare per-sample reference
prepare_ref() {
  local S="$1"
  local REF_FA="${REF_DIR}/${S}.fa"
  if [[ ! -f "$REF_FA" ]]; then
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

  echo "[RUN] PanGenie+Beagle: $SAMPLE $TAG"

  local WORK
  WORK=$(mktemp -d "${SAMPLE}_${TAG}_pangenie.XXXXXX")
  trap 'rm -rf "${WORK}"' RETURN
  mkdir -p "${WORK}"/{reads,panel,logs}

  # Prepare PanGenie input VCF (rename contig)
  local PANEL_LOCAL="${WORK}/panel/${SAMPLE}.pangenie.vcf"
  printf "CHM13#0#0\t%s\n" "$SAMPLE" > "${WORK}/panel/rename.tsv"
  bcftools annotate --rename-chrs "${WORK}/panel/rename.tsv" -Ov -o "$PANEL_LOCAL" "$PANEL_VCF"

  # Merge R1+R2 into single FASTQ for PanGenie
  local READS_FQ="${WORK}/reads/${SAMPLE}.${TAG}.fastq"
  cat "$R1" "$R2" > "$READS_FQ"

  # Build PanGenie index (once per sample)
  local INDEX_PREFIX="PanGenie_index/${SAMPLE}"
  if [[ ! -s "${INDEX_PREFIX}_UniqueKmersMap.cereal" ]]; then
    echo "[PanGenie-index] Building index for $SAMPLE..."
    ${NUMA} "$PANGENIE_INDEX" \
      -r "$REF" -v "$PANEL_LOCAL" -t "$THREADS" -e "$JELLY_HASH" -o "$INDEX_PREFIX" \
      > "${WORK}/logs/pangenie-index.stdout.log" 2> "${WORK}/logs/pangenie-index.stderr.log"
  fi

  # Genotype
  local PG_PREFIX="${WORK}/${SAMPLE}.${TAG}"
  ${NUMA} "$PANGENIE" \
    -f "$INDEX_PREFIX" -i "$READS_FQ" -s "$SAMPLE" \
    -j "$THREADS" -t "$THREADS" -e "$JELLY_HASH" -o "$PG_PREFIX" \
    > "${WORK}/logs/pangenie.stdout.log" 2> "${WORK}/logs/pangenie.stderr.log"

  local FINAL_VCF="${PG_PREFIX}_genotyping.vcf"
  [[ -s "$FINAL_VCF" ]] || { echo "[ERROR] PanGenie failed for $SAMPLE $TAG"; return; }

  # Fix ploidy + normalize
  local PG_FIXED="${WORK}/${SAMPLE}.pangenie.fixed.vcf.gz"
  bcftools +fixploidy "$FINAL_VCF" -Ou -- -f 2 \
  | bcftools norm -m -any -Ou \
  | bcftools view -m 2 -M 2 -Oz -o "$PG_FIXED"
  bcftools index -f "$PG_FIXED"

  # Prepare biallelic panel for Beagle
  local PANEL_GZ="${WORK}/panel/panel.phased.vcf.gz"
  bcftools norm -m -any "$PANEL_LOCAL" -Ou \
  | bcftools view -m 2 -M 2 -Oz -o "$PANEL_GZ"
  bcftools index -f "$PANEL_GZ"

  # Beagle phasing (impute=false)
  local BEAGLE_PREFIX="${WORK}/${SAMPLE}.${TAG}.phased"
  local GMAP_PLINK="${WORK}/${SAMPLE}.beagle.map"
  awk -v chr="$SAMPLE" 'BEGIN{OFS="\t"} NR>1 {print chr, ".", $3, $1+1}' "$GMAP" > "$GMAP_PLINK"

  ${NUMA} beagle -Xmx16g \
    ref="$PANEL_GZ" gt="$PG_FIXED" out="$BEAGLE_PREFIX" \
    map="$GMAP_PLINK" impute=false nthreads="$THREADS" ne=5000 \
    > "${WORK}/logs/beagle.stdout.log" 2> "${WORK}/logs/beagle.stderr.log" \
  || { echo "[ERROR] Beagle failed for $SAMPLE $TAG"; return; }

  local PHASED_VCF="${BEAGLE_PREFIX}.vcf.gz"
  bcftools index -f "$PHASED_VCF"

  cp -f "$PHASED_VCF" "${OUT_DIR}/MHC_${SAMPLE}_${TAG}.vcf.gz"
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
