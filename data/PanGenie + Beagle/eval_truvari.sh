#!/bin/bash
set -euo pipefail

# Evaluate PanGenie+Beagle results: truvari bench (F1/Precision/Recall)

REF_DIR="Reference"
OUTROOT="SV_Evaluation"
TRUTH_DIR="../Truth"
TEST_DIR="Results"

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)

mkdir -p "$OUTROOT"
shopt -s nullglob

for S in "${SAMPLES[@]}"; do
  REF="${REF_DIR}/${S}.fa"
  [[ -f "$REF" ]] || { echo "[WARN] Missing reference: $REF"; continue; }
  [[ -f "${REF}.fai" ]] || samtools faidx "$REF"
  CONTIG_LEN=$(awk '{print $2}' "${REF}.fai")

  TRUTH_ORIG="${TRUTH_DIR}/${S}/MHC_${S}.vcf.gz"
  [[ -f "$TRUTH_ORIG" ]] || { echo "[WARN] Missing truth VCF: $TRUTH_ORIG"; continue; }
  mkdir -p "${OUTROOT}/${S}"
  TRUTH_RENAMED="${OUTROOT}/${S}/truth_renamed.vcf.gz"
  bcftools annotate --rename-chrs <(echo "0 ${S}") "$TRUTH_ORIG" -Oz -o "$TRUTH_RENAMED"
  bcftools index -t -f "$TRUTH_RENAMED"

  for T in "${TAGS[@]}"; do
    TEST="${TEST_DIR}/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    OUTDIR="${OUTROOT}/${S}/${S}_${T}"
    TMPD="${OUTDIR}/tmp"
    mkdir -p "$OUTDIR" "$TMPD"

    [[ -f "$TEST" ]] || { echo "[WARN] Missing: $TEST"; continue; }

    TEST_SPLIT="${OUTDIR}/test.split.vcf.gz"
    TRUTH_SPLIT="${OUTDIR}/truth.split.vcf.gz"

    fix_contig_hdr() {
      local invcf="$1"
      bcftools view -h "$invcf" \
        | sed "s|##contig=<ID=${S}>|##contig=<ID=${S},length=${CONTIG_LEN}>|"
    }

    echo "[NORM] $S $T (test)"
    bcftools reheader -h <(fix_contig_hdr "$TEST") "$TEST" \
      | bcftools norm -m -both -f "$REF" \
      | bcftools view -i 'GT="alt"' \
      | bcftools sort -T "$TMPD" -Oz -o "$TEST_SPLIT"

    echo "[NORM] $S $T (truth)"
    bcftools reheader -h <(fix_contig_hdr "$TRUTH_RENAMED") "$TRUTH_RENAMED" \
      | bcftools norm -m -both -f "$REF" \
      | bcftools sort -T "$TMPD" -Oz -o "$TRUTH_SPLIT"

    tabix -f -p vcf "$TEST_SPLIT"
    tabix -f -p vcf "$TRUTH_SPLIT"

    BENCH_DIR="${OUTDIR}/bench"
    rm -rf "$BENCH_DIR"
    echo "[TRUVARI] $S $T"
    truvari bench -f "$REF" -b "$TRUTH_SPLIT" -c "$TEST_SPLIT" -o "$BENCH_DIR"

    rm -rf "$TMPD"
  done
done

echo "Done."
