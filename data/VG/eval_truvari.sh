#!/bin/bash
set -euo pipefail

# Evaluate VG results: truvari bench (F1/Precision/Recall)

REF="$(readlink -f ../hprc_haps/MHC-CHM13.0.fa)"
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
  for T in "${TAGS[@]}"; do
    BASE_VCF="${TEST_DIR}/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    TRUTH_VCF="${TRUTH_DIR}/${S}/MHC_${S}.vcf.gz"
    OUTDIR="${OUTROOT}/${S}/${S}_${T}"
    TMPD="${OUTDIR}/tmp"
    mkdir -p "$OUTDIR" "$TMPD"

    [[ -f "$BASE_VCF" ]]  || { echo "[WARN] Missing: $BASE_VCF"; continue; }
    [[ -f "$TRUTH_VCF" ]] || { echo "[WARN] Missing: $TRUTH_VCF"; continue; }

    BASE_SPLIT="${OUTDIR}/base.split.vcf.gz"
    TRUTH_SPLIT="${OUTDIR}/truth.split.vcf.gz"

    echo "[NORM] $S $T (base)"
    bcftools norm -m -both -f "$REF" "$BASE_VCF" \
      | bcftools view -i 'GT="alt"' \
      | bcftools sort -T "$TMPD" -Oz -o "$BASE_SPLIT"

    echo "[NORM] $S $T (truth)"
    bcftools norm -m -both -f "$REF" "$TRUTH_VCF" \
      | bcftools sort -T "$TMPD" -Oz -o "$TRUTH_SPLIT"

    tabix -f -p vcf "$BASE_SPLIT"
    tabix -f -p vcf "$TRUTH_SPLIT"

    BENCH_DIR="${OUTDIR}/bench"
    rm -rf "$BENCH_DIR"
    echo "[TRUVARI] $S $T"
    truvari bench -f "$REF" -b "$BASE_SPLIT" -c "$TRUTH_SPLIT" -o "$BENCH_DIR"

    rm -rf "$TMPD"
  done
done

echo "Done. Results in ${OUTROOT}/<SAMPLE>/<SAMPLE>_<TAG>/bench/"
