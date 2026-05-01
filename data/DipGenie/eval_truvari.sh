#!/bin/bash
set -euo pipefail

# Evaluate DipGenie results: truvari bench (F1/Precision/Recall)

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)

HPRC_HAPS="../hprc_haps"
REF="${HPRC_HAPS}/MHC-CHM13.0.fa"
OUTROOT="SV_Evaluation"
TRUTH_DIR="../Truth"
TEST_DIR="Results"

mkdir -p "$OUTROOT"
shopt -s nullglob

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

for S in "${SAMPLES[@]}"; do
  TRUTH="${TRUTH_DIR}/${S}/MHC_${S}.vcf.gz"
  [[ -f "$TRUTH" ]] || { echo "[WARN] Missing truth VCF: $TRUTH"; continue; }

  for T in "${TAGS[@]}"; do
    TEST="${TEST_DIR}/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    OUTDIR="${OUTROOT}/${S}/${S}_${T}"
    TMPD="${OUTDIR}/tmp"
    mkdir -p "$OUTDIR" "$TMPD"

    [[ -f "$TEST" ]] || { echo "[WARN] Missing: $TEST"; continue; }

    TEST_SPLIT="${OUTDIR}/test.split.vcf.gz"
    TRUTH_SPLIT="${OUTDIR}/truth.split.vcf.gz"

    echo "[NORM] $S $T (test)"
    bcftools norm -m -both -f "$REF" "$TEST" \
      | bcftools view -i 'GT="alt"' \
      | bcftools sort -T "$TMPD" -Oz -o "$TEST_SPLIT"

    echo "[NORM] $S $T (truth)"
    bcftools norm -m -both -f "$REF" "$TRUTH" \
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
