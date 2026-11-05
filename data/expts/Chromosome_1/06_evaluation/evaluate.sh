#!/usr/bin/env bash
# Evaluate one call set against GIAB v5.0q chr1 truth (same protocol for every method).
# Usage: evaluate.sh calls.vcf.gz output_dir
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/06_evaluation
dl=$WORK/download
mkdir -p "$(dirname "$2")"
"$PYTHON" "$here/evaluate.py" --query "$1" --reference "$dl/reference/CHM13v2.fa" \
  --truth-small "$WORK/truth/smvar.vcf.gz" --truth-sv "$WORK/truth/stvar.vcf.gz" \
  --bed-small "$WORK/truth/smvar.bed" --bed-sv "$WORK/truth/stvar.bed" \
  --sdf "$dl/reference/CHM13v2.sdf" --stratification "$dl/stratifications/CHM13-all-stratifications.tsv" \
  --sample "$SAMPLE" --gender male --threads "$THREADS" --rtg "$RTG" --truvari "$TRUVARI" --out "$2" > "$2.log" 2>&1
