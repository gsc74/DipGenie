#!/usr/bin/env bash
# hap.py SNP error classes for every evaluated run of one mode.
# Usage: snp_errors.sh [mixed|SR]   -> $WORK/snp-errors/<mode>.tsv
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/07_snp_errors
mode=${1:-mixed}
suffix=$([[ $mode == SR ]] && echo -SR || true)
runs=()
for depth in 2x 4x 8x 16x 32x Total; do
  runs+=("$WORK/runs/VG-HapCUT2$suffix-$depth" "$WORK/runs/PanGenie-HapCUT2$suffix-$depth"
         "$WORK/runs/${LABEL:-DipGenie}/$mode/$depth")
done
mkdir -p "$WORK/snp-errors"
"$PYTHON" "$here/snp_errors.py" --bcftools "$BCFTOOLS" --jobs 18 "${runs[@]}" > "$WORK/snp-errors/$mode.tsv"
