#!/bin/bash
set -euo pipefail
shopt -s nullglob

echo -e "Sample\tDepth\tSwitchErrorRate"

for depth in 2x 4x full; do
  for ser in Evaluation/*/*_${depth}/SER.txt; do
    sample="$(basename "$(dirname "$(dirname "$ser")")")"   # e.g. HG002
    rate="$(grep -m1 'switch error rate' "$ser" | awk '{print $NF}')"
    [[ -n "$rate" ]] || rate="NA"
    echo -e "${sample}\t${depth}\t${rate}"
  done
done