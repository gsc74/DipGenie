#!/bin/bash
set -euo pipefail
shopt -s nullglob

echo -e "Sample\tDepth\tF1_Score"

for depth in 2x 4x full; do
  for log in SV_Evaluation/*/*_${depth}/bench/log.txt; do
    sample="$(basename "$(dirname "$(dirname "$log")")")"   # e.g. HG002
    f1="$(grep -m1 '"recall"' "$log" | awk -F': ' '{print $2}' | tr -d ', ')"
    [[ -n "$f1" ]] || f1="NA"
    echo -e "${sample}\t${depth}\t${f1}"
  done
done
