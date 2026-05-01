#!/bin/bash
set -euo pipefail

# Collect DipGenie results into CSV files

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)
TOOL="DipGenie"

mkdir -p output

# --- F1-score ---
echo "Sample,Depth,Value,Tool" > output/F1.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG="SV_Evaluation/${S}/${S}_${T}/bench/summary.json"
    if [[ -f "$LOG" ]]; then
      f1=$(python3 -c "import json; d=json.load(open('$LOG')); print(d.get('f1', 'NA'))")
    else
      f1="--"
    fi
    echo "${S},${T},${f1},${TOOL}" >> output/F1.csv
  done
done

# --- Precision ---
echo "Sample,Depth,Value,Tool" > output/Precision.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG="SV_Evaluation/${S}/${S}_${T}/bench/summary.json"
    if [[ -f "$LOG" ]]; then
      val=$(python3 -c "import json; d=json.load(open('$LOG')); print(d.get('precision', 'NA'))")
    else
      val="--"
    fi
    echo "${S},${T},${val},${TOOL}" >> output/Precision.csv
  done
done

# --- Recall ---
echo "Sample,Depth,Value,Tool" > output/Recall.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG="SV_Evaluation/${S}/${S}_${T}/bench/summary.json"
    if [[ -f "$LOG" ]]; then
      val=$(python3 -c "import json; d=json.load(open('$LOG')); print(d.get('recall', 'NA'))")
    else
      val="--"
    fi
    echo "${S},${T},${val},${TOOL}" >> output/Recall.csv
  done
done

# --- SER ---
echo "Sample,Depth,Value,Tool" > output/SER.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    SER_FILE="Evaluation/${S}/${S}_${T}/SER.txt"
    if [[ -f "$SER_FILE" ]]; then
      rate=$(grep -m1 'switch error rate' "$SER_FILE" | awk '{print $NF}') || true
      [[ -n "$rate" ]] || rate="--"
    else
      rate="--"
    fi
    echo "${S},${T},${rate},${TOOL}" >> output/SER.csv
  done
done

# --- Mismatch Rate ---
echo "Sample,Depth,Value,Tool" > output/Mismatch.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    SER_FILE="Evaluation/${S}/${S}_${T}/SER.txt"
    if [[ -f "$SER_FILE" ]]; then
      rate=$(grep -m1 'Hamming distance \[%\]' "$SER_FILE" | grep 'Block-wise' | awk '{print $NF}') || true
      [[ -n "$rate" ]] || rate="--"
    else
      rate="--"
    fi
    echo "${S},${T},${rate},${TOOL}" >> output/Mismatch.csv
  done
done

# --- Recombinations (from DipGenie stderr) ---
echo "Sample,Depth,Value,Tool" > output/Recombinations.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG="logs/${S}_${T}.err"
    if [[ -f "$LOG" ]]; then
      val=$(grep -m1 'recombinations' "$LOG" | awk '{print $NF}') || true
      [[ -n "$val" ]] || val="--"
    else
      val="--"
    fi
    echo "${S},${T},${val},${TOOL}" >> output/Recombinations.csv
  done
done

echo "All results saved to output/"
