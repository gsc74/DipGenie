#!/bin/bash
set -euo pipefail

# Collect VG results into CSV files

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)

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
    echo "${S},${T},${f1},VG" >> output/F1.csv
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
    echo "${S},${T},${val},VG" >> output/Precision.csv
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
    echo "${S},${T},${val},VG" >> output/Recall.csv
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
    echo "${S},${T},${rate},VG" >> output/SER.csv
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
    echo "${S},${T},${rate},VG" >> output/Mismatch.csv
  done
done

# --- SV Counts ---
echo "Sample,Depth,Value,Tool" > output/SV_counts.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    VCF="Results/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    if [[ -f "$VCF" ]]; then
      count=$(bcftools view -i 'GT="alt"' "$VCF" | bcftools query -f '%REF\t%ALT\n' \
        | awk '{split($2,a,","); for(i in a) if (length(a[i])-length($1)>=50 || length($1)-length(a[i])>=50) c++} END{print c+0}')
    else
      count="--"
    fi
    echo "${S},${T},${count},VG" >> output/SV_counts.csv
  done
done

# --- Runtime and Memory from logs (parsed from /usr/bin/time -v output) ---
echo "Sample,Depth,Tool,Value" > output/Runtime.csv
echo "Sample,Depth,Tool,Value" > output/Memory.csv
for T in "${TAGS[@]}"; do
  for S in "${SAMPLES[@]}"; do
    LOG="logs/${S}_${T}.err"
    if [[ -f "$LOG" ]]; then
      # /usr/bin/time -v format
      rt=$(grep 'Elapsed (wall clock)' "$LOG" | sed -E 's/.*: (.*)/\1/' | tail -1) || true
      mem=$(grep 'Maximum resident set size' "$LOG" | awk '{print $NF}' | tail -1) || true
      # Convert wall clock h:mm:ss or m:ss to seconds
      if [[ -n "$rt" ]]; then
        rt=$(echo "$rt" | awk -F: '{if(NF==3) print $1*3600+$2*60+$3; else if(NF==2) print $1*60+$2; else print $1}')
      else
        rt="--"
      fi
      # Convert KB to GB
      if [[ -n "$mem" && "$mem" != "--" ]]; then
        mem=$(awk "BEGIN{printf \"%.3f\", $mem/1048576}")
      else
        mem="--"
      fi
    else
      rt="--"
      mem="--"
    fi
    echo "${S},${T},VG,${rt}" >> output/Runtime.csv
    echo "${S},${T},VG,${mem}" >> output/Memory.csv
  done
done

echo "All results saved to output/"
