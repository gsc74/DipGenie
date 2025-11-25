#!/bin/bash
set -euo pipefail

THREADS=128
GRAPH_DIR="Graph"
READS_DIR="Reads"
OUT_BASE="Results"
DipGenie="./dipGenie/DipGenie"

mkdir -p logs "$OUT_BASE"
command -v numactl >/dev/null || { echo "[FATAL] numactl not found"; exit 127; }
[ -x "$DipGenie" ] || { echo "[FATAL] DipGenie not executable at $DipGenie"; exit 127; }

shopt -s nullglob

while read -r SAMPLE; do
  [[ -z "$SAMPLE" ]] && continue
  echo "[INFO] $SAMPLE"

  # Resolve inputs
  GFA=( "$GRAPH_DIR"/*_"$SAMPLE".gfa )
  R2x=( "$READS_DIR/2x/${SAMPLE}"*.fastq )
  R4x=( "$READS_DIR/4x/${SAMPLE}"*.fastq )
  Rfull=( "$READS_DIR/Total/${SAMPLE}"*.fastq )

  [[ ${#GFA[@]} -ge 1 ]] || { echo "[WARN] GFA missing for $SAMPLE"; continue; }
  [[ -f "${R2x[0]:-}" ]] || echo "[WARN] 2x reads missing for $SAMPLE"
  [[ -f "${R4x[0]:-}" ]] || echo "[WARN] 4x reads missing for $SAMPLE"
  [[ -f "${Rfull[0]:-}" ]] || echo "[WARN] Total reads missing for $SAMPLE"

  run_one () {
    local READ="$1" TAG="$2"
    [[ -f "$READ" ]] || return 0
    local OUT_DIR="$OUT_BASE/$SAMPLE/${SAMPLE}_${TAG}"
    local OUT_FA="$OUT_DIR/full.fa"
    mkdir -p "$OUT_DIR"
    echo "[RUN] $SAMPLE $TAG -> $OUT_FA"
    numactl --cpunodebind=0 --membind=0 \
      "$DipGenie" -t"$THREADS" -a1 -R18 -g "${GFA[0]}" -r "$READ" -o "$OUT_FA" \
      > "logs/${SAMPLE}_${TAG}.out" \
      2> "logs/${SAMPLE}_${TAG}.err"
  }

  # run_one "${R2x[0]:-}"  "2x"
  # run_one "${R4x[0]:-}"  "4x"
  run_one "${Rfull[0]:-}" "full"

done < "${1:-/dev/stdin}"

