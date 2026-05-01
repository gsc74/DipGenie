#!/bin/bash
set -euo pipefail

# Usage: bash run_tool.sh

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)

THREADS=128
GRAPH_DIR="../Graphs"
READS_DIR="../Reads"
OUT_BASE="Results"
DIPGENIE="../../DipGenie/DipGenie"

mkdir -p logs "$OUT_BASE"
command -v numactl >/dev/null || { echo "[FATAL] numactl not found"; exit 127; }
[[ -x "$DIPGENIE" ]] || { echo "[FATAL] DipGenie not executable at $DIPGENIE"; exit 127; }

NUMA="numactl --cpunodebind=0 --membind=0"

shopt -s nullglob

for SAMPLE in "${SAMPLES[@]}"; do
  echo "[INFO] Processing $SAMPLE"

  GFA="${GRAPH_DIR}/MHC_left_${SAMPLE}.gfa"
  [[ -f "$GFA" ]] || { echo "[WARN] GFA missing for $SAMPLE: $GFA"; continue; }

  run_one () {
    local READ="$1" TAG="$2"
    [[ -f "$READ" ]] || { echo "[WARN] Reads missing: $READ"; return 0; }
    local OUT_DIR="$OUT_BASE/$SAMPLE/${SAMPLE}_${TAG}"
    local OUT_FA="$OUT_DIR/full.fa"
    mkdir -p "$OUT_DIR"
    echo "[RUN] $SAMPLE $TAG -> $OUT_FA"
    ${NUMA} "$DIPGENIE" -t"$THREADS" -p2 -R18 \
      -g "$GFA" \
      -r "$READ" \
      -o "$OUT_FA" \
      > "logs/${SAMPLE}_${TAG}.out" \
      2> "logs/${SAMPLE}_${TAG}.err"
  }

  run_one "${READS_DIR}/1x/${SAMPLE}.mhc.1x.fastq"       "1x"
  run_one "${READS_DIR}/2x/${SAMPLE}.mhc.2x.fastq"       "2x"
  run_one "${READS_DIR}/Total/${SAMPLE}.mhc.Total.fastq"  "Total"

done

echo "Done."
