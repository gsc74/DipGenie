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
LOG_DIR="logs"
VG_HAP="./vg_haplotypes.py"
TMP_BASE="tmp_vg_hap"

mkdir -p "$OUT_BASE" "$LOG_DIR"

NUMA="numactl --cpunodebind=0 --membind=0"
command -v numactl >/dev/null || { echo "[FATAL] numactl not found"; exit 127; }

shopt -s nullglob

for SAMPLE in "${SAMPLES[@]}"; do
  echo "[INFO] Processing $SAMPLE"

  TMP_DIR="${TMP_BASE}_${SAMPLE}"
  mkdir -p "$TMP_DIR"

  GFA="${GRAPH_DIR}/MHC_left_${SAMPLE}.gfa"
  [[ -f "$GFA" ]] || { echo "[WARN] missing $GFA"; continue; }

  PREFIX="${GRAPH_DIR}/MHC_left_${SAMPLE}"
  XG="${PREFIX}.xg"
  GBWT="${PREFIX}.gbwt"
  GBZ="${PREFIX}.gbz"
  DIST="${PREFIX}.dist"
  RI="${PREFIX}.ri"

  # Build VG indexes if GBZ missing
  if [[ ! -f "$GBZ" ]]; then
    echo "[IDX] ${SAMPLE}: GFA -> GBWT + XG -> GBZ"
    ${NUMA} gfa2gbwt -b "$PREFIX" < "$GFA" \
      > "$LOG_DIR/${SAMPLE}_gfa2gbwt.out" 2> "$LOG_DIR/${SAMPLE}_gfa2gbwt.err"
    ${NUMA} vg convert -g "$GFA" -x > "$XG"
    ${NUMA} vg gbwt -x "$XG" "$GBWT" --gbz-format -g "$GBZ" \
      > "$LOG_DIR/${SAMPLE}_gbz.out" 2> "$LOG_DIR/${SAMPLE}_gbz.err"
  fi

  [[ -f "$DIST" ]] || ${NUMA} vg index -j "$DIST" "$GBZ" \
    > "$LOG_DIR/${SAMPLE}_dist.out" 2> "$LOG_DIR/${SAMPLE}_dist.err"
  [[ -f "$RI" ]] || ${NUMA} vg gbwt -p --num-threads "$THREADS" -r "$RI" -Z "$GBZ" \
    > "$LOG_DIR/${SAMPLE}_ri.out" 2> "$LOG_DIR/${SAMPLE}_ri.err"

  run_hap () {
    local READ="$1" TAG="$2"
    [[ -f "$READ" ]] || { echo "[WARN] $TAG reads missing for $SAMPLE: $READ"; return; }
    local OUT_DIR="$OUT_BASE/$SAMPLE/${SAMPLE}_${TAG}"
    mkdir -p "$OUT_DIR"
    local OUT_FA="$OUT_DIR/full.fa"
    echo "[RUN] $SAMPLE $TAG -> $OUT_FA"
    ${NUMA} /usr/bin/time -v \
      "$VG_HAP" -g "$GBZ" -r "$READ" -d "$TMP_DIR" -t "$THREADS" -o "$OUT_FA" \
      > "$LOG_DIR/${SAMPLE}_${TAG}.out" 2> "$LOG_DIR/${SAMPLE}_${TAG}.err"
  }

  run_hap "${READS_DIR}/1x/${SAMPLE}.mhc.1x.fastq"       "1x"
  run_hap "${READS_DIR}/2x/${SAMPLE}.mhc.2x.fastq"       "2x"
  run_hap "${READS_DIR}/Total/${SAMPLE}.mhc.Total.fastq"  "Total"

  rm -rf "$TMP_DIR"

done

echo "Done."
