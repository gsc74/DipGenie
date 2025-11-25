#!/bin/bash
set -euo pipefail

THREADS=128
GRAPH_DIR="Graph"
READS_DIR="Reads"
OUT_BASE="Results_VG"
LOG_DIR="logs_VG"
VG_HAP="./vg_haplotypes.py"
TMP_DIR=tmp_vg_hap


mkdir -p "$TMP_DIR"
mkdir -p "$OUT_BASE" "$LOG_DIR"
shopt -s nullglob

while read -r SAMPLE; do
  [[ -z "$SAMPLE" ]] && continue
  echo "[INFO] $SAMPLE"

  GFA="$GRAPH_DIR/MHC_left_${SAMPLE}.gfa"
  [[ -f "$GFA" ]] || { echo "[WARN] missing $GFA"; continue; }

  # All index files stay in Graph/
  PREFIX="${GRAPH_DIR}/MHC_left_${SAMPLE}"   # e.g. Graph/MHC_left_HG002
  XG="${PREFIX}.xg"
  GBWT="${PREFIX}.gbwt"                       # from gfa2gbwt
  COMBINED_GBWT="${PREFIX}.combined.gbwt"    # from gfa2gbwt
  GBZ="${PREFIX}.gbz"
  DIST="${PREFIX}.dist"
  RI="${PREFIX}.ri"

  # Build GBWT+XG+GBZ if GBZ missing
  if [[ ! -f "$GBZ" ]]; then
    echo "[IDX] ${SAMPLE}: GFA -> GBWT + XG -> GBZ (in Graph/)"
    numactl --cpunodebind=0 --membind=0 gfa2gbwt -b "$PREFIX" < "$GFA" \
      > "$LOG_DIR/${SAMPLE}_gfa2gbwt.out" 2> "$LOG_DIR/${SAMPLE}_gfa2gbwt.err"

    numactl --cpunodebind=0 --membind=0 vg convert -g "$GFA" -x > "$XG"

    numactl --cpunodebind=0 --membind=0 vg gbwt -x "$XG" "$GBWT" --gbz-format -g "$GBZ" \
      > "$LOG_DIR/${SAMPLE}_gbz.out" 2> "$LOG_DIR/${SAMPLE}_gbz.err"
  fi

  # Build DIST and RI if missing (still in Graph/)
  [[ -f "$DIST" ]] || numactl --cpunodebind=0 --membind=0 vg index -j "$DIST" "$GBZ" \
    > "$LOG_DIR/${SAMPLE}_dist.out" 2> "$LOG_DIR/${SAMPLE}_dist.err"

  [[ -f "$RI" ]]   || numactl --cpunodebind=0 --membind=0 vg gbwt -p --num-threads "$THREADS" -r "$RI" -Z "$GBZ" \
    > "$LOG_DIR/${SAMPLE}_ri.out" 2> "$LOG_DIR/${SAMPLE}_ri.err"

  # Haplotype generation
  run_hap () {
    local READ="$1" TAG="$2"
    [[ -f "$READ" ]] || { echo "[WARN] $TAG reads missing for $SAMPLE"; return; }
    local OUT_DIR="$OUT_BASE/$SAMPLE/${SAMPLE}_${TAG}"
    mkdir -p "$OUT_DIR"
    local OUT_FA="$OUT_DIR/full.fa"
    echo "[RUN] $SAMPLE $TAG -> $OUT_FA"
    numactl --cpunodebind=0 --membind=0 /usr/bin/time -v \
      "$VG_HAP" -g "$GBZ" -r "$READ" -d "$TMP_DIR" -t "$THREADS" -o "$OUT_FA" \
      > "$LOG_DIR/${SAMPLE}_${TAG}.out" 2> "$LOG_DIR/${SAMPLE}_${TAG}.err"
  }

  run_hap "$READS_DIR/2x/${SAMPLE}.mhc.2x.fastq"       "2x"
  run_hap "$READS_DIR/4x/${SAMPLE}.mhc.4x.fastq"       "4x"
  run_hap "$READS_DIR/Total/${SAMPLE}.mhc.all.fastq"   "full"

done < "${1:-/dev/stdin}"

