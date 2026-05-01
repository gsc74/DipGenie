#!/bin/bash
set -euo pipefail

# Evaluate VG results: whatshap compare (SER + mismatch rate)
# Requires: conda activate whatshap

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate whatshap

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)

TRUTH_DIR="../Truth"
TEST_DIR="Results"

mkdir -p Evaluation
shopt -s nullglob

ensure_tbi () {
  local v="$1"
  [[ -f "${v}.tbi" ]] || tabix -p vcf "$v"
}

for S in "${SAMPLES[@]}"; do
  TRUTH="${TRUTH_DIR}/${S}/MHC_${S}.vcf.gz"
  [[ -f "$TRUTH" ]] || { echo "[WARN] Missing truth VCF: $TRUTH"; continue; }
  ensure_tbi "$TRUTH"

  for T in "${TAGS[@]}"; do
    TEST="${TEST_DIR}/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    [[ -f "$TEST" ]] || { echo "[WARN] Missing test VCF: $TEST"; continue; }
    ensure_tbi "$TEST"

    OUTDIR="Evaluation/${S}/${S}_${T}"
    mkdir -p "$OUTDIR"

    echo "[RUN] whatshap compare $S $T"
    whatshap compare "$TRUTH" "$TEST" > "${OUTDIR}/SER.txt"
  done
done

echo "Done. Results in Evaluation/<SAMPLE>/<SAMPLE>_<TAG>/SER.txt"
