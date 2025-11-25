set -euo pipefail

# Activate whatshap env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate whatshap

REF="../BKP/hprc_haps/MHC-CHM13.0.fa"
SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540 HG01106 HG01109 HG01123 HG01258 HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(2x 4x full)

mkdir -p Evaluation_VG
shopt -s nullglob

ensure_tbi () {
  local v="$1"
  [[ -f "${v}.tbi" ]] || tabix -p vcf "$v"
}

for S in "${SAMPLES[@]}"; do
  TRUTH="Truth/${S}/MHC_${S}.vcf.gz"
  if [[ ! -f "$TRUTH" ]]; then
    echo "[WARN] Missing truth VCF for $S: $TRUTH"; continue
  fi
  ensure_tbi "$TRUTH"

  for T in "${TAGS[@]}"; do
    TEST="Results_VG/${S}/${S}_${T}/MHC_${S}_${T}.vcf.gz"
    if [[ ! -f "$TEST" ]]; then
      echo "[WARN] Missing test VCF for $S $T: $TEST"; continue
    fi
    ensure_tbi "$TEST"

    OUTDIR="Evaluation/${S}/${S}_${T}"
    mkdir -p "$OUTDIR"

    echo "[RUN] $S $T"
    whatshap compare "$TRUTH" "$TEST" > "${OUTDIR}/SER.txt"
  done
done

echo "Done. See Evaluation/<SAMPLE>/<SAMPLE_TAG>/SER.txt"