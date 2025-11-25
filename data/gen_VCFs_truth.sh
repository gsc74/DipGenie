set -euo pipefail

# Environment setup
set +u
source ~/.bashrc 2>/dev/null || true
source ~/.bash_profile 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cactus
set -u

# Input paths
HAPS_ROOT_REL="./hprc_haps"
REF_NAME="CHM13.0"
REF_FILE="MHC-CHM13.0.fa"
MAX_CORES=256

SAMPLES=(
  HG002 HG00438 HG00621 HG00741 HG03540
  HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
  HG01891 HG01928 HG01952 HG01978
  HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886
)

# Absolute paths
HAPS_ROOT="$(readlink -f "$HAPS_ROOT_REL")"
REF_ABS="$HAPS_ROOT/$REF_FILE"

mkdir -p Truth
shopt -s nullglob

for S in "${SAMPLES[@]}"; do
  OUTDIR="Truth/${S}"
  mkdir -p "$OUTDIR"

  SEQFILE="${OUTDIR}/MHC_${S}.seqfile"
  H1="$HAPS_ROOT/MHC-${S}.1.fa"
  H2="$HAPS_ROOT/MHC-${S}.2.fa"

  # Seqfile with absolute paths
  cat > "$SEQFILE" <<EOF
${REF_NAME} ${REF_ABS}
${S}.1       ${H1}
${S}.2       ${H2}
EOF

  echo "[BUILD] ${S} -> ${OUTDIR}"

  # Run cactus-pangenome
  rm -rf "${OUTDIR}/js"
  cactus-pangenome "${OUTDIR}/js" "$SEQFILE" \
    --outDir "$OUTDIR" \
    --outName "MHC_${S}" \
    --reference "${REF_NAME}" \
    --vcf \
    --maxCores ${MAX_CORES} --indexCores 32 --mapCores 8 \
    --batchSystem single_machine
  rm -rf "${OUTDIR}/js"
done

echo "Done. Truth VCFs are in ./Truth/<SAMPLE>/"