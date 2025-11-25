set -euo pipefail

# Activate tools (bcftools/tabix/truvari must be in this env)
set +u
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate truvari || true
set -u

REF="./hprc_haps/MHC-CHM13.0.fa"
OUTROOT="SV_Evaluation_VG"
SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540 HG01106 HG01109 HG01123 HG01258 HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
DEPTHS=(2x 4x full)

mkdir -p "$OUTROOT"
shopt -s nullglob

for S in "${SAMPLES[@]}"; do
  for D in "${DEPTHS[@]}"; do
    BASE_VCF="Results_VG/${S}/${S}_${D}/MHC_${S}_${D}.vcf.gz"
    TRUTH_VCF="Truth/${S}/MHC_${S}.vcf.gz"
    OUTDIR="${OUTROOT}/${S}/${S}_${D}"
    TMPD="${OUTDIR}/tmp"
    mkdir -p "$OUTDIR" "$TMPD"

    if [[ ! -f "$BASE_VCF" ]]; then
      echo "[WARN] Missing base VCF: $BASE_VCF"; continue
    fi
    if [[ ! -f "$TRUTH_VCF" ]]; then
      echo "[WARN] Missing truth VCF: $TRUTH_VCF"; continue
    fi

    # Normalize, split multiallelics, THEN sort before indexing
    BASE_SPLIT="${OUTDIR}/MHC_${S}_${D}.base.split.vcf.gz"
    TRUTH_SPLIT="${OUTDIR}/MHC_${S}.truth.split.vcf.gz"

    echo "[NORM] $S $D (base)"
    bcftools norm -m -both -f "$REF" "$BASE_VCF" \
      | bcftools sort -T "$TMPD" -Oz -o "$BASE_SPLIT"
    echo "[NORM] $S $D (truth)"
    bcftools norm -m -both -f "$REF" "$TRUTH_VCF" \
      | bcftools sort -T "$TMPD" -Oz -o "$TRUTH_SPLIT"

    tabix -f -p vcf "$BASE_SPLIT"
    tabix -f -p vcf "$TRUTH_SPLIT"

    # Truvari benchmark (Precision/Recall/F1 in summary.txt)
    BENCH_DIR="${OUTDIR}/bench"
    rm -rf "$BENCH_DIR"
    echo "[TRUVARI] $S $D -> $BENCH_DIR"
    truvari bench \
      -f "$REF" \
      -b "$BASE_SPLIT" \
      -c "$TRUTH_SPLIT" \
      -o "$BENCH_DIR"

    [[ -f "${BENCH_DIR}/summary.txt" ]] && cp "${BENCH_DIR}/summary.txt" "${OUTDIR}/summary.metrics.txt"

    rm -rf "$TMPD"
  done
done

echo "Done. See ${OUTROOT}/<SAMPLE>/<SAMPLE_DEPTH>/{*.split.vcf.gz, bench/, summary.metrics.txt}"