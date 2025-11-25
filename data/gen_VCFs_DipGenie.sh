set -euo pipefail

# ---- env (be lenient for system rc files) ----
set +u
source ~/.bashrc 2>/dev/null || true
source ~/.bash_profile 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cactus
set -u

ROOT="Results"
REF_ABS="./hprc_haps/MHC-CHM13.0.fa"   # absolute path
MAX_CORES=256

shopt -s nullglob

# Split full.fa into hap1/hap2 and rewrite headers to >SAMPLE#1 / >SAMPLE#2
split_full() {
  local full="$1" dir="$2" sample="$3"
  local f1="$dir/full_1.fa" f2="$dir/full_2.fa"
  : > "$f1"; : > "$f2"
  awk -v F1="$f1" -v F2="$f2" -v S="$sample" '
    BEGIN{which=2}
    /^>/{
      which = 3 - which;               # toggle 1 <-> 2 per record
      print ">" S "#" which > (which==1?F1:F2)
      next
    }
    { print > (which==1?F1:F2) }
  ' "$full"
}

build_one() {
  local dir="$1" sample="$2" tag="$3"
  local seqfile="$dir/MHC_${sample}_${tag}.seqfile"
  local out_prefix="MHC_${sample}_${tag}"

  # Split + rename headers
  if [[ -f "$dir/full.fa" ]]; then
    split_full "$dir/full.fa" "$dir" "$sample"
  fi
  [[ -s "$dir/full_1.fa" && -s "$dir/full_2.fa" ]] || { echo "[WARN] missing split FASTAs in $dir"; return; }

  # Absolute paths
  local f1_abs f2_abs
  f1_abs="$(readlink -f "$dir/full_1.fa")"
  f2_abs="$(readlink -f "$dir/full_2.fa")"

  # Seqfile (absolute paths)
  cat > "$seqfile" <<EOF
CHM13.0 $REF_ABS
${sample}.1       $f1_abs
${sample}.2       $f2_abs
EOF

  echo "[BUILD] $sample/$tag"
  rm -rf "$dir/js"
  cactus-pangenome "$dir/js" "$seqfile" \
    --outDir "$dir" \
    --outName "$out_prefix" \
    --reference CHM13.0 \
    --vcf \
    --maxCores $MAX_CORES --indexCores 32 --mapCores 8 \
    --batchSystem single_machine
  rm -rf "$dir/js"
}

# Traverse Results/<SAMPLE>/<SAMPLE>_{2x,4x,full}/
for d in "$ROOT"/*/*/; do
  base="$(basename "$d")"
  sample="$(basename "$(dirname "$d")")"
  case "$base" in
    "${sample}_2x"| "${sample}_4x"| "${sample}_full")
      tag="${base#${sample}_}"
      build_one "$d" "$sample" "$tag"
      ;;
  esac
done

echo "Done."
