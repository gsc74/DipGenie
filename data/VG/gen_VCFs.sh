#!/bin/bash
set -euo pipefail

# Generate test VCFs from VG output using cactus-pangenome
# Requires: conda activate cactus

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cactus

ROOT="Results"
REF_ABS="$(readlink -f ../hprc_haps/MHC-CHM13.0.fa)"
MAX_CORES=128

shopt -s nullglob

split_full() {
  local full="$1" dir="$2" sample="$3"
  local f1="$dir/full_1.fa" f2="$dir/full_2.fa"
  : > "$f1"; : > "$f2"
  awk -v F1="$f1" -v F2="$f2" -v S="$sample" '
    BEGIN{which=2}
    /^>/{
      which = 3 - which;
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

  if [[ -f "$dir/full.fa" ]]; then
    split_full "$dir/full.fa" "$dir" "$sample"
  fi
  [[ -s "$dir/full_1.fa" && -s "$dir/full_2.fa" ]] || { echo "[WARN] missing split FASTAs in $dir"; return; }

  local f1_abs f2_abs
  f1_abs="$(readlink -f "$dir/full_1.fa")"
  f2_abs="$(readlink -f "$dir/full_2.fa")"

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

for d in "$ROOT"/*/*/; do
  base="$(basename "$d")"
  sample="$(basename "$(dirname "$d")")"
  case "$base" in
    "${sample}_1x"|"${sample}_2x"|"${sample}_Total")
      tag="${base#${sample}_}"
      build_one "$d" "$sample" "$tag"
      ;;
  esac
done

echo "Done."
