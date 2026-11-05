#!/usr/bin/env bash
# Leave-one-out chr1 panel, DipGenie graph index and chr1 truth subsets.
#   1. remove HG002 from the HPRC v2 bubble VCF (231 -> 230 samples);
#   2. extract chr1 and recompute AC/AN/AF from the retained genotypes;
#   3. drop ALT alleles no retained haplotype carries, and records left without an ALT;
#   4. audit sample names and ALT support; build the DipGenie index; subset the truth to chr1.
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/02_panel
export PATH="$(dirname "$(command -v "$BCFTOOLS")"):$PATH"
dl=$WORK/download
out=$WORK/panel
mkdir "$out"
"$PYTHON" "$here/panel_filter_test.py" "$out/toy-filter"

"$SAMTOOLS" faidx "$dl/reference/CHM13v2.fa" chr1 > "$out/chr1.fa"
"$SAMTOOLS" faidx "$out/chr1.fa"

[[ -f $dl/mc_filtered_ids.vcf.gz.tbi || -f $dl/mc_filtered_ids.vcf.gz.csi ]] || "$BCFTOOLS" index --threads 8 "$dl/mc_filtered_ids.vcf.gz"
"$BCFTOOLS" view --no-version --threads 8 -s ^"$SAMPLE" -r chr1 "$dl/mc_filtered_ids.vcf.gz" -Oz -o "$out/chr1.raw.vcf.gz"
"$BCFTOOLS" reheader -f "$out/chr1.fa.fai" -o "$out/chr1.raw.header.vcf.gz" "$out/chr1.raw.vcf.gz"
mv "$out/chr1.raw.header.vcf.gz" "$out/chr1.raw.vcf.gz"
"$BCFTOOLS" index --threads 8 "$out/chr1.raw.vcf.gz"
"$PYTHON" "$here/audit_panel.py" samples "$out/chr1.raw.vcf.gz" "$out"
bash "$here/filter_panel.sh" "$out/chr1.raw.vcf.gz" "$out/chr1.noHG002.vcf.gz" "$out/keep.samples.txt"
"$PYTHON" "$here/audit_panel.py" verify "$out/chr1.noHG002.vcf.gz" "$out" > "$out/validation.log"
"$BCFTOOLS" index --tbi "$out/chr1.noHG002.vcf.gz"

"$DIPGENIE" --vcf "$out/chr1.noHG002.vcf.gz" --reference "$out/chr1.fa" \
  --write-index "$out/chr1.gindex" -t "$THREADS" > "$out/index.stdout" 2> "$out/index.log"

mkdir "$WORK/truth"
for type in smvar stvar; do
  [[ -f $dl/HG002_CHM13v2.0_v5.0q_$type.vcf.gz.tbi ]] || "$BCFTOOLS" index -t "$dl/HG002_CHM13v2.0_v5.0q_$type.vcf.gz"
  "$BCFTOOLS" view -r chr1 "$dl/HG002_CHM13v2.0_v5.0q_$type.vcf.gz" -Oz -o "$WORK/truth/$type.vcf.gz"
  "$BCFTOOLS" index "$WORK/truth/$type.vcf.gz"
  awk '$1=="chr1"' "$dl/HG002_CHM13v2.0_v5.0q_$type.benchmark.bed" > "$WORK/truth/$type.bed"
done
sha256sum "$out/chr1.noHG002.vcf.gz" "$out/chr1.fa" "$WORK/truth"/*.vcf.gz "$WORK/truth"/*.bed > "$out/outputs.sha256"
touch "$out/COMPLETE"
