#!/usr/bin/env bash
# Haplotype-level phasing as for haplotype-resolved assemblies (hifiasm; Cheng et al. 2021, 2022):
# yak trioeval switch and Hamming error rates against the parental 31-mers of the HG002 Q100 v1.1
# chr1 (maternal and paternal sequences, so no parental reads are needed).
#   actual:   haplotypes = bcftools consensus -H 1/-H 2 of the call set on CHM13 chr1;
#   roofline: oracle phase: heterozygous calls matching a truth heterozygote take the truth phase,
#             every other heterozygous call is flipped with the preceding matched call, and no
#             genotype changes (oracle_phase.py); the remaining error is due to genotypes.
# Bases outside the GIAB v5.0q benchmark regions are masked to N (--mask), and yak skips every
# 31-mer that contains N, so the truth VCF itself scores 0/0.
# Usage: phasing.sh calls.vcf.gz output_dir
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/06_evaluation
calls=$(readlink -f "${1:?calls.vcf.gz}"); out=${2:?output_dir}
t=$WORK/truth/yak
# Keep a shared lock for the entire scoring run; publication/recovery requires
# exclusive ownership. Completed caches are reused, never replaced under readers.
source "$here/locked_cache.sh"
build_yak_cache() {
  local tmp=$1
  for p in MATERNAL PATERNAL; do
    zcat "$WORK/download/hg002v1.1.fasta.gz" | awk -v n=">chr1_$p" '/^>/{k=($1==n)}k' > "$tmp/$p.fa"
  done
  "$YAK" count -t "$THREADS" -o "$tmp/mat.yak" "$tmp/MATERNAL.fa"
  "$YAK" count -t "$THREADS" -o "$tmp/pat.yak" "$tmp/PATERNAL.fa"
  "$SAMTOOLS" faidx "$WORK/panel/chr1.fa" chr1 > "$tmp/chm13.chr1.fa"
  "$SAMTOOLS" faidx "$tmp/chm13.chr1.fa"
  # Mask = chr1 minus the small-variant benchmark regions.
  sort -k1,1 -k2,2n "$WORK/truth/smvar.bed" | awk -v L="$(cut -f2 "$tmp/chm13.chr1.fa.fai")" \
    'BEGIN{OFS="\t";e=0} {if($2>e)print "chr1",e,$2; if($3>e)e=$3} END{if(e<L)print "chr1",e,L}' > "$tmp/mask.bed"
  "$BCFTOOLS" concat -a "$WORK/truth/smvar.vcf.gz" "$WORK/truth/stvar.vcf.gz" -Ou \
    | "$BCFTOOLS" norm -f "$tmp/chm13.chr1.fa" -m -any -Ou | "$BCFTOOLS" view -i 'GT!="mis"' -Ou \
    | "$BCFTOOLS" sort -Oz -o "$tmp/truth.norm.vcf.gz"
  "$BCFTOOLS" index -t "$tmp/truth.norm.vcf.gz"
  "$YAK" version > "$tmp/yak.version"
}
acquire_read_cache "$t" build_yak_cache
mkdir -p "$out"
score() { # label vcf
  for h in 1 2; do
    "$BCFTOOLS" consensus -s "$SAMPLE" -H $h --mask "$t/mask.bed" -f "$t/chm13.chr1.fa" "$2" 2>> "$out/consensus.log" \
      | sed "s/^>.*/>hap$h/"
  done > "$out/$1.haplotypes.fa"
  # Both haplotypes in one file: yak sums the counts of the two sequences.
  "$YAK" trioeval -c1 -d1 -t "$THREADS" "$t/pat.yak" "$t/mat.yak" "$out/$1.haplotypes.fa" > "$out/$1.trioeval.txt"
  rm -f "$out/$1.haplotypes.fa"
}
score actual "$calls"
"$BCFTOOLS" norm -f "$t/chm13.chr1.fa" -m -any "$calls" -Ou | "$BCFTOOLS" view -i 'GT!="mis"' -Ou \
  | "$BCFTOOLS" sort -Oz -o "$out/calls.norm.vcf.gz"
"$BCFTOOLS" index -f -t "$out/calls.norm.vcf.gz"
"$PYTHON" "$here/oracle_phase.py" "$t/truth.norm.vcf.gz" "$out/calls.norm.vcf.gz" "$out/roofline.vcf" 2> "$out/roofline.log"
"$BCFTOOLS" view -Oz -o "$out/roofline.vcf.gz" "$out/roofline.vcf" && rm -f "$out/roofline.vcf"
"$BCFTOOLS" index -f -t "$out/roofline.vcf.gz"
score roofline "$out/roofline.vcf.gz"
# W = switch error rate, H = Hamming error rate (yak output lines).
for k in actual roofline; do awk -v k=$k '$1=="W"||$1=="H"{printf "%s\t%s\t%.4f\n", k, ($1=="W"?"switch":"hamming"), 100*$4}' "$out/$k.trioeval.txt"; done > "$out/phasing.tsv"
cat "$out/phasing.tsv"
