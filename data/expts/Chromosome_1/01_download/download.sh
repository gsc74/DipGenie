#!/usr/bin/env bash
# Download every external input of the chromosome 1 benchmark into $WORK/download and
# verify sizes and published checksums (sources.tsv). Completed files are reused.
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/01_download
out=$WORK/download
mkdir -p "$out"
cd "$out"

verify() {
  local name=$1 bytes=$2 checksum=$3
  [[ $(stat -c %s "$name") == "$bytes" ]] || { echo "size mismatch: $name" >&2; return 1; }
  case $checksum in
    sha256:*) echo "${checksum#sha256:}  $name" | sha256sum --check --quiet ;;
    md5:*) echo "${checksum#md5:}  $name" | md5sum --check --quiet ;;
  esac
}

# 1) Reference, panel VCF, GIAB truth, stratifications and Hi-C reads (HTTPS, no credentials).
tail -n +2 "$here/sources.tsv" | while IFS=$'\t' read -r name url bytes checksum; do
  if [[ -f $name ]] && verify "$name" "$bytes" "$checksum"; then echo "reuse $name"; continue; fi
  curl -fL --retry 10 --retry-delay 30 -C - -o "$name.partial" "$url"
  mv "$name.partial" "$name"
  verify "$name" "$bytes" "$checksum"
  echo "downloaded $name"
done

# 2) Evaluation reference (CHM13v2.fa with RTG SDF) and GIAB v3.6 stratifications.
[[ -f reference/CHM13v2.fa ]] || { mkdir -p reference; tar -xzf CHM13v2.tar.gz -C reference --strip-components=1; }
[[ -f reference/CHM13v2.fa.fai ]] || "$SAMTOOLS" faidx reference/CHM13v2.fa
[[ -d reference/CHM13v2.sdf ]] || "$RTG" format -o reference/CHM13v2.sdf reference/CHM13v2.fa
[[ -d stratifications ]] || { mkdir stratifications; tar -xzf 'genome-stratifications-CHM13@all.tar.gz' -C stratifications --strip-components=1; }
[[ -f chm13v2.0.fa.gz.fai ]] || "$SAMTOOLS" faidx chm13v2.0.fa.gz
[[ -f chm13v2.0.fa.gz.bwt.2bit.64 ]] || "$BWA_MEM2" index chm13v2.0.fa.gz

# 3) Paired short reads, SRA SRR11321732 (HG002, Illumina 2x250 bp).
acc=SRR11321732
if [[ ! -f SR/COMPLETE ]]; then
  mkdir -p SR
  "$PREFETCH" --transport http --verify yes --max-size 500G -O SR "$acc"
  "$FASTERQ_DUMP" --split-files --threads "$THREADS" --temp SR/tmp -O SR "SR/$acc"
  for mate in 1 2; do
    "$PIGZ" -p "$THREADS" -1 "SR/${acc}_$mate.fastq"
    "$PIGZ" -t "SR/${acc}_$mate.fastq.gz"
  done
  rm -rf SR/tmp "SR/$acc"
  touch SR/COMPLETE
fi
echo "PASS all downloads in $out"
