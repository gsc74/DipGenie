#!/usr/bin/env bash
# Nested chr1 read sets for SR and Hi-C at 2x, 4x, 8x, 16x, 32x and Total.
#
# SR (SRR11321732): bwa-mem2 to the full CHM13 v2.0 reference; keep pairs whose two primary
# alignments lie entirely in chr1 (paired_region.py extract); rank pairs by a seeded blake2b
# hash of the read ID and take prefixes until the aligned M/=/X bases reach depth x chr1 length
# (prepare_chr1_sr.py). Every smaller set is a subset of every larger one.
#
# Hi-C: bwa-mem2 -5SP -Y to the full reference; keep pairs whose primary alignments lie in chr1
# (regions.py extract), restore all their alignments, build contacts with pairtools
# (parse --min-mapq 30 --walks-policy 5unique, sort, dedup --max-mismatch 3) and draw nested
# read and contact subsets with the same rule (regions.py sample).
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/03_reads
dl=$WORK/download
ref=$dl/chm13v2.0.fa.gz
out=$WORK/reads
mkdir -p "$out"
export PAIRTOOLS
(cd "$here" && "$PYTHON" -m unittest test_paired_region test_regions test_prepare_chr1_sr test_audit_inputs) > "$out/toy-tests.log" 2>&1

# SR
src=$out/SR-source
if [[ ! -f $src/PASS ]]; then
  mkdir "$src"
  align=("$BWA_MEM2" mem -K 100000000 -t "$THREADS" -R '@RG\tID:HG002.SRR11321732\tSM:HG002\tPL:ILLUMINA'
         "$ref" "$dl/SR/SRR11321732_1.fastq.gz" "$dl/SR/SRR11321732_2.fastq.gz")
  extract=("$PYTHON" "$here/paired_region.py" extract --sam - --region chr1:1-248387328
           --output "$src/extracted" --accession SRR11321732)
  printf '%s | %s\n' "${align[*]}" "${extract[*]}" > "$src/extraction.command"
  "${align[@]}" 2> "$src/alignment.log" | "${extract[@]}" 2> "$src/extraction.log"
  touch "$src/PASS"
fi
[[ -f $out/SR/PASS.json ]] || "$PYTHON" "$here/prepare_chr1_sr.py" --source "$src" --output "$out/SR" --seed 42
"$PYTHON" "$here/audit_inputs.py" "$out/SR" chr1 --output "$out/SR-audit.json"

# Hi-C
hic=$out/HiC-source
mkdir -p "$hic"
r1=$dl/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fastq.gz
r2=$dl/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fastq.gz
if [[ ! -f $hic/alignment.PASS ]]; then
  "$BWA_MEM2" mem -5SP -Y -K 100000000 -t "$THREADS" \
    -R '@RG\tID:HG002.HiC_2_NovaSeq_rep1_run2_S1_L001\tSM:HG002\tPL:ILLUMINA\tLB:HiC_rep1_run2' \
    "$ref" "$r1" "$r2" 2> "$hic/alignment.log" | "$SAMTOOLS" view -@ 8 -b -o "$hic/full.namegrouped.bam" -
  touch "$hic/alignment.PASS"
fi
if [[ ! -f $hic/extraction.PASS ]]; then
  rm -rf "$hic/extracted"
  printf 'chr1\t0\t248387328\n' > "$hic/chr1.bed"
  "$SAMTOOLS" view -@ 8 -h -L "$hic/chr1.bed" "$hic/full.namegrouped.bam" |
    "$PYTHON" "$here/regions.py" extract --bam - --output "$hic/extracted" 2> "$hic/extraction.log"
  # All alignments of accepted pairs, including supplementary ones, for split-read parsing.
  "$SAMTOOLS" view -@ 8 -b -N "$hic/extracted/chr1/read_ids.txt" \
    -o "$hic/extracted/chr1/all-alignments.bam" "$hic/full.namegrouped.bam"
  touch "$hic/extraction.PASS"
fi
d=$hic/extracted/chr1
if [[ ! -f $d/contacts.PASS ]]; then
  "$SAMTOOLS" faidx "$ref"
  cut -f1,2 "$ref.fai" > "$hic/chrom.sizes"
  "$PAIRTOOLS" parse -c "$hic/chrom.sizes" --assembly CHM13v2.0 --min-mapq 30 --walks-policy 5unique \
    --drop-sam --add-columns mapq,read_side,algn_idx,cigar --output-stats "$d/parse.stats" \
    -o "$d/parsed.pairs.gz" "$d/all-alignments.bam"
  "$PAIRTOOLS" sort --nproc 16 --memory 32G --tmpdir "$d" -o "$d/sorted.pairs.gz" "$d/parsed.pairs.gz"
  "$PAIRTOOLS" dedup --max-mismatch 3 --backend scipy --output "$d/nodups.pairs.gz" \
    --output-dups "$d/duplicates.pairs.gz" --output-unmapped "$d/unmapped.pairs.gz" \
    --output-stats "$d/dedup.stats" "$d/sorted.pairs.gz"
  touch "$d/contacts.PASS"
fi
[[ -f $out/HiC/PASS.json ]] || "$PYTHON" "$here/regions.py" sample --input "$d" --output "$out/HiC" --region chr1 --seed 42
"$PYTHON" "$here/audit_inputs.py" "$out/HiC" chr1 --output "$out/HiC-audit.json"

# Manifest for 04_dipgenie/run.sh
{ printf 'depth\tsr_r1\tsr_r2\thic_r1\thic_r2\tcontacts\n'
  for x in 2x 4x 8x 16x 32x Total; do
    printf '%s\tSR/%s/HG002.R1.fastq.gz\tSR/%s/HG002.R2.fastq.gz\tHiC/%s/HG002.R1.fastq.gz\tHiC/%s/HG002.R2.fastq.gz\tHiC/%s/contacts.pairs.gz\n' \
      "$x" "$x" "$x" "$x" "$x" "$x"
  done; } > "$out/manifest.tsv"
cat "$out/SR/coverage.csv" "$out/HiC/coverage.csv"
touch "$out/COMPLETE"
