#!/usr/bin/env bash
# VG+HapCUT2 and PanGenie+HapCUT2 at every depth (one node, depths in turn).
# Usage: bash run_baselines.sh VG|PanGenie [depth ...]   (default: all six depths)
#
# Genotyping: VG = KMC k-mer counts + vg haplotypes --diploid-sampling + call_fasta (the in-process
#             dipcall v0.3 stages of DipGenie: minimap2 v2.31, HTSbox r352);
#             PanGenie = PanGenie v4.2.2 on the uncompressed concatenated SR mates.
# Phasing (HapCUT2 v1.3.4, documented defaults): calls normalised identically
# (bcftools norm -m -any, non-reference GTs, bcftools sort); Hi-C fragments from extractHAIRS --hic 1 on a
# bwa-mem2 -5SPM alignment; SR fragments from extractHAIRS --nf 1 on a bwa-mem2 alignment;
# fragments concatenated and assembled with HAPCUT2 --hic 1. Genotypes are never changed.
# SR-only output (appendix): the same normalised calls phased by HapCUT2 from the SR alignment
# alone with its basic defaults (extractHAIRS --bam --VCF --out; HAPCUT2 --fragments --VCF --output),
# written to $WORK/runs/<method>-HapCUT2-SR-<depth>.
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/05_baselines
method=${1:?VG or PanGenie}; shift
depths=("$@"); [[ $# -gt 0 ]] || depths=(2x 4x 8x 16x 32x Total)
ref=$WORK/download/chm13v2.0.fa.gz
evalref=$WORK/download/reference/CHM13v2.fa
reads=$WORK/reads
set +e  # each run is a subshell with its own set -e; record its exit code and continue
for d in "${depths[@]}"; do
  out=$WORK/runs/$method-HapCUT2-$d; mkdir -p "$out"; rm -f "$out/COMPLETE"
  ( set -Eeuo pipefail
    sr=$reads/SR/$d; hc=$reads/HiC/$d
    timed() { local n=$1; shift; /usr/bin/time -v -o "$out/$n.time" "$@"; }
    if [[ $method == VG ]]; then
      mkdir -p "$out/kmc"; printf '%s\n' "$sr/HG002.R1.fastq.gz" "$sr/HG002.R2.fastq.gz" > "$out/reads.list"
      timed counting "$KMC" -k29 -m16 -okff -t"$THREADS" -hp "@$out/reads.list" "$out/kmc/sample" "$out/kmc" > "$out/kmc.log"
      timed genotyping "$VG" haplotypes --diploid-sampling -v2 -t"$THREADS" -i "$WORK/vg-index/panel.hapl" \
        -k "$out/kmc/sample.kff" -g "$out/sample.gbz" "$WORK/vg-index/panel.gbz" 2> "$out/vg.log"
      "$VG" paths -x "$out/sample.gbz" -F -S recombination > "$out/haplotypes.fa"
      timed calling "$WORK/vg-index/call_fasta" "$out/haplotypes.fa" "$WORK/panel/chr1.fa" "$out/genotypes.vcf.gz" > "$out/calling.log"
      calls=$out/genotypes.vcf.gz
    else
      zcat "$sr/HG002.R1.fastq.gz" "$sr/HG002.R2.fastq.gz" > "$out/reads.fq"
      timed genotyping "$PANGENIE_DIR/PanGenie" -f "$WORK/pangenie-index/panel" -i "$out/reads.fq" -s "$SAMPLE" \
        -j "$THREADS" -t "$THREADS" -e 3000000000 -o "$out/pangenie" > "$out/pangenie.log" 2>&1
      rm -f "$out/reads.fq"
      calls=$out/pangenie_genotyping.vcf
    fi
    # Same normalisation for both callers (PanGenie writes no ##contig line).
    "$BCFTOOLS" annotate --header-lines <(echo '##contig=<ID=chr1,length=248387328>') "$calls" 2> "$out/norm.log" \
      | "$BCFTOOLS" norm -f "$evalref" -m -any 2>> "$out/norm.log" \
      | "$BCFTOOLS" view -i 'GT="alt"' -Ou \
      | "$BCFTOOLS" sort -T "$out/sort-tmp" -Ov -o "$out/calls.unphased.vcf"  # left-alignment can move long deletions out of order
    for lib in HiC SR; do
      bam=$out/${lib,,}.chr1.bam
      [[ -f $bam ]] && continue
      if [[ $lib == HiC ]]; then opt=(-5SPM); in=$hc; else opt=(-R "'@RG\\tID:HG002\\tSM:HG002'"); in=$sr; fi
      st=$SAMTOOLS; [[ $lib == HiC ]] && st=${SAMTOOLS_HIC:-$SAMTOOLS}
      timed align-$lib bash -o pipefail -c "\"$BWA_MEM2\" mem ${opt[*]} -t $THREADS '$ref' '$in/HG002.R1.fastq.gz' '$in/HG002.R2.fastq.gz' 2> '$out/align-$lib.log' \
        | \"$st\" fixmate -m -@8 - - | \"$st\" sort -@32 -m 2G -o '$out/$lib.sorted.bam' -"
      "$st" markdup -@16 "$out/$lib.sorted.bam" "$out/$lib.md.bam" 2> "$out/markdup-$lib.log"
      "$st" index -@16 "$out/$lib.md.bam"
      "$st" view -@16 -b -o "$bam" "$out/$lib.md.bam" chr1 && "$st" index "$bam"
      rm -f "$out/$lib.sorted.bam" "$out/$lib.md.bam" "$out/$lib.md.bam.bai"
    done
    timed fragments-HiC "$HAPCUT2_DIR/extractHAIRS" --hic 1 --bam "$out/hic.chr1.bam" --VCF "$out/calls.unphased.vcf" --out "$out/hic.fragments" > "$out/extracthairs.hic.log" 2>&1
    timed fragments-SR "$HAPCUT2_DIR/extractHAIRS" --nf 1 --bam "$out/sr.chr1.bam" --VCF "$out/calls.unphased.vcf" --out "$out/sr.fragments" > "$out/extracthairs.sr.log" 2>&1
    cat "$out/sr.fragments" "$out/hic.fragments" > "$out/all.fragments"
    timed phasing "$HAPCUT2_DIR/HAPCUT2" --hic 1 --fragments "$out/all.fragments" --VCF "$out/calls.unphased.vcf" \
      --output "$out/hapcut2.haplotype" --outvcf 1 > "$out/hapcut2.log" 2>&1
    "$PYTHON" "$here/hapcut2_vcf.py" "$out/calls.unphased.vcf" "$out/hapcut2.haplotype.phased.VCF" "$out/calls.vcf" > "$out/phase-conversion.json"
    "$BCFTOOLS" view -Oz -o "$out/calls.vcf.gz" "$out/calls.vcf" && "$BCFTOOLS" index -f "$out/calls.vcf.gz" && rm "$out/calls.vcf"
    bash "$here/../06_evaluation/evaluate.sh" "$out/calls.vcf.gz" "$out/evaluation"
    bash "$here/../06_evaluation/phasing.sh" "$out/calls.vcf.gz" "$out/phasing"
    touch "$out/COMPLETE" )
  rc=$?
  echo "[$(date +%F_%T)] $method+HapCUT2 $d exit $rc"
  # SR only: same genotypes, HapCUT2 on SR fragments alone (default settings).
  src=$out; out=$WORK/runs/$method-HapCUT2-SR-$d; rm -rf "$out"; mkdir -p "$out"
  ( set -Eeuo pipefail
    test -f "$src/COMPLETE"
    timed() { local n=$1; shift; /usr/bin/time -v -o "$out/$n.time" "$@"; }
    cp "$src/calls.unphased.vcf" "$out/calls.unphased.vcf"
    timed fragments "$HAPCUT2_DIR/extractHAIRS" --bam "$src/sr.chr1.bam" --VCF "$out/calls.unphased.vcf" --out "$out/sr.fragments" > "$out/extracthairs.sr.log" 2>&1
    timed phasing "$HAPCUT2_DIR/HAPCUT2" --fragments "$out/sr.fragments" --VCF "$out/calls.unphased.vcf" \
      --output "$out/hapcut2.haplotype" --outvcf 1 > "$out/hapcut2.log" 2>&1
    "$PYTHON" "$here/hapcut2_vcf.py" "$out/calls.unphased.vcf" "$out/hapcut2.haplotype.phased.VCF" "$out/calls.vcf" > "$out/phase-conversion.json"
    "$BCFTOOLS" view -Oz -o "$out/calls.vcf.gz" "$out/calls.vcf" && "$BCFTOOLS" index -f "$out/calls.vcf.gz" && rm "$out/calls.vcf"
    bash "$here/../06_evaluation/evaluate.sh" "$out/calls.vcf.gz" "$out/evaluation"
    bash "$here/../06_evaluation/phasing.sh" "$out/calls.vcf.gz" "$out/phasing"
    touch "$out/COMPLETE" )
  rc=$?
  echo "[$(date +%F_%T)] $method+HapCUT2 SR-only $d exit $rc"
done
