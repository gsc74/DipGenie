# Chromosome 1 benchmark (HG002, HPRC v2 panel, SR + Hi-C)

Scripts for chromosome 1 experiments: DipGenie v2,
VG + HapCUT2 and PanGenie + HapCUT2 on nested short-read (SR) and Hi-C read sets
of HG002 at 2x, 4x, 8x, 16x, 32x and the full chr1 data (Total), all evaluated
against GIAB v5.0q on CHM13 v2.0.

The launcher runs DipGenie v2 with its default settings.

Stages run in order; each one writes under `$WORK` and refuses to overwrite
completed outputs.

| Stage | Script | Output in `$WORK` |
|---|---|---|
| 1. Download | [01_download/download.sh](01_download/download.sh) | `download/` |
| 2. Panel | [02_panel/prepare_panel.sh](02_panel/prepare_panel.sh) | `panel/`, `truth/` |
| 3. Reads | [03_reads/prepare_reads.sh](03_reads/prepare_reads.sh) | `reads/{SR,HiC}/<depth>/`, `reads/manifest.tsv` |
| 4. DipGenie v2 | [04_dipgenie/run.sh](04_dipgenie/run.sh) | `runs/<LABEL>/<mixed,SR,HiC>/<depth>/`; default label `DipGenie` |
| 5. Baselines | [05_baselines/prepare_indexes.sh](05_baselines/prepare_indexes.sh), [05_baselines/run_baselines.sh](05_baselines/run_baselines.sh) | `vg-index/`, `pangenie-index/`, `runs/<tool>-HapCUT2-<depth>/` (SR + Hi-C), `runs/<tool>-HapCUT2-SR-<depth>/` (SR only) |
| 6. Evaluation | [06_evaluation/evaluate.sh](06_evaluation/evaluate.sh), [06_evaluation/phasing.sh](06_evaluation/phasing.sh) | `<run>/evaluation/`, `<run>/phasing/` |
| 7. SNP error analysis | [07_snp_errors/snp_errors.sh](07_snp_errors/snp_errors.sh) | `snp-errors/<mode>.tsv` |

## Setup

Build DipGenie (`make` in the repository root), then:

```bash
cp data/expts/Chromosome_1/config.example.sh /abs/path/config.sh   # edit REPO, WORK and tool paths
export CONFIG=/abs/path/config.sh
```

`config.example.sh` lists every external tool with the version used in the
benchmark. Python scripts need `pysam`; evaluation additionally needs hap.py,
RTG Tools, Truvari and yak. DipGenie and the VG `call_fasta` helper use the
same in-process assembly caller (dipcall v0.3 stages with minimap2 v2.31 and
HTSbox r352); neither invokes `run-dipcall` or bcftools. The comparison/evaluation pipelines still require their listed
tools. All scripts run with plain `bash`.

```bash
E=data/expts/Chromosome_1
bash $E/01_download/download.sh
bash $E/02_panel/prepare_panel.sh
bash $E/03_reads/prepare_reads.sh
for i in 0 1 2 3 4 5; do bash $E/04_dipgenie/run.sh $i; done           # 2x 4x 8x 16x 32x Total; SR + Hi-C
for i in 0 1 2 3 4 5; do MODE=SR bash $E/04_dipgenie/run.sh $i; done   # short reads only
bash $E/05_baselines/prepare_indexes.sh
bash $E/05_baselines/run_baselines.sh VG          # optional depth list, default all six
bash $E/05_baselines/run_baselines.sh PanGenie
```

## 1. Inputs

[01_download/sources.tsv](01_download/sources.tsv) lists every URL with its size
and checksum; `download.sh` verifies both.

- **Pangenome panel**: HPRC release 2 Minigraph-Cactus bubble VCF,
  [mc_filtered_ids.vcf.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/scratch/2026_03_30_pangenie/mc_filtered_ids.vcf.gz)
  (231 samples, CHM13 coordinates).
- **References**: the T2T CHM13 v2.0 analysis set
  ([chm13v2.0.fa.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz))
  for read alignment, and the GIAB CHM13v2 bundle (`CHM13v2.tar.gz`)
  for evaluation. Their chr1 sequences are identical.
- **Truth**: GIAB HG002 v5.0q small-variant and structural-variant VCFs with
  benchmark BEDs, GIAB v3.6 CHM13 stratifications, and the HG002 Q100 v1.1
  diploid assembly ([hg002v1.1.fasta.gz](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz)),
  whose maternal and paternal chr1 give the parental k-mers for phasing.
- **Short reads**: SRA [SRR11321732](https://www.ncbi.nlm.nih.gov/sra/SRR11321732)
  (HG002, Illumina 2x250 bp), via `prefetch` + `fasterq-dump --split-files`.
- **Hi-C**: HPP HG002 `HiC_2_NovaSeq_rep1_run2_S1_L001` R1/R2.

## 2. Leave-one-out chr1 panel

`prepare_panel.sh`:

1. removes the target sample, `bcftools view -s ^HG002` (231 -> 230 samples);
2. keeps chr1 (`-r chr1`), recomputes `AC/AN/AF` from the retained genotypes
   (`+fill-tags`), trims ALT alleles carried by no retained haplotype
   (`view -a`), and drops records left without an ALT (`-c 1:nref`)
   ([filter_panel.sh](02_panel/filter_panel.sh));
3. audits sample names and ALT support ([audit_panel.py](02_panel/audit_panel.py));
4. builds the DipGenie graph index (230 samples + the CHM13 reference path) and
   restricts the truth VCFs/BEDs to chr1.

The resulting panel `chr1.noHG002.vcf.gz` has 2,899,869 records and is the only
panel used by DipGenie, VG and PanGenie.

## 3. Nested read sets

Depth is defined by aligned bases on chr1 (CHM13 chr1 length 248,387,328 bp).

- **SR**: bwa-mem2 against the full reference; pairs whose two primary
  alignments lie entirely in chr1 are kept ([paired_region.py](03_reads/paired_region.py)).
  Pairs are ranked by a seeded (`--seed 42`) blake2b hash of the read name and
  prefixes are taken until the aligned bases reach `depth x chr1 length`
  ([prepare_chr1_sr.py](03_reads/prepare_chr1_sr.py)).
- **Hi-C**: bwa-mem2 `-5SP -Y`; chr1 pairs are extracted and all their alignments
  restored ([regions.py](03_reads/regions.py) `extract`), then pairtools
  `parse --min-mapq 30 --walks-policy 5unique`, `sort`, `dedup --max-mismatch 3`;
  duplicate pairs are removed. Reads are subsampled with the same hash rule
  (`regions.py sample`).

Every smaller set is a subset of every larger one, and all methods read the same
files ([audit_inputs.py](03_reads/audit_inputs.py) checks this).
`reads/manifest.tsv` lists the files per depth.

These alignments only select the benchmark reads and define depth; DipGenie
uses the paired FASTQs alone and localizes the reads on the graph itself.

## 4. DipGenie v2

`run.sh` runs [run.py](04_dipgenie/run.py) at one depth with SR + Hi-C,
64 threads bound to one NUMA node and recombination budget `-R 2140`; it
records input/binary hashes and `/usr/bin/time -v` resources, and then evaluates.
`MODE=SR` or `MODE=HiC` selects one library. The launcher refuses an allocation
with fewer than 64 physical cores on one NUMA node, rather than silently using
SMT or another NUMA node.

Give every set of runs a fresh `LABEL`; existing outputs are never overwritten. `inputs.json` records the exact command, binary
hash, input hashes and allocated CPUs. Component `stages.tsv` records timing
boundaries; `native-observations.tsv.gz` and the factor files record read
placement and linkage weights.

## 5. Baselines

[prepare_indexes.sh](05_baselines/prepare_indexes.sh) builds, from the same panel:

- VG 1.78.0: `construct`, GBWT with the reference path, GBZ, distance and
  r-index, `vg haplotypes` index; `make build/call_fasta` links the same complete
  caller as DipGenie, and its binary hash/source revision are saved with the
  indexes.
- PanGenie v4.2.2: `PanGenie-index`.

[run_baselines.sh](05_baselines/run_baselines.sh) genotypes and phases each depth:

- **VG**: KMC k=29 counts of the SR reads, `vg haplotypes --diploid-sampling`,
  then the same caller (`call_fasta`) on the two sampled haplotypes.
- **PanGenie**: `PanGenie -s HG002` on the SR reads (its own genotyping VCF).
- **HapCUT2 v1.3.4** at its documented settings, identical for both tools:
  - calls are split with `bcftools norm -m -any`, non-reference GTs are kept, and
    records are re-sorted (`bcftools sort`; left-alignment can move long deletions);
  - Hi-C read pairs are aligned with bwa-mem2 `-5SPM`, SR pairs with bwa-mem2;
  - fragments: `extractHAIRS --hic 1` (Hi-C) and `extractHAIRS --nf 1` (SR),
    concatenated;
  - `HAPCUT2 --hic 1`, converted to a phased VCF with
    [hapcut2_vcf.py](05_baselines/hapcut2_vcf.py). Genotypes are never changed.
- **SR only** (same script, written to `runs/<tool>-HapCUT2-SR-<depth>/`): the same
  normalised calls, phased by HapCUT2 from the SR alignment alone with its basic
  defaults (`extractHAIRS --bam --VCF --out`, `HAPCUT2 --fragments --VCF --output`).
  Genotypes, and therefore SNP/SV accuracy, equal the SR + Hi-C run; only phasing differs.

## 6. Evaluation

[evaluate.sh](06_evaluation/evaluate.sh) applies one protocol to every call set:

- **SNPs**: hap.py against v5.0q small variants inside the benchmark BED (alleles
  longer than 10 kb are excluded for every method because hap.py cannot
  preprocess them).
- **SVs (>= 50 bp)**: Truvari 5.5.0 against v5.0q structural variants.
- **Phasing** ([phasing.sh](06_evaluation/phasing.sh)), as for haplotype-resolved
  assemblies (hifiasm): the two haplotypes of each call set are built with
  `bcftools consensus -H 1/-H 2` on CHM13 chr1, with bases outside the
  small-variant benchmark regions masked to N, and `yak trioeval -c1 -d1`
  classifies their 31-mers against the maternal and paternal chr1 of the HG002
  Q100 v1.1 assembly (`yak count` on each; every k-mer occurs once, hence
  `-c1 -d1`). Both haplotypes are evaluated in one file, so yak sums their
  counts. Its `W` line is the switch error rate and its `H` line the Hamming
  error rate, `sum_i min(p_i, m_i) / sum_i (p_i + m_i)` over the paternal and
  maternal 31-mers `p_i`, `m_i` of haplotype `i`. Masked 31-mers are skipped, so
  the truth VCF scores 0/0.
- **Truth-guided rephasing** (diagnostic, rows `roofline` in `phasing.tsv`): the same
  evaluation after [oracle_phase.py](06_evaluation/oracle_phase.py) gives every heterozygous
  call that matches a truth heterozygote (same CHROM/POS/REF/ALT after
  `bcftools norm -m -any`) the truth orientation and flips every other heterozygous call
  exactly when the preceding matched call was flipped; genotypes are never changed. It
  shows how much of the error a better phase assignment of the same calls could remove.
  It also repairs switches within blocks and phases unmatched calls heuristically, so it
  is neither a lower bound nor an exact split into genotype and orientation error.

```bash
bash $E/06_evaluation/evaluate.sh $WORK/runs/DipGenie/mixed/8x/calls.vcf.gz $WORK/runs/DipGenie/mixed/8x/evaluation
bash $E/06_evaluation/phasing.sh  $WORK/runs/DipGenie/mixed/8x/calls.vcf.gz $WORK/runs/DipGenie/mixed/8x/phasing
```

[summarize.py](06_evaluation/summarize.py) prints every column of the benchmark
results table for one depth as TSV, read from the tools' own output files: SNP
and SV precision/recall/F1 (hap.py, Truvari refine), switch and Hamming error
rates (yak; the rephased rates are printed alongside), and, for DipGenie, the approximation ratios
rho, rho_cov and rho_ph from `inference/run.tsv`.

```bash
python3 $E/06_evaluation/summarize.py 8x $WORK/runs/VG-HapCUT2-8x \
  $WORK/runs/PanGenie-HapCUT2-8x $WORK/runs/DipGenie/mixed/8x          # SR + Hi-C
```

## 7. SNP error analysis

[snp_errors.py](07_snp_errors/snp_errors.py) counts the SNP TP, FN and FP of evaluated
runs from hap.py's per-record output (`evaluation/happy.vcf.gz`) and splits FN and FP
by the hap.py match type `BK` (`am`: the allele is called with a different genotype;
`lm`: a nearby call; `miss`: no matching call).

The `roofline` (oracle-phase) rows of `phasing/phasing.tsv` (stage 6) complement this:
equal actual and oracle error rates mean that the remaining switch and Hamming error is due to
genotype errors rather than to the orientation of phase blocks.

```bash
bash $E/07_snp_errors/snp_errors.sh mixed    # -> $WORK/snp-errors/mixed.tsv
bash $E/07_snp_errors/snp_errors.sh SR       # -> $WORK/snp-errors/SR.tsv
```

## Tests (no data needed)

```bash
cd data/expts/Chromosome_1
python3 04_dipgenie/test_run.py
python3 06_evaluation/test_evaluate.py
(cd 07_snp_errors && python3 test_snp_errors.py)
(cd 06_evaluation && python3 test_oracle_phase.py)                                # needs pysam
(cd 02_panel && python3 panel_filter_test.py "$(mktemp -u)")                      # needs bcftools
(cd 03_reads && python3 test_paired_region.py && PAIRTOOLS=pairtools \
   python3 -m unittest test_prepare_chr1_sr test_regions test_audit_inputs)        # needs samtools, pairtools
```
