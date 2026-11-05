# Usage

Build instructions and a quick start are in the [README](../../README.md).

## Input

- A sequence-labelled GFA (v1.1) with complete panel haplotypes as `P` (path) or `W` (walk, PanSN `sample#hap#contig`) lines and zero-overlap links; acyclic and cyclic graphs are supported. Use `--REF` to identify the reference path or sample prefix, resolving to one path per component.
- A phased, non-overlapping, sequence-resolved panel VCF and its reference FASTA. Symbolic alleles and overlapping records must be resolved during panel preparation.
- A saved graph index. The complete reference is stored in the index; no separate reference argument is needed for inference or dipcall calling.

We recommend building an index before inference:

```bash
./DipGenie --gfa panel.gfa --REF CHM13#0 --write-index panel.gindex
./DipGenie --vcf panel.vcf.gz --reference reference.fa --write-index panel.gindex
```

These are alternative commands; choose a new index directory for each panel. The index stores graph shards and reference sequences.

## Inference

Short reads:

```bash
./DipGenie --gindex panel.gindex --sr -r1 SR.R1.fastq.gz -r2 SR.R2.fastq.gz -t8 -R18 --out run
```

Hi-C:

```bash
./DipGenie --gindex panel.gindex --hic -h1 HiC.R1.fastq.gz -h2 HiC.R2.fastq.gz -t8 -R18 --out run
```

Both libraries:

```bash
./DipGenie --gindex panel.gindex --sr -r1 SR.R1.fastq.gz -r2 SR.R2.fastq.gz --hic -h1 HiC.R1.fastq.gz -h2 HiC.R2.fastq.gz -t8 -R18 --out run
```

Use separate `-t 8 -R 18` arguments if preferred. `-R` is a **shared recombination budget for the two haplotypes of each component**, not a per-haplotype budget; the emitted haplotypes can contain more switches than `-R`, because exchanging the two haplotypes is not counted. The default is 18; chromosome-scale runs should set a budget for their chromosome. Long-read modes are not supported.

DipGenie places paired reads on the graph itself; no read mapper, SAM/BAM or
contact file is needed.

Both mates must be in matching order and have matching molecule identifiers. FASTQ sequences and qualities are preserved. Two read inputs are supported:

- **One component** (e.g. a single-chromosome panel): supply reads already restricted to that chromosome; they are used as given, without binning.
- **Several components** (e.g. a whole-genome panel): supply all reads. The index stores the complete reference sequence of every component, and both mates are localized against one combined index of these references. Within-component pairs enter that component's bin; cross-component and unassigned pairs are retained separately and are not used as within-component phase links. Components with no assigned pairs are reported as `NO_READ_PAIRS`, not fabricated haplotypes.

No external reference or read mapper is needed in either case.

## FASTA and VCF output

Add `--output-vcf calls.vcf.gz --sample HG002` to the inference command to generate FASTA and phased variant calls in one invocation:

```bash
./DipGenie --gindex panel.gindex --sr -r1 SR.R1.fastq.gz -r2 SR.R2.fastq.gz -t8 -R18 --out run --sample HG002 --output-vcf calls.vcf.gz
```

The dipcall 0.3 assembly stages are linked into the binary: minimap2 2.31
assembly alignment, dipcall alignment filtering, HTSbox pileup, diploid GT
conversion, and callable BED generation. `make` builds these pinned sources.
There is no runtime dependency on run-dipcall, Perl, k8, samtools, bedtk or an
external minimap2/HTSbox executable. Normalization, PASS/called-allele filtering,
multi-component joining and CSI indexing also run in process, without bcftools.
The normalizer handles the internal caller's literal DNA alleles and GT/AD
fields; it is not an arbitrary third-party VCF transformation interface.
Haplotype order and sample identity are preserved. Raw calls,
callable BEDs, backend versions and logs remain in a private sibling
`.dipgenie-dipcall-*` directory. Calling never uses a truth BED. This is the
diploid all-contigs policy; no sex-specific PAR policy is inferred silently.

The primary outputs are `all-components.fa`, `run.tsv` (per component for multi-component inputs), and the optional VCF. Output directories must be new. Logs and intermediate inputs are retained on failure.

## Performance

Work within each component is parallel, and the result is the same for any
number of threads. Components are processed one after another. The search is
heuristic; `approximation_ratio` in `run.tsv` bounds how far the returned pair
can be from the optimum of the objective (`NULL` if no finite ratio is
available). It is not a guarantee about F1 or phasing accuracy.

## Tests and experiments

`make check` runs exact-score and bound tests, serial/parallel parity, raw SR/Hi-C/mixed DAG and cyclic toys, paired component routing, VCF/GFA/index round trips, embedded assembly-caller fixtures and invalid-input guards. It also runs a complete raw-input-to-VCF test with an empty `PATH`:

```bash
python3 test/joint/raw_toy_test.py --executable ./DipGenie --libraries mixed --recombinant --dipcall
```

Optional caller checks compare against the dipcall kit and bcftools, used only
as test oracles:

```bash
python3 test/joint/assembly_parity.py build/test_assembly_calling --kit /path/to/dipcall.kit
python3 test/joint/variant_output_parity.py build/test_assembly_calling --bcftools /path/to/bcftools
```

Experiment scripts are under [data/expts](../expts). [Chromosome_1](../expts/Chromosome_1) runs the HG002 chromosome-1 benchmark (HG002 chr1, HPRC v2 panel, SR + Hi-C, DipGenie v2 versus VG/PanGenie + HapCUT2), from download to evaluation.
