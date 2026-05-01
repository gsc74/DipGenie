# Reproducible Workflow for MHC Haplotype Phasing Benchmarks

This repository contains scripts to reproduce the benchmarks in our paper comparing four haplotype phasing tools on the MHC region.

## Tools Compared

| Folder | Method |
|--------|--------|
| `DipGenie/` | DipGenie (graph-based diplotype assembly) |
| `VG/` | vg haplotypes (graph-based haplotype sampling) |
| `Paragraph + Beagle/` | Paragraph SV genotyping + Beagle phasing |
| `PanGenie + Beagle/` | PanGenie genotyping + Beagle phasing |

## Data

Download `Data.zip` from Zenodo and place it in this directory. Then run:

```bash
bash prepare_data.sh
```

This will unpack:
- `Graphs/` — GFA pangenome graphs (one per leave-one-out sample)
- `hprc_haps/` — HPRC haplotype FASTA files
- `MHC_CHM13.gmap` — Genetic map for the MHC region
- `Reads/` — Simulated reads at 1x, 2x, and Total coverage
- `Truth/` — Ground-truth phased VCFs per sample
- `VCFs/` — Panel VCFs for linear tools (one per leave-one-out sample)

## Samples

22 HPRC samples are used in a leave-one-out evaluation:

HG002, HG00438, HG00621, HG00741, HG03540, HG01106, HG01109, HG01123, HG01258, HG01358, HG01361, HG01891, HG01928, HG01952, HG01978, HG02080, HG02257, HG02486, HG02559, HG02622, HG02717, HG02886

The sample list is defined directly in each `run_tool.sh` script.

## Prerequisites

### Software

- **DipGenie**: Compiled binary (included or built from source)
- **vg** (≥ 1.50): `vg`, `gfa2gbwt`, `kmc`, `seqtk`
- **Paragraph**: Install via conda (`conda activate paragraph`)
- **PanGenie**: Build from source (`pangenie/build/src/PanGenie`)
- **Beagle** (≥ 5.4): `beagle` JAR in PATH
- **Cactus** (≥ 2.6): `cactus-pangenome` for VCF generation
- **bcftools**, **samtools**, **tabix**, **bwa-mem2**
- **whatshap** (for evaluation): `conda activate whatshap`
- **truvari** (for SV evaluation)
- **numactl** (for NUMA-aware execution)
- Python 3 with standard libraries

## Workflow (per tool)

Each tool folder contains independent scripts to be run in order:

### 1. Run the tool

```bash
cd DipGenie/  # or VG/, Paragraph + Beagle/, PanGenie + Beagle/
bash run_tool.sh
```

### 2. Generate VCFs (DipGenie and VG only)

DipGenie and VG output FASTA haplotypes that must be converted to VCF:

```bash
bash gen_VCFs.sh
```

### 3. Evaluate phasing accuracy (whatshap compare)

```bash
bash eval_whatshap.sh
```

Produces `Evaluation/{SAMPLE}/{SAMPLE}_{TAG}/SER.txt` with switch error rate and mismatch rate.

### 4. Evaluate SV accuracy (truvari bench)

```bash
bash eval_truvari.sh
```

Produces `SV_Evaluation/{SAMPLE}/{SAMPLE}_{TAG}/bench/summary.json` with F1, precision, recall.

### 5. Collect results

```bash
bash plot.sh
```

Produces CSV files in `output/` (F1, Precision, Recall, SER, Mismatch).

## Coverage Levels

| Label | Description |
|-------|-------------|
| 1x | ~1x coverage of MHC region |
| 2x | ~2x coverage of MHC region |
| Total | Full coverage (max: ~12x) |

## Notes

- Linear tools (Paragraph+Beagle, PanGenie+Beagle) use paired-end reads (R1/R2).
- Graph tools (DipGenie, VG) use interleaved single-file reads.
- All scripts assume they are run from within their respective tool folder.
- The `prepare_data.sh` script must be run first to unpack data into the shared directories.

