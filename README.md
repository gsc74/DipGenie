[![BioConda Install](https://img.shields.io/conda/dn/bioconda/dipgenie?label=BioConda%20Installs)](https://anaconda.org/bioconda/dipgenie)


## <div align="center"><span style="color:red;"><b>DipGenie</b></span> (<span style="color:red;"><b>D</b></span>iploid <span style="color:red;"><b>G</b></span>enome <span style="color:red;"><b>i</b></span>nference)</div>

## <a name="started"></a>Getting Started

### Prerequisites
- **Zlib** - [zlib](https://zlib.net/)
- A C++17 compiler with OpenMP, Make and Boost headers; Python 3 for the tests

Or install them with `conda env create -f environment.yml`.


## Table of Contents

- [Get DipGenie](#get_dipgenie)
- [Usage](#usage)
- [Introduction](#intro)
- [Results](#results)
- [Future work](#future)
- [Publications](#pub)


## <a name="get_dipgenie"></a>Get DipGenie

```bash
git clone --recursive https://github.com/gsc74/DipGenie
cd DipGenie
make -j8
make check

# Test run (MHC, HG002 paired short reads at 2x)
# 1. Build a graph index from a phased VCF panel (once per panel)
./DipGenie --vcf test/MHC_4.vcf.gz --reference test/MHC-CHM13.0.fa.gz --write-index test/MHC_4.gindex -t32
# 2. Infer the two haplotypes (FASTA) and phased variants (VCF)
./DipGenie --gindex test/MHC_4.gindex --sr -r1 test/HG002.mhc.R1.fq.gz -r2 test/HG002.mhc.R2.fq.gz \
  -t32 -R18 --out HG002 --sample HG002 --output-vcf HG002.vcf.gz
```


## <a name="usage"></a>Usage

```bash
# Short reads only
./DipGenie --gindex panel.gindex --sr -r1 SR.R1.fq.gz -r2 SR.R2.fq.gz -t64 -R2140 --out run
# Hi-C only
./DipGenie --gindex panel.gindex --hic -h1 HiC.R1.fq.gz -h2 HiC.R2.fq.gz -t64 -R2140 --out run
# Short reads + Hi-C, with phased VCF output
./DipGenie --gindex panel.gindex --sr -r1 SR.R1.fq.gz -r2 SR.R2.fq.gz \
  --hic -h1 HiC.R1.fq.gz -h2 HiC.R2.fq.gz -t64 -R2140 --out run --sample HG002 --output-vcf calls.vcf.gz
```

**Pangenome input** (indexed once with `--write-index`, then reused with `--gindex`):
- `--vcf panel.vcf.gz --reference ref.fa`: a phased, multi-allelic VCF panel and its reference.
- `--gfa panel.gfa --REF CHM13#0`: a GFA v1.1 graph with haplotypes as `P` or `W` lines (acyclic or cyclic); `--REF` names the reference path.

**Reads** are paired four-line FASTQ files (plain or gzipped), one file per mate:
- `-r1`/`-r2` (short reads) and `-h1`/`-h2` (Hi-C) list the same pairs in the same order.
- Mates have the same name, optionally ending in `/1` and `/2` (or with Illumina comments `1:N:...`, `2:N:...`); names are unique within a file.
- Reads are at most 1,024 bp. Interleaved or single-end files are not accepted.
- For a one-chromosome panel, give reads of that chromosome; for a multi-chromosome panel, give all reads (DipGenie assigns them to chromosomes).

**Options**

| Option | Meaning |
|---|---|
| `--out DIR` | New output directory |
| `-t N` | Threads (default 1) |
| `-R N` | Recombination budget per chromosome (default 18, see below) |
| `--sample NAME` | Sample name in the FASTA and VCF (default `SAMPLE`); the reference is not a sample |
| `--output-vcf FILE` | Also write phased variant calls (`.vcf` or `.vcf.gz`; gzipped output is indexed) |

**Choosing `-R`.** The default suits short regions such as the MHC test; for a whole chromosome set `-R` yourself. We use the expected number of switches between panel haplotypes under the Li and Stephens model, doubled:
E[R] = (L − 1)(H − 1)/H · (1 − exp(−4·Ne·r/H)), R = 2·round(E[R]), with chromosome length L, H = 51 panel haplotypes, Ne = 5,000 and recombination rate r per bp.
[scripts/compute_R.py](scripts/compute_R.py) computes it; for CHM13 chromosome 1 it gives 2140:

```bash
python3 scripts/compute_R.py --length 248387328 --rate 1.12e-8    # 2140
```

Here r = 1.12 cM/Mb is the human chromosome 1 rate from Table 1 of [Jensen-Seaman et al. 2004](https://doi.org/10.1101/gr.1970304), which lists the rate of every human chromosome (divide cM/Mb by 10⁸ for r per bp). For other chromosomes or species, use their length and recombination rate. A multi-chromosome index uses the same `-R` for every chromosome.

**Outputs**: `run/all-components.fa` (two haplotypes per chromosome), `run/run.tsv` (scores and approximation ratio), and with `--output-vcf` the phased VCF. See [docs/usage.md](data/docs/usage.md) for details.


## <a name="intro"></a>Introduction
DipGenie infers the two haplotypes of a diploid genome from a pangenome graph and paired-end short reads, Hi-C reads, or both. Read coverage determines which alleles are present and in what dosage, while read pairs and Hi-C contacts link alleles and so determine their phase. Both are scored together while the two haplotypes are chosen as mosaics of the panel haplotypes.

DipGenie implements Frequency-aware Diploid Path Inference (FDPI). A joint dynamic program over segments shared by all panel haplotypes extends both haplotypes together, supports acyclic and cyclic graphs, and reports an instance-specific approximation ratio from an upper bound on the optimum. It outputs both haplotypes as FASTA and, optionally, phased variant calls; variant calling is built in, so no other programs are needed at run time.


## <a name="results"></a>Results
We evaluated DipGenie v2 on chromosome 1 of HG002 with the HPRC release 2 panel (HG002 left out), using Illumina short reads and Hi-C reads at six depths per library (2× to 58.7× short reads and to 56.2× Hi-C), against VG and PanGenie genotypes phased by HapCUT2.

**Accuracy.** With short reads and Hi-C, DipGenie v2 has the lowest switch and Hamming error rates and the highest SV F1 at every depth, and the highest SNP F1 up to 8×. With all reads, its switch and Hamming error rates are 5.6% and 4.7%, against at best 16.3% and 12.3% for the pipelines. With short reads only, it keeps the lowest switch error rate and the highest SV F1 at every depth, but no method orients the phase blocks along the chromosome.

**Runtime.** DipGenie v2 infers both haplotypes of chromosome 1 in 1.7–3.6 h on 64 cores with at most 131 GiB of memory, against 9–42 min for the pipelines.

See the full performance data, with accuracy and runtime tables for both read combinations, in [perf.md](perf.md).


## <a name="future"></a>Future Work
- Improve the accuracy of genome inference.
- Design faster algorithms for genome inference.


## <a name="pub"></a>Publications
- **Ghanshyam Chandra, William T. Doan and Daniel Gibney**. Haplotype-resolved diploid genome inference on pangenome graphs. *bioRxiv* (doi: [10.1101/2025.11.26.690754](https://doi.org/10.1101/2025.11.26.690754))
- **William T. Doan and Ghanshyam Chandra**. Linkage-aware diploid genome inference using pangenome graphs. *Manuscript in preparation.*


## License

See [LICENSE](LICENSE) for the DipGenie notice and third-party licenses.
