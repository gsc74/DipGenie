# DipGenie v2 performance

## Setup

We evaluated DipGenie v2 on chromosome 1 of HG002 using a leave-one-out design: HG002 was removed from the [HPRC release 2](https://humanpangenome.org) panel VCF, leaving 230 samples (461 haplotypes with the CHM13 reference).

- **Short reads:** Illumina NovaSeq 6000, 2×250 bp, SRA [SRR11321732](https://www.ncbi.nlm.nih.gov/sra/SRR11321732).
- **Hi-C reads:** 2×251 bp, HPRC HG002 [R1](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/hic/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R1_001.fastq.gz) and [R2](https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/hic/HG002.HiC_2_NovaSeq_rep1_run2_S1_L001_R2_001.fastq.gz).
- Reads on chromosome 1 were subsampled to nested sets of 2×, 4×, 8×, 16× and 32× per library; Total is all reads (58.7× short reads, 56.2× Hi-C).
- DipGenie v2 is compared with [VG](https://github.com/vgteam/vg) and [PanGenie](https://github.com/eblerjana/pangenie) genotypes phased by [HapCUT2](https://github.com/vibansal/HapCUT2) from the same reads.

Metrics:
- SNP F1: [`hap.py`](https://github.com/Illumina/hap.py) with vcfeval against GIAB v5.0q.
- SV F1: [`truvari`](https://github.com/ACEnglish/truvari) bench and refine against GIAB v5.0q.
- Switch and Hamming error rates (%) of the two haplotypes: [`yak`](https://github.com/lh3/yak) `trioeval` against the [HG002 Q100](https://github.com/marbl/hg002) assembly. The switch error rate measures local phase and the Hamming error rate the orientation of phase blocks along the chromosome (about 50% for random orientation); lower is better.

Best value per depth in bold.

## Accuracy, short reads + Hi-C

DipGenie v2 has the lowest switch and Hamming error rates and the highest SV F1 at every depth, and the highest SNP F1 up to 8×.

| Depth | Method | SNP F1 | SV F1 | Switch error (%) | Hamming error (%) |
|:--|:--|--:|--:|--:|--:|
| 2× | VG+HapCUT2 | 0.3702 | 0.5959 | 16.69 | 49.93 |
| 2× | PanGenie+HapCUT2 | 0.7690 | 0.5583 | 24.68 | 49.64 |
| 2× | DipGenie v2 | **0.9165** | **0.7886** | **6.49** | **48.47** |
| 4× | VG+HapCUT2 | 0.7603 | 0.7025 | 24.53 | 49.59 |
| 4× | PanGenie+HapCUT2 | 0.8606 | 0.5867 | 25.77 | 49.87 |
| 4× | DipGenie v2 | **0.9311** | **0.8042** | **6.29** | **42.76** |
| 8× | VG+HapCUT2 | 0.9040 | 0.7528 | 24.69 | 49.60 |
| 8× | PanGenie+HapCUT2 | 0.9275 | 0.6131 | 24.59 | 36.96 |
| 8× | DipGenie v2 | **0.9372** | **0.7960** | **6.37** | **11.21** |
| 16× | VG+HapCUT2 | 0.9423 | 0.7750 | 21.45 | 22.77 |
| 16× | PanGenie+HapCUT2 | **0.9528** | 0.6290 | 20.79 | 22.32 |
| 16× | DipGenie v2 | 0.9407 | **0.7976** | **6.25** | **7.49** |
| 32× | VG+HapCUT2 | 0.9530 | 0.7808 | 18.15 | 14.45 |
| 32× | PanGenie+HapCUT2 | **0.9565** | 0.6361 | 17.56 | 15.04 |
| 32× | DipGenie v2 | 0.9420 | **0.7893** | **6.05** | **5.70** |
| Total | VG+HapCUT2 | 0.9571 | 0.7846 | 16.56 | 12.34 |
| Total | PanGenie+HapCUT2 | **0.9595** | 0.6390 | 16.31 | 12.55 |
| Total | DipGenie v2 | 0.9435 | **0.7970** | **5.63** | **4.73** |

## Accuracy, short reads only

DipGenie v2 has the lowest switch error rate and the highest SV F1 at every depth, and the highest SNP F1 up to 8×. Without Hi-C, no method orients the phase blocks along the chromosome (Hamming error rates of 47–50%).

| Depth | Method | SNP F1 | SV F1 | Switch error (%) | Hamming error (%) |
|:--|:--|--:|--:|--:|--:|
| 2× | VG+HapCUT2 | 0.3702 | 0.5959 | 16.69 | 49.93 |
| 2× | PanGenie+HapCUT2 | 0.7690 | 0.5601 | 23.93 | 49.79 |
| 2× | DipGenie v2 | **0.9175** | **0.7919** | **6.14** | **48.18** |
| 4× | VG+HapCUT2 | 0.7603 | 0.7029 | 23.93 | 49.79 |
| 4× | PanGenie+HapCUT2 | 0.8606 | 0.5886 | 24.85 | 49.88 |
| 4× | DipGenie v2 | **0.9310** | **0.8032** | **5.56** | **49.60** |
| 8× | VG+HapCUT2 | 0.9040 | 0.7511 | 23.89 | 49.60 |
| 8× | PanGenie+HapCUT2 | 0.9275 | 0.6113 | 23.66 | 49.83 |
| 8× | DipGenie v2 | **0.9364** | **0.7983** | **5.51** | **48.22** |
| 16× | VG+HapCUT2 | 0.9423 | 0.7754 | 22.62 | 49.79 |
| 16× | PanGenie+HapCUT2 | **0.9528** | 0.6290 | 22.44 | **49.29** |
| 16× | DipGenie v2 | 0.9406 | **0.8012** | **5.39** | 49.58 |
| 32× | VG+HapCUT2 | 0.9530 | 0.7801 | 22.06 | 49.83 |
| 32× | PanGenie+HapCUT2 | **0.9565** | 0.6345 | 21.40 | 49.62 |
| 32× | DipGenie v2 | 0.9420 | **0.7918** | **5.39** | **47.20** |
| Total | VG+HapCUT2 | 0.9571 | 0.7850 | 21.38 | 49.33 |
| Total | PanGenie+HapCUT2 | **0.9595** | 0.6354 | 21.14 | 49.30 |
| Total | DipGenie v2 | 0.9435 | **0.8001** | **5.36** | **48.35** |

## Runtime and memory, short reads + Hi-C

| Depth | Method | Wall-clock time (h:mm:ss) | Peak memory (GiB) |
|:--|:--|--:|--:|
| 2× | VG+HapCUT2 | 0:13:47 | 39.6 |
| 2× | PanGenie+HapCUT2 | 0:09:45 | 39.5 |
| 2× | DipGenie v2 | 1:44:47 | 103.1 |
| 4× | VG+HapCUT2 | 0:14:22 | 48.8 |
| 4× | PanGenie+HapCUT2 | 0:11:31 | 48.7 |
| 4× | DipGenie v2 | 1:49:51 | 105.7 |
| 8× | VG+HapCUT2 | 0:16:45 | 61.7 |
| 8× | PanGenie+HapCUT2 | 0:12:21 | 61.7 |
| 8× | DipGenie v2 | 2:04:07 | 108.7 |
| 16× | VG+HapCUT2 | 0:20:22 | 85.8 |
| 16× | PanGenie+HapCUT2 | 0:16:56 | 86.1 |
| 16× | DipGenie v2 | 2:20:37 | 113.9 |
| 32× | VG+HapCUT2 | 0:30:18 | 131.5 |
| 32× | PanGenie+HapCUT2 | 0:28:04 | 131.5 |
| 32× | DipGenie v2 | 2:49:58 | 122.2 |
| Total | VG+HapCUT2 | 0:42:27 | 136.4 |
| Total | PanGenie+HapCUT2 | 0:37:17 | 136.0 |
| Total | DipGenie v2 | 3:33:39 | 130.8 |

## Runtime and memory, short reads only

| Depth | Method | Wall-clock time (h:mm:ss) | Peak memory (GiB) |
|:--|:--|--:|--:|
| 2× | VG+HapCUT2 | 0:12:38 | 39.6 |
| 2× | PanGenie+HapCUT2 | 0:08:57 | 39.5 |
| 2× | DipGenie v2 | 1:52:15 | 103.9 |
| 4× | VG+HapCUT2 | 0:13:06 | 48.8 |
| 4× | PanGenie+HapCUT2 | 0:10:03 | 48.7 |
| 4× | DipGenie v2 | 1:58:02 | 105.1 |
| 8× | VG+HapCUT2 | 0:13:38 | 61.7 |
| 8× | PanGenie+HapCUT2 | 0:09:34 | 61.7 |
| 8× | DipGenie v2 | 2:02:08 | 109.0 |
| 16× | VG+HapCUT2 | 0:15:04 | 85.8 |
| 16× | PanGenie+HapCUT2 | 0:11:19 | 86.1 |
| 16× | DipGenie v2 | 2:08:00 | 110.6 |
| 32× | VG+HapCUT2 | 0:18:30 | 131.5 |
| 32× | PanGenie+HapCUT2 | 0:14:59 | 131.5 |
| 32× | DipGenie v2 | 2:22:02 | 115.2 |
| Total | VG+HapCUT2 | 0:23:52 | 136.4 |
| Total | PanGenie+HapCUT2 | 0:20:19 | 136.0 |
| Total | DipGenie v2 | 2:47:52 | 122.0 |

All runs used 64 physical cores of an Intel Xeon Platinum 8592+ processor with 512 GB of memory. Times are wall-clock times per run, from the reads to the phased call set; index construction and evaluation are excluded.

The scripts to reproduce these results, from data download to evaluation, are in [data/expts/Chromosome_1](data/expts/Chromosome_1).
