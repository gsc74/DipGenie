# Experiments

Build and test DipGenie using the [repository instructions](../../README.md).

- [Chromosome_1](Chromosome_1/README.md): HG002 chromosome 1 benchmark.
  It covers the download of the HPRC v2 panel, reads, references and GIAB truth;
  leave-one-out panel filtering; nested SR and Hi-C subsampling (2x to Total);
  DipGenie v2; VG + HapCUT2 and PanGenie + HapCUT2; evaluation; and the
  SNP error analysis.

Full runs need the external data listed in
[sources.tsv](Chromosome_1/01_download/sources.tsv), the comparison tools in
[config.example.sh](Chromosome_1/config.example.sh), and a Slurm allocation.
Large inputs and outputs are not committed.
