# Copy to config.sh, edit, then `export CONFIG=/abs/path/config.sh` before running any stage.
REPO=/path/to/DipGenie                      # this repository (built with make)
WORK=/path/to/chr1-study                    # downloads, intermediate files and results
THREADS=64
PYTHON=python3                              # Python 3 with pysam 0.24.1
SAMPLE=HG002
DIPGENIE=${DIPGENIE:-$REPO/DipGenie}
LABEL=${LABEL:-DipGenie}                    # run directory under $WORK/runs

# Tools used by the benchmark (versions of the evaluated runs). DipGenie itself needs none of them;
# it calls variants in process (dipcall v0.3 stages, minimap2 v2.31, HTSbox r352).
BCFTOOLS=bcftools                           # 1.25
SAMTOOLS=samtools                           # 1.25
SAMTOOLS_HIC=$SAMTOOLS
BWA_MEM2=bwa-mem2                           # 2.3 (baseline alignments)
PAIRTOOLS=pairtools                         # 1.1.2 (Hi-C duplicate marking in read preparation)
VG=vg                                       # 1.78.0
KMC=kmc                                     # 3.2.4
PANGENIE_DIR=/path/to/pangenie             # PanGenie v4.2.2 release binaries (PanGenie, PanGenie-index)
HAPCUT2_DIR=/path/to/HapCUT2/build          # HapCUT2 v1.3.4 (extractHAIRS, HAPCUT2)
RTG=rtg                                     # RTG Tools 3.13
TRUVARI=truvari                             # 5.5.0
YAK=yak                                     # yak r93 (https://github.com/lh3/yak)
PREFETCH=prefetch                           # SRA Toolkit (download only)
FASTERQ_DUMP=fasterq-dump
PIGZ=pigz                                   # 2.8
# hap.py 0.3.15 and MAFFT 7.526 must be on PATH for evaluation.
