#!/bin/bash
set -euo pipefail

# Download and unpack all datasets from Zenodo

ZENODO_URL="https://zenodo.org/records/19939755/files/Data.zip?download=1"

if [[ ! -f Data.zip ]]; then
  echo "[0/7] Downloading Data.zip from Zenodo..."
  wget -q --show-progress -O "Data.zip" "$ZENODO_URL"
fi

echo "[1/7] Unpacking Data.zip..."
unzip -o Data.zip

echo "[2/7] Unpacking Graphs..."
unzip -o Graphs.zip -d Graphs

echo "[3/7] Unpacking Truth VCFs..."
unzip -o Truth.zip

echo "[4/7] Unpacking haplotype assemblies..."
unzip -o hprc_haps.zip

echo "[5/7] Unpacking Reads..."
unzip -o Reads.zip

echo "[6/7] Unpacking VCFs (panel VCFs)..."
unzip -o VCFs.zip -d VCFs

echo "[7/7] Unpacking genetic map..."
unzip -o MHC_CHM13.gmap.zip

echo ""
echo "=== Data layout ==="
echo "Graphs/MHC_left_{SAMPLE}.gfa   — Leave-one-out pangenome graphs (22 samples)"
echo "Truth/{SAMPLE}/MHC_{SAMPLE}.vcf.gz — Ground truth VCFs"
echo "hprc_haps/MHC-{SAMPLE}.{1,2}.fa   — Haplotype assemblies"
echo "hprc_haps/MHC-CHM13.0.fa          — CHM13 reference"
echo "Reads/{1x,2x,Total}/              — Illumina paired-end reads"
echo "VCFs/{SAMPLE}.vcf                  — Panel VCFs (vcfbub-decomposed)"
echo "MHC_CHM13.gmap                     — Genetic map for MHC region"
echo ""
echo "Done."
