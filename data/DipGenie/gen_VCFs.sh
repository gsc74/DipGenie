#!/bin/bash
set -euo pipefail

# Generate VCFs from DipGenie haplotype FASTA output via cactus-pangenome

SAMPLES=(HG002 HG00438 HG00621 HG00741 HG03540
         HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
         HG01891 HG01928 HG01952 HG01978
         HG02080 HG02257 HG02486 HG02559 HG02622 HG02717 HG02886)
TAGS=(1x 2x Total)

RESULTS="Results"
HPRC_HAPS="../hprc_haps"
NUMA="numactl --cpunodebind=0 --membind=0"

shopt -s nullglob

for S in "${SAMPLES[@]}"; do
  for T in "${TAGS[@]}"; do
    FULL_FA="${RESULTS}/${S}/${S}_${T}/full.fa"
    [[ -f "$FULL_FA" ]] || { echo "[WARN] Missing: $FULL_FA"; continue; }

    OUT_DIR="${RESULTS}/${S}/${S}_${T}"
    VCF_OUT="${OUT_DIR}/MHC_${S}_${T}.vcf.gz"
    [[ -f "$VCF_OUT" ]] && { echo "[SKIP] $VCF_OUT exists"; continue; }

    echo "[VCF] $S $T"

    # Split full.fa into two haplotype FASTA files
    HAP1="${OUT_DIR}/hap1.fa"
    HAP2="${OUT_DIR}/hap2.fa"
    awk '/^>/{n++} n==1{print > "'"$HAP1"'"} n==2{print > "'"$HAP2"'"}' "$FULL_FA"

    # Rename headers for cactus
    sed -i "s/^>.*/>$S.1/" "$HAP1"
    sed -i "s/^>.*/>$S.2/" "$HAP2"

    # Prepare seqfile for cactus-pangenome
    SEQFILE="${OUT_DIR}/cactus.seqfile"
    cat > "$SEQFILE" <<EOF
CHM13.0 ${HPRC_HAPS}/MHC-CHM13.0.fa
${S}.1 $(realpath "$HAP1")
${S}.2 $(realpath "$HAP2")
EOF

    JOBSTORE="${OUT_DIR}/jobstore"
    rm -rf "$JOBSTORE"

    ${NUMA} cactus-pangenome "$JOBSTORE" "$SEQFILE" \
      --outDir "$OUT_DIR" --outName "MHC_${S}_${T}" \
      --reference CHM13.0 --vcf --gfa \
      --maxCores "$((128))" \
      2> "${OUT_DIR}/cactus.err" || { echo "[ERROR] cactus failed for $S $T"; continue; }

    # Index VCF
    if [[ -f "$VCF_OUT" ]]; then
      bcftools index -t -f "$VCF_OUT"
    fi

    # Cleanup
    rm -rf "$JOBSTORE" "$HAP1" "$HAP2" "$SEQFILE"
  done
done

echo "Done."
