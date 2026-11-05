#!/usr/bin/env bash
# VG and PanGenie indexes of the same HG002-excluded chr1 panel (built once).
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/05_baselines
# Build the same embedded caller used by DipGenie, including its pinned static
# libraries. No standalone run-dipcall/bcftools dependency for this helper.
make -C "$REPO" -j"$THREADS" build/call_fasta
panel=$WORK/panel/chr1.noHG002.vcf.gz
ref=$WORK/panel/chr1.fa
vgdir=$WORK/vg-index; mkdir "$vgdir"
"$VG" construct -v "$panel" -r "$ref" -aS -f -i -p > "$vgdir/panel.vg"
"$VG" convert -x "$vgdir/panel.vg" > "$vgdir/panel.xg"
"$VG" gbwt -x "$vgdir/panel.xg" -v "$panel" --num-jobs "$THREADS" -o "$vgdir/haplotypes.gbwt"
"$VG" gbwt -x "$vgdir/panel.xg" -E -o "$vgdir/reference.gbwt"
"$VG" gbwt -m "$vgdir/haplotypes.gbwt" "$vgdir/reference.gbwt" -o "$vgdir/panel.gbwt"
"$VG" gbwt -x "$vgdir/panel.xg" "$vgdir/panel.gbwt" -g "$vgdir/panel.gbz"
"$VG" index -t "$THREADS" --no-nested-distance -j "$vgdir/panel.dist" "$vgdir/panel.gbz"
"$VG" gbwt --num-threads "$THREADS" -r "$vgdir/panel.ri" -Z "$vgdir/panel.gbz"
"$VG" haplotypes -v2 -t "$THREADS" -d "$vgdir/panel.dist" -H "$vgdir/panel.hapl" "$vgdir/panel.gbz"
# Freeze the compiled helper with the indexes; record its exact producing code.
install -m755 "$REPO/build/call_fasta" "$vgdir/call_fasta"
sha256sum "$vgdir/call_fasta" > "$vgdir/caller.sha256"
git -C "$REPO" rev-parse HEAD > "$vgdir/caller.revision"
touch "$vgdir/COMPLETE"

pgdir=$WORK/pangenie-index; mkdir "$pgdir"
"$BCFTOOLS" view -Ov -o "$pgdir/panel.vcf" "$panel"
"$PANGENIE_DIR/PanGenie-index" -r "$ref" -v "$pgdir/panel.vcf" -t "$THREADS" -e 3000000000 -o "$pgdir/panel" \
  > "$pgdir/index.stdout" 2> "$pgdir/index.log"
rm -f "$pgdir/panel.vcf"
touch "$pgdir/COMPLETE"
