#!/usr/bin/env bash
# DipGenie v2 at one depth (manifest row), then evaluation.
# MODE=mixed (SR + Hi-C; default), MODE=SR or MODE=HiC.
# Paper settings: 64 threads on one NUMA node, R = 2140, DipGenie defaults.
# Usage: [MODE=SR] bash run.sh <index>   (index 0..5 = 2x 4x 8x 16x 32x Total)
# Output: $WORK/runs/<LABEL>/<mode>/<depth> (default LABEL=DipGenie).
set -Eeuo pipefail
source "${CONFIG:?export CONFIG=/path/to/config.sh (copy of config.example.sh)}"
here=$REPO/data/expts/Chromosome_1/04_dipgenie
command -v numactl >/dev/null
task=${1:?usage: run.sh <index 0-5>}
mode=${MODE:-mixed}; [[ $mode == mixed || $mode == SR || $mode == HiC ]]
cpus=$("$PYTHON" "$here/run.py" --cpu-list --threads 64)
node=$(lscpu -p=CPU,NODE | awk -F, -v cpu="${cpus%%,*}" '$1==cpu {print $2}')
export OMP_NUM_THREADS=64 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
numactl --physcpubind="$cpus" --membind="$node" "$PYTHON" "$here/run.py" \
  --binary "$DIPGENIE" --gindex "$WORK/panel/chr1.gindex" --inputs "$WORK/reads/manifest.tsv" \
  --output "$WORK/runs/${LABEL:-DipGenie}" --task "$task" --mode "$mode" --threads 64 --budget 2140 \
  --sample "$SAMPLE"
depth=$(awk -F'\t' -v t="$task" 'NR==t+2 {print $1}' "$WORK/reads/manifest.tsv")
run=$WORK/runs/${LABEL:-DipGenie}/$mode/$depth
bind=(numactl --physcpubind="$cpus" --membind="$node")
"${bind[@]}" --show > "$run/placement.txt"
# Independently replay the final dosage, phase factors and objective certificate.
# Validation/evaluation time is separate from inference + calling resources.txt.
/usr/bin/time -v -o "$run/replay.resources.txt" "${bind[@]}" "$PYTHON" \
  "$REPO/test/joint/audit_model.py" "$run/inference" --output "$run/score-replay.json" \
  > "$run/score-replay.log" 2>&1
/usr/bin/time -v -o "$run/evaluation.resources.txt" "${bind[@]}" bash \
  "$here/../06_evaluation/evaluate.sh" "$run/calls.vcf.gz" "$run/evaluation"
/usr/bin/time -v -o "$run/phasing.resources.txt" "${bind[@]}" bash \
  "$here/../06_evaluation/phasing.sh" "$run/calls.vcf.gz" "$run/phasing"
touch "$run/EVALUATED"
