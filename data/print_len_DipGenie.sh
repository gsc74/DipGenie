#!/bin/bash

samples=(
  HG002 HG00438 HG00621 HG00741 HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
  HG01891 HG01928 HG01952 HG01978 HG02080 HG02257 HG02486 HG02559 HG02622 HG02717
  HG02886 HG03540
)

coverages=(2x 4x full)

echo -e "Sample\tCoverage\tAsm_lens_Mb"

for s in "${samples[@]}"; do
  for cov in "${coverages[@]}"; do
    dir="Results/${s}/${s}_${cov}"

    asm_lens="NA"
    if compgen -G "${dir}/full_*.fa" > /dev/null; then
      asm_lens=$(seqkit stats "${dir}"/full_*.fa -T 2>/dev/null \
        | awk 'NR>1 {printf("%.2f|", $5/1e6)}' | sed 's/|$//')
      [[ -z "$asm_lens" ]] && asm_lens="NA"
    fi

    echo -e "${s}\t${cov}\t${asm_lens}"
  done
done
