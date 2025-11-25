#!/bin/bash

samples=(
HG002 HG00438 HG00621 HG00741 HG01106 HG01109 HG01123 HG01258 HG01358 HG01361
HG01891 HG01928 HG01952 HG01978 HG02080 HG02257 HG02486 HG02559 HG02622 HG02717
HG02886 HG03540
)

coverages=(2x 4x full)

echo -e "Sample\tCoverage\tSV_count"

for s in "${samples[@]}"; do
  for cov in "${coverages[@]}"; do
    vcf="Results_VG/${s}/${s}_${cov}/MHC_${s}_${cov}.vcf.gz"
    if [[ -f "$vcf" ]]; then
      count=$(bcftools query -f '%REF\t%ALT\n' "$vcf" 2>/dev/null \
        | awk '{split($2,a,","); for(i in a) if (length(a[i]) - length($1) >= 50 || length($1) - length(a[i]) >= 50) c++} END{print c+0}')
      echo -e "${s}\t${cov}\t${count}"
    else
      echo -e "${s}\t${cov}\tNA"
    fi
  done
done

