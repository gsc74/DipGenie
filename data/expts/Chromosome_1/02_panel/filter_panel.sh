#!/usr/bin/env bash
# Input has already been region-restricted. Keep sample filtering, trimming
# and count filtering in separate stages: bcftools filters can precede -S.
set -Eeuo pipefail
input=${1:?}; output=${2:?}; keep=${3:?}
# Repair even internally inconsistent source counts before view can consult
# them. Recompute again after subsetting/trimming; GT remains authoritative.
bcftools +fill-tags "$input" --no-version -Ou -- -t AC,AN,AF |
  bcftools view --no-version -S "$keep" -Ou |
  bcftools view --no-version -a -Ou |
  bcftools +fill-tags --no-version -Ou -- -t AC,AN,AF |
  bcftools view --no-version -c 1:nref --threads 8 -Oz -o "$output"
bcftools index --threads 8 "$output"
