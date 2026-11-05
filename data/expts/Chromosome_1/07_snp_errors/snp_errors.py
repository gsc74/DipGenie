#!/usr/bin/env python3
"""Classify the SNP errors of evaluated runs.

For each run directory (with evaluation/happy.vcf.gz from evaluate.sh) this
prints one TSV row with the hap.py SNP decisions TP, FN and FP, with FN and FP
split by the hap.py match type BK (am: the allele is called but the genotype
differs; lm: a call is nearby; miss: no matching call).

Usage: snp_errors.py RUN [RUN ...]
"""
import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

HAPPY_FORMAT = '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%BD\t%BK\t%BVT]\n'


def alt_indices(gt):
    return {int(a) for a in re.split('[|/]', gt) if a not in ('.', '0')}


def parse_happy(lines):
    """Count SNP decisions and collect FN/FP SNP alleles from hap.py records (TRUTH, QUERY)."""
    counts = dict.fromkeys(['TP', 'FN', 'FP', 'FN_am', 'FN_lm', 'FN_miss', 'FP_am', 'FP_lm', 'FP_miss'], 0)
    fn, fp = set(), set()
    for line in lines:
        f = line.rstrip('\n').split('\t')
        chrom, pos, ref, alt = f[:4]
        (tgt, tbd, tbk, tbvt), (qgt, qbd, qbk, qbvt) = f[4:8], f[8:12]
        alleles = [ref] + alt.split(',')
        truth, query = alt_indices(tgt), alt_indices(qgt)
        if tbvt == 'SNP':
            if tbd == 'TP':
                counts['TP'] += 1
            elif tbd == 'FN':
                counts['FN'] += 1
                counts['FN_' + (tbk if tbk in ('am', 'lm') else 'miss')] += 1
                fn.update((chrom, int(pos), ref, alleles[i]) for i in truth - query
                          if len(ref) == 1 and len(alleles[i]) == 1)
        if qbvt == 'SNP' and qbd == 'FP':
            counts['FP'] += 1
            counts['FP_' + (qbk if qbk in ('am', 'lm') else 'miss')] += 1
            fp.update((chrom, int(pos), ref, alleles[i]) for i in query - truth
                      if len(ref) == 1 and len(alleles[i]) == 1)
    return counts, fn, fp


def classify(args, run):
    happy = os.path.join(run, 'evaluation', 'happy.vcf.gz')
    out = subprocess.run([args.bcftools, 'query', '-f', HAPPY_FORMAT, happy],
                         check=True, capture_output=True, text=True).stdout.splitlines()
    counts, _, _ = parse_happy(out)
    row = {'run': run, **counts}
    return row


COLUMNS = ['run', 'TP', 'FN', 'FP', 'FN_am', 'FN_lm', 'FN_miss', 'FP_am', 'FP_lm', 'FP_miss']


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bcftools', default='bcftools')
    p.add_argument('--jobs', type=int, default=8)
    p.add_argument('runs', nargs='+', help='run directories containing evaluation/happy.vcf.gz')
    args = p.parse_args()
    print('\t'.join(COLUMNS))
    with ProcessPoolExecutor(args.jobs) as pool:
        for row in pool.map(classify, [args] * len(args.runs), args.runs):
            print('\t'.join(str(row[c]) for c in COLUMNS), flush=True)


if __name__ == '__main__':
    sys.exit(main())
