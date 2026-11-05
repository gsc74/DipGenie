#!/usr/bin/env python3
"""One depth of the paper's results table as TSV, from the tools' own output files.

SNP P/R/F1: hap.py (evaluation/happy.summary.csv, SNP PASS); SV P/R/F1: Truvari refine
(evaluation/truvari/refine.variant_summary.json); switch and Hamming error rates (paper) and
the oracle-phasing reference (diagnostic only): yak trioeval (phasing/{actual,roofline}.trioeval.txt, W and H lines); approximation
ratios rho, rho_cov, rho_ph for DipGenie runs (inference/run.tsv).

usage: summarize.py <depth> <run dir> [<run dir> ...]   (each run has evaluation/ and phasing/)
"""
import csv
import json
from pathlib import Path
import sys

RATIOS = ('approximation_ratio', 'genotype_approximation_ratio', 'phase_approximation_ratio')


def yak(path, line):
    rows = [l.split('\t') for l in open(path) if l.startswith(line + '\t')]
    return '%.2f' % (100 * float(rows[0][3]))


def row(run):
    snp = [r for r in csv.DictReader(open(run/'evaluation/happy.summary.csv')) if r['Type'] == 'SNP' and r['Filter'] == 'PASS'][0]
    tv = run/'evaluation/truvari'
    sv = json.load(open(tv/'refine.variant_summary.json' if (tv/'refine.variant_summary.json').exists() else tv/'summary.json'))
    out = ['%.4f' % float(snp[k]) for k in ('METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score')]
    out += ['%.4f' % sv[k] for k in ('precision', 'recall', 'f1')]
    for kind in ('actual', 'roofline'):
        out += [yak(run/f'phasing/{kind}.trioeval.txt', 'W'), yak(run/f'phasing/{kind}.trioeval.txt', 'H')]
    tsv = run/'inference/run.tsv'
    kv = dict(l.rstrip('\n').split('\t', 1) for l in open(tsv) if '\t' in l) if tsv.exists() else {}
    return out + [kv.get(k, 'NA') for k in RATIOS]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    print('\t'.join(['depth', 'run', 'SNP_P', 'SNP_R', 'SNP_F1', 'SV_P', 'SV_R', 'SV_F1', 'switch', 'hamming',
                     'switch_roofline', 'hamming_roofline', 'rho', 'rho_cov', 'rho_ph']))
    for r in map(Path, sys.argv[2:]):
        print('\t'.join([sys.argv[1], str(r), *row(r)]))


if __name__ == '__main__':
    main()
