#!/usr/bin/env python3
"""Combine an input caller's genotypes with HapCUT2 phase for the shared evaluator.

usage: hapcut2_vcf.py original.vcf hapcut2.phased.VCF out.vcf

Genotypes always come from the ORIGINAL call set (HapCUT2 only phases; it may
write ./. for variants it drops). A heterozygous record takes HapCUT2's phased
GT and PS only when HapCUT2 phased it (PS present, '|' separator) with the same
allele pair. Every other heterozygous record is written UNPHASED (no PS), so
the input caller's own haplotype labels never enter the phase metric.
"""
import json
import sys

original, phased, dst = sys.argv[1:4]


def records(path):
    with open(path) as f:
        for line in f:
            if not line.startswith('#'):
                p = line.rstrip('\n').split('\t')
                yield p


hap = {}
for p in records(phased):
    field = dict(zip(p[8].split(':'), p[9].split(':')))
    hap[(p[0], p[1], p[3], p[4])] = (field.get('GT', '.'), field.get('PS', '.'))

counts = dict(records=0, het=0, het_phased=0, het_unphased=0, missing_from_hapcut2=0, phase_sets=set())
with open(original) as f, open(dst, 'w') as o:
    for line in f:
        if line.startswith('##'):
            if line.startswith('##FORMAT=') and 'ID=GT,' not in line:
                continue
            o.write(line)
            continue
        if line.startswith('#'):
            o.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="HapCUT2 phase set">\n')
            o.write(line)
            continue
        p = line.rstrip('\n').split('\t')
        counts['records'] += 1
        field = dict(zip(p[8].split(':'), p[9].split(':')))
        alleles = field['GT'].replace('|', '/').split('/')
        het = len(alleles) == 2 and alleles[0] != alleles[1] and '.' not in alleles
        if not het:
            p[8], p[9] = 'GT', '/'.join(alleles)
        else:
            counts['het'] += 1
            key = (p[0], p[1], p[3], p[4])
            if key not in hap:
                counts['missing_from_hapcut2'] += 1
            gt, ps = hap.get(key, ('.', '.'))
            ok = ps not in ('.', '') and '|' in gt and sorted(gt.split('|')) == sorted(alleles)
            if ok:
                counts['het_phased'] += 1
                counts['phase_sets'].add(ps)
                p[8], p[9] = 'GT:PS', f'{gt}:{ps}'
            else:
                counts['het_unphased'] += 1
                p[8], p[9] = 'GT', '/'.join(sorted(alleles))
        o.write('\t'.join(p) + '\n')
counts['phase_sets'] = len(counts['phase_sets'])
json.dump(counts, sys.stdout, indent=2)
print()
