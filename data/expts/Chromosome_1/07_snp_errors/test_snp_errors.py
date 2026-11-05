#!/usr/bin/env python3
"""Unit test of the hap.py record classification (no data or bcftools needed)."""
import snp_errors


def record(pos, ref, alt, truth, query):
    return '\t'.join(['chr1', str(pos), ref, alt, *truth, *query])


lines = [
    record(10, 'A', 'C', ('0|1', 'TP', 'am', 'SNP'), ('1|0', 'TP', 'am', 'SNP')),
    record(20, 'A', 'G', ('1|1', 'FN', 'am', 'SNP'), ('0|1', 'FP', 'am', 'SNP')),   # genotype error
    record(30, 'C', 'T', ('0|1', 'FN', '.', 'SNP'), ('.', 'N', '.', 'NOCALL')),       # missed allele
    record(40, 'G', 'A,T', ('0|1', 'FN', 'lm', 'SNP'), ('0|2', 'FP', 'lm', 'SNP')),   # wrong allele
    record(50, 'T', 'TA', ('0|1', 'FN', '.', 'INDEL'), ('.', 'N', '.', 'NOCALL')),    # not a SNP
]
counts, fn, fp = snp_errors.parse_happy(lines)
assert counts == {'TP': 1, 'FN': 3, 'FP': 2, 'FN_am': 1, 'FN_lm': 1, 'FN_miss': 1,
                  'FP_am': 1, 'FP_lm': 1, 'FP_miss': 0}, counts
assert fn == {('chr1', 30, 'C', 'T'), ('chr1', 40, 'G', 'A')}, fn
assert fp == {('chr1', 40, 'G', 'T')}, fp
print('test_snp_errors: OK')
