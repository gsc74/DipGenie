#!/usr/bin/env python3
"""Oracle-phase input: a call set with the truth phase on its correctly genotyped heterozygous calls.

Both VCFs must be split to one ALT per record (bcftools norm -m -any). A query
heterozygote whose (CHROM, POS, REF, ALT) is a truth heterozygote takes the
truth's orientation. Every other heterozygous call is flipped exactly when the
preceding matched call was flipped, so its phase relative to that call is the
method's own and the oracle adds no switch the method did not make. Genotypes
are never changed (false calls, wrong genotypes and homozygous calls keep
their alleles); only phase differs from the input. Missing alleles in the truth
GT (e.g. 0|.) are read as reference.

usage: oracle_phase.py truth.norm.vcf.gz query.norm.vcf.gz out.vcf
"""
import sys

import pysam


def gt_of(record):
    return tuple(0 if a is None else a for a in record.samples[0]['GT'])


def het(gt):
    return len(gt) == 2 and sorted(gt) == [0, 1]


def truth_orientation(path):
    with pysam.VariantFile(path) as truth:
        return {(r.chrom, r.pos, r.ref, r.alts[0]): gt_of(r) for r in truth if r.alts and het(gt_of(r))}


def transfer(truth_path, query_path, out_path):
    truth = truth_orientation(truth_path)
    calls = matched = carried = 0
    flipped = False
    with pysam.VariantFile(query_path) as query, pysam.VariantFile(out_path, 'w', header=query.header) as out:
        for r in query:
            sample = r.samples[0]
            gt = gt_of(r)
            if het(gt):
                calls += 1
                key = (r.chrom, r.pos, r.ref, r.alts[0])
                if key in truth:
                    flipped = truth[key] != gt
                    gt = truth[key]
                    matched += 1
                elif flipped:
                    gt = gt[::-1]
                    carried += 1
            sample['GT'] = gt
            sample.phased = True
            out.write(r)
    return calls, matched, carried


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    n, m, c = transfer(*sys.argv[1:])
    print(f'heterozygous calls {n}; given the truth phase {m}; flipped with the preceding match {c}',
          file=sys.stderr)
