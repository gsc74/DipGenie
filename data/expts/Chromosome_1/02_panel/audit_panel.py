#!/usr/bin/env python3
"""Explicit sample alias exclusion and retained-genotype ALT validation."""
import collections
import json
from pathlib import Path
import re
import sys
import pysam


def is_target_sample(name):
    # HG00232 is a distinct registered sample, not an HG002 haplotype.
    # Accept no literal HG002/NA24385/GM24385 alias, including assembly,
    # maternal/paternal, numbered-haplotype and delimiter-qualified forms.
    return re.search(r'(?:HG002|NA24385|GM24385)(?![0-9])', name, re.I) is not None


def samples(vcf, directory):
    with pysam.VariantFile(str(vcf)) as source:
        names = list(source.header.samples)
    assert len(names) == len(set(names)) and names
    removed = [s for s in names if is_target_sample(s)]
    kept = [s for s in names if not is_target_sample(s)]
    assert kept
    for filename, values in [('keep.samples.txt', kept), ('removed.samples.txt', removed)]:
        (directory / filename).write_text(''.join(s + '\n' for s in values))
    (directory / 'sample-audit.json').write_text(json.dumps(dict(
        input_samples=len(names), retained_samples=len(kept), excluded_samples=removed,
        excluded_identity_aliases=['HG002', 'NA24385', 'GM24385'],
        separate_samples_with_longer_numeric_ids=[s for s in kept if s.startswith('HG002')]), indent=2) + '\n')
    print(f'PASS sample audit: input={len(names)} retained={len(kept)} removed={removed}')


def verify(vcf, directory):
    expected = (directory / 'keep.samples.txt').read_text().splitlines()
    records = alleles = 0
    with pysam.VariantFile(str(vcf), threads=4) as source:
        assert list(source.header.samples) == expected and not any(map(is_target_sample, expected))
        assert list(source.header.contigs) == ['chr1']
        for record in source:
            assert record.contig == 'chr1' and record.alts, (record.contig, record.pos)
            # AC/AN have been recomputed from retained GTs by fill-tags;
            # verify every ALT has a carrier and the haploid denominator is sane.
            ac = record.info['AC']
            an = record.info['AN']
            assert len(ac) == len(record.alts) and all(v is not None and v > 0 for v in ac)
            assert 0 < an <= 2 * len(expected) and sum(ac) <= an
            records += 1
            alleles += len(record.alts)
    assert records > 0
    result = dict(records=records, alternate_alleles=alleles, samples=len(expected),
                  contigs=['chr1'], HG002_aliases_present=False,
                  unsupported_ALT_alleles=0, AC_AN_source='bcftools fill-tags from retained genotypes',
                  shared_alleles_carried_by_remaining_samples='retained')
    (directory / 'panel-audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    command, vcf, directory = sys.argv[1:]
    (samples if command == 'samples' else verify)(Path(vcf), Path(directory))
