#!/usr/bin/env python3
"""Ensure target aliases/private ALTs are removed without dropping other IDs."""
import os
from pathlib import Path
import subprocess
import sys
import pysam
from audit_panel import is_target_sample, samples, verify

directory = Path(sys.argv[1])
directory.mkdir()
assert all(map(is_target_sample, ['HG002', 'HG002#1', 'HG002_2', 'HG002maternal',
                                 'NA24385', 'GM24385.1', 'asm.HG002.paternal']))
assert not any(map(is_target_sample, ['HG00232', 'HG00235', 'NA243850']))
header = pysam.VariantHeader()
header.add_meta('fileformat', value='VCFv4.2')
header.contigs.add('chr1', length=100)
header.formats.add('GT', 1, 'String', 'Genotype')
header.info.add('AC', 'A', 'Integer', 'Allele count')
header.info.add('AN', 1, 'Integer', 'Allele number')
for name in ['HG002', 'HG002#1', 'NA24385', 'GM24385_2', 'HG00232', 'CONTROL']:
    header.add_sample(name)
with pysam.VariantFile(str(directory / 'input.vcf'), 'w', header=header) as output:
    for pos, alleles, retained in [(10, ('A','C','G'), (2,2)),
                                  (20, ('A','C'), (0,0)),
                                  (30, ('A','C'), (1,1))]:
        record = output.new_record(contig='chr1', start=pos-1, alleles=alleles)
        record.info['AC'] = tuple(999 for _ in alleles[1:])  # stale on purpose
        record.info['AN'] = 999
        for name in header.samples:
            record.samples[name]['GT'] = (1,1) if is_target_sample(name) else retained
            record.samples[name].phased = True
        output.write(record)
samples(directory / 'input.vcf', directory)
subprocess.run(['bash', str(Path(__file__).with_name('filter_panel.sh')),
                str(directory / 'input.vcf'), str(directory / 'filtered.vcf.gz'),
                str(directory / 'keep.samples.txt')], check=True)
verify(directory / 'filtered.vcf.gz', directory)
with pysam.VariantFile(str(directory / 'filtered.vcf.gz')) as source:
    records = list(source)
    assert list(source.header.samples) == ['HG00232', 'CONTROL']
    assert [(r.pos, r.alts) for r in records] == [(10, ('G',)), (30, ('C',))]
    assert all(r.info['AC'] == (4,) and r.info['AN'] == 4 for r in records)
    assert all(s['GT'] == (1,1) and s.phased for r in records for s in r.samples.values())
print('PASS aliases/private ALT removal, GT remapping, shared ALT retention, stale AC repair')
