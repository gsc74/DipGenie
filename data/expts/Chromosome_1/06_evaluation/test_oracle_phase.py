#!/usr/bin/env python3
"""oracle_phase.py changes only the orientation of matching heterozygous calls (needs pysam)."""
from pathlib import Path
import tempfile
import unittest

import pysam

import oracle_phase

HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n')


def write(path, rows):
    path.write_text(HEADER + ''.join(f'chr1\t{p}\t.\t{r}\t{a}\t.\t.\t.\tGT\t{g}\n' for p, r, a, g in rows))


class Oracle(unittest.TestCase):
    def test_transfer(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write(d/'truth.vcf', [(10, 'A', 'C', '1|0'), (20, 'G', 'T', '0|1'), (30, 'C', 'G', '1|1'),
                                  (40, 'T', 'A', '.|1'), (60, 'C', 'T', '0|1')])
            write(d/'query.vcf', [(10, 'A', 'C', '0|1'), (20, 'G', 'A', '1|0'), (30, 'C', 'G', '0/1'),
                                  (40, 'T', 'A', '1/0'), (50, 'A', 'G', '1|1'), (60, 'C', 'T', '0|1'),
                                  (70, 'G', 'C', '1|0')])
            n, m, c = oracle_phase.transfer(str(d/'truth.vcf'), str(d/'query.vcf'), str(d/'out.vcf'))
            self.assertEqual((n, m, c), (6, 3, 2))
            got = [(r.pos, r.alts[0], r.samples[0]['GT'], r.samples[0].phased) for r in pysam.VariantFile(str(d/'out.vcf'))]
            self.assertEqual(got, [(10, 'C', (1, 0), True),   # truth orientation (flipped)
                                   (20, 'A', (0, 1), True),   # unmatched: follows the flip at 10
                                   (30, 'G', (1, 0), True),   # truth homozygous: follows the flip at 10
                                   (40, 'A', (0, 1), True),   # truth .|1 read as 0|1 (flipped)
                                   (50, 'G', (1, 1), True),   # homozygous call unchanged
                                   (60, 'T', (0, 1), True),   # matched, not flipped
                                   (70, 'C', (1, 0), True)])  # unmatched: keeps its phase


if __name__ == '__main__':
    unittest.main()
