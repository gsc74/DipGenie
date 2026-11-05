import argparse
import json
from pathlib import Path
import tempfile
import unittest
import audit_inputs as a
import regions as r


class AuditTests(unittest.TestCase):
    def test_nested_files_and_contact_failure(self):
        with tempfile.TemporaryDirectory(prefix='hic-nesting-audit-') as temp:
            root = Path(temp)
            source = root / 'source'
            source.mkdir()
            # Small artificial region exercises every nested coverage without
            # allocating a chromosome-scale toy. The actual pipeline is not changed.
            original = r.REGIONS['MHC']
            r.REGIONS['MHC'] = ('chr6', 100, 200)
            try:
                import hashlib
                with (source / 'pairs.tsv').open('w') as meta, (source / 'R1.fastq').open('w') as one, (source / 'R2.fastq').open('w') as two:
                    meta.write('readID\tinside_bases\traw_bases\trank\n')
                    for n in range(80):
                        name = f'pair{n}'
                        rank = hashlib.blake2b(f'42:{name}'.encode(), digest_size=16).hexdigest()
                        meta.write(f'{name}\t50\t50\t{rank}\n')
                        one.write(f'@{name}/1\n' + 'A'*25 + '\n+\n' + 'I'*25 + '\n')
                        two.write(f'@{name}/2\n' + 'C'*25 + '\n+\n' + 'I'*25 + '\n')
                import gzip
                header = '#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2 read_side1 read_side2\n'
                with gzip.open(source / 'duplicates.pairs.gz', 'wt') as out:
                    out.write(header)
                with gzip.open(source / 'nodups.pairs.gz', 'wt') as out:
                    out.write(header)
                    for n in range(80):
                        out.write(f'pair{n}\tchr6\t111\tchr6\t190\t+\t-\tUU\t60\t60\t1\t2\n')
                r.sample(argparse.Namespace(input=source, output=root/'reads', region='MHC', seed=42))
                result = a.audit(root/'reads', 'MHC')
                self.assertEqual(result['status'], 'PASS')
                self.assertEqual(len(result['rows']), 6)
                self.assertEqual(result['rows'][-1]['pairs'], 80)
                r.sample(argparse.Namespace(input=source, output=root/'repeat', region='MHC', seed=42))
                r.sample(argparse.Namespace(input=source, output=root/'reseed', region='MHC', seed=43))
                self.assertEqual((root/'reads/membership.tsv').read_bytes(),
                                 (root/'repeat/membership.tsv').read_bytes())
                self.assertNotEqual((root/'reads/membership.tsv').read_bytes(),
                                    (root/'reseed/membership.tsv').read_bytes())
                self.assertEqual(a.audit(root/'reseed', 'MHC')['seed'], 43)
                for mate in (1,2):
                    self.assertEqual(list(r.fastq(root/f'reads/Total/HG002.R{mate}.fastq.gz')),
                                     list(r.fastq(root/f'reseed/Total/HG002.R{mate}.fastq.gz')))
                # The independent audit must detect a misleading seed receipt.
                receipt=root/'reseed/PASS.json'
                contents=json.loads(receipt.read_text());contents['seed']=42
                receipt.write_text(json.dumps(contents))
                with self.assertRaises(AssertionError): a.audit(root/'reseed', 'MHC')
                # Payload changes must fail even when IDs/nesting are preserved.
                path=root/'reads/4x/HG002.R1.fastq.gz'
                with gzip.open(path,'rt') as stream: original_payload=stream.read()
                with gzip.open(path,'wt') as stream:stream.write(original_payload.replace('A'*25,'G'+'A'*24,1))
                with self.assertRaises(AssertionError):a.audit(root/'reads', 'MHC')
                with gzip.open(path,'wt') as stream:stream.write(original_payload)
                with gzip.open(root/'reads/2x/contacts.pairs.gz', 'at') as out:
                    out.write('outside\tchr6\t111\tchr6\t190\t+\t-\tUU\t60\t60\t1\t2\n')
                with self.assertRaises(AssertionError):
                    a.audit(root/'reads', 'MHC')
            finally:
                r.REGIONS['MHC'] = original


if __name__ == '__main__':
    unittest.main()
