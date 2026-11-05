#!/usr/bin/env python3
"""Protocol guards that do not require the external callers/evaluators."""
import json
from pathlib import Path
import tempfile
import unittest
import evaluate


class Protocol(unittest.TestCase):
    def test_commands(self):
        bench, refine = evaluate.bench_commands('truth', 'query', 'ref', 'GIAB.bed', 'out', 8)
        self.assertEqual(bench[bench.index('--includebed')+1], 'GIAB.bed')
        self.assertEqual(bench[bench.index('--pick')+1], 'ac')
        self.assertEqual(bench[bench.index('-r')+1], '2000')
        self.assertIn('--passonly', bench)
        self.assertIn('--use-original-vcfs', refine)
        self.assertEqual(refine[refine.index('--align')+1], 'mafft')
        for cmd in (bench, refine):
            self.assertNotIn('--subset', cmd)
            self.assertNotIn('--write-phab', cmd)

    def test_refined_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            before = {'base cnt': 10, 'comp cnt': 12, 'TP-base': 8, 'TP-comp': 8, 'FN': 2, 'FP': 4}
            after = dict(before, **{'TP-base': 9, 'TP-comp': 11, 'FN': 1, 'FP': 1})
            (p/'summary.json').write_text(json.dumps(before))
            (p/'refine.variant_summary.json').write_text(json.dumps(after))
            (p/'refine.log.txt').write_text('Finished refine\n')
            self.assertEqual(evaluate.refined_summary(p)['refined'], after)
            after['base cnt'] = 9
            (p/'refine.variant_summary.json').write_text(json.dumps(after))
            with self.assertRaises(RuntimeError):
                evaluate.refined_summary(p)
            (p/'refine.log.txt').write_text('0 regions to be refined\nNo regions to be refined. Exiting\n')
            self.assertEqual(evaluate.refined_summary(p)['refined'], before)
            (p/'refine.log.txt').write_text('[ERROR] alignment failed\nFinished refine\n')
            with self.assertRaises(RuntimeError):
                evaluate.refined_summary(p)


if __name__ == '__main__':
    unittest.main()
