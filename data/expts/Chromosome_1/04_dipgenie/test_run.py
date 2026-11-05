#!/usr/bin/env python3
"""Portable single/mixed-library launcher guards; no scheduler or caller needed."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import run


class PairedRunTests(unittest.TestCase):
    def test_physical_core_selection_stays_on_one_numa_node(self):
        topology='# CPU,Core,Socket,Node,Online\n0,0,0,0,Y\n1,1,0,0,Y\n2,2,0,1,Y\n3,3,0,1,Y\n4,0,0,0,Y\n5,1,0,0,Y\n6,4,1,2,N\n'
        with patch.object(run.os,'sched_getaffinity',return_value=set(range(7))), patch.object(run.subprocess,'check_output',return_value=topology):
            self.assertEqual(run.physical_cpus(2),[0,1])
            with self.assertRaises(RuntimeError):run.physical_cpus(3)
        with patch.object(run.os,'sched_getaffinity',return_value={2,3,4}), patch.object(run.subprocess,'check_output',return_value=topology):
            self.assertEqual(run.physical_cpus(2),[2,3])

    def fixture(self, root, mode):
        selected = run.input_fields(mode)
        for name in selected:
            (root/name).write_bytes(b'toy input\n')
        table = root/'reads.tsv'
        with table.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=run.FIELDS,
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for depth in run.DEPTHS:
                writer.writerow(dict(depth=depth, **{
                    name: name if name in selected else '.' for name in run.FIELDS[1:]}))
        return table

    def test_each_mode_requires_only_its_own_inputs(self):
        for mode, count in (('SR', 2), ('HiC', 2), ('mixed', 4)):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); table = self.fixture(root, mode)
                rows = run.rows(table, mode)
                self.assertEqual([r['depth'] for r in rows], list(run.DEPTHS))
                self.assertEqual(len(run.input_fields(mode)), count)
                self.assertTrue(all(Path(rows[0][f]).is_absolute() for f in run.input_fields(mode)))
                for field in set(run.FIELDS[1:]) - set(run.input_fields(mode)):
                    self.assertEqual(rows[0][field], '.')
                if mode != 'mixed':
                    with self.assertRaises(ValueError): run.rows(table, 'mixed')

    def test_selected_missing_empty_and_aliased_mate(self):
        for problem in ('missing', 'empty', 'alias'):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); table = self.fixture(root, 'SR')
                mate = root/'sr_r2'
                if problem == 'missing': mate.unlink()
                elif problem == 'empty': mate.write_bytes(b'')
                else: mate.unlink(); mate.symlink_to(root/'sr_r1')
                with self.assertRaises((ValueError, FileNotFoundError)): run.rows(table, 'SR')

    def test_contacts_column_is_ignored(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);table=self.fixture(root,'mixed')
            self.assertNotIn('contacts',run.input_fields('mixed'))
            self.assertEqual(run.rows(table,'mixed')[0]['contacts'],'.')

    def test_depth_order_and_extra_columns(self):
        for problem in ('depth', 'columns'):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); table = self.fixture(root, 'SR')
                text = table.read_text()
                if problem == 'depth': text = text.replace('2x\t', '4x\t', 1)
                else: text = text.replace('\n4x\t', '\textra\n4x\t', 1)
                table.write_text(text)
                with self.assertRaises(ValueError): run.rows(table, 'SR')

    def test_post_run_integrity(self):
        for changed in ('none', 'binary', 'graph', 'input', 'runner'):
            with self.subTest(changed=changed), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); binary = root/'binary'; index = root/'index'
                index.mkdir(); manifest = index/'manifest.dgi'; read = root/'reads'
                for path in (binary, manifest, read): path.write_bytes(b'original')
                receipt = dict(binary_sha256=run.digest(binary),
                               index_manifest_sha256=run.digest(manifest),
                               source_sha256=run.digest(run.__file__),
                               input_sha256={str(read): run.digest(read)})
                if changed == 'runner': receipt['source_sha256'] = 'changed'
                elif changed != 'none':
                    {'binary': binary, 'graph': manifest, 'input': read}[changed].write_bytes(b'changed')
                if changed == 'none': run.verify_inputs(receipt, binary, index)
                else:
                    with self.assertRaises(RuntimeError): run.verify_inputs(receipt, binary, index)

    def test_actual_launcher_selects_only_requested_libraries(self):
        for mode in ('SR','HiC','mixed'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); table = self.fixture(root, mode)
                binary = root/'binary'; binary.write_bytes(b'toy executable')
                index = root/'index'; index.mkdir(); (index/'manifest.dgi').write_bytes(b'toy index')
                calls = []
                def mock_run(command, **kwargs):
                    calls.append(command)
                    inference = Path(command[command.index('--out')+1]); inference.mkdir()
                    for name in ('all-components.fa', 'run.tsv'):
                        (inference/name).write_bytes(b'toy output\n')
                    vcf = Path(command[command.index('--output-vcf')+1])
                    vcf.write_bytes(b'toy caller output\n')
                    Path(str(vcf)+'.csi').write_bytes(b'toy index\n')
                    return SimpleNamespace(returncode=0)
                args = ['run.py', '--binary', str(binary), '--gindex', str(index),
                        '--inputs', str(table), '--output', str(root/'out'),
                        '--task', '0', '--mode', mode]
                with patch('sys.argv', args), patch.object(run.subprocess, 'run', side_effect=mock_run):
                    run.main()
                command = calls[0]
                self.assertEqual('--sr' in command, mode in ('SR', 'mixed'))
                self.assertEqual('--hic' in command, mode in ('HiC', 'mixed'))
                options = {x for x in command[command.index('--gindex'):] if x.startswith('-')}
                self.assertLessEqual(options, {'--gindex', '--out', '-t', '-R', '--sample', '--output-vcf',
                                               '--sr', '-r1', '-r2', '--hic', '-h1', '-h2'})
                dest = root/'out'/mode/'2x'
                receipt = json.loads((dest/'inputs.json').read_text())
                self.assertEqual(set(receipt['input_sha256']),
                                 {str(root/name) for name in run.input_fields(mode)})
                self.assertTrue((dest/'COMPLETE').is_file())


if __name__ == '__main__': unittest.main()
