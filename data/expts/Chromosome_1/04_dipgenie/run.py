#!/usr/bin/env python3
"""DipGenie v2 inference and variant calling for one manifest row; inputs supplied by a portable TSV.
The manifest's contacts column is not used."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

DEPTHS = ('2x', '4x', '8x', '16x', '32x', 'Total')
FIELDS = ('depth', 'sr_r1', 'sr_r2', 'hic_r1', 'hic_r2', 'contacts')

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''): h.update(block)
    return h.hexdigest()

def input_fields(mode):
    if mode not in ('SR', 'HiC', 'mixed'): raise ValueError('Unknown read mode')
    return ((('sr_r1', 'sr_r2') if mode in ('SR', 'mixed') else ()) +
            (('hic_r1', 'hic_r2') if mode in ('HiC', 'mixed') else ()))

def rows(path, mode='mixed'):
    path = Path(path).resolve()
    with path.open(newline='') as stream:
        table = csv.DictReader(stream, delimiter='\t')
        if tuple(table.fieldnames or ()) != FIELDS: raise ValueError('Invalid input TSV header')
        result = list(table)
    if tuple(row['depth'] for row in result) != DEPTHS:
        raise ValueError('Expected six unique rows in order: 2x 4x 8x 16x 32x Total')
    selected = input_fields(mode)
    for row in result:
        if None in row: raise ValueError('Extra input TSV columns')
        for field in selected:
            if not row.get(field) or row[field].strip() in ('', '.'):
                raise ValueError(f'Missing {field} for {mode}')
            value = Path(row[field])
            if not value.is_absolute(): value = path.parent / value
            value = value.resolve(strict=True)
            if not value.is_file() or value.stat().st_size == 0: raise ValueError(f'Missing/empty {field}')
            row[field] = str(value)
        mates = [row[f] for f in selected]
        if len(set(mates)) != len(mates): raise ValueError('A mate FASTQ is supplied more than once')
    return result

def verify_inputs(receipt, binary, index):
    if digest(binary) != receipt['binary_sha256']: raise RuntimeError('Binary changed during the run')
    if digest(index/'manifest.dgi') != receipt['index_manifest_sha256']:
        raise RuntimeError('Graph manifest changed during the run')
    if digest(__file__) != receipt['source_sha256']: raise RuntimeError('Runner changed during the run')
    for path, expected in receipt['input_sha256'].items():
        if digest(path) != expected: raise RuntimeError(f'Input changed during the run: {path}')

def physical_cpus(count):
    allowed = os.sched_getaffinity(0)
    table = subprocess.check_output(['lscpu', '-p=CPU,CORE,SOCKET,NODE,ONLINE'], text=True)
    nodes = {}
    for line in table.splitlines():
        if line.startswith('#'): continue
        cpu, core, socket, node, online = line.split(',')
        cpu, core, socket, node = int(cpu), int(core), int(socket), int(node)
        if online == 'Y' and cpu in allowed:
            nodes.setdefault(node, {}).setdefault((socket, core), cpu)
    for node in sorted(nodes):
        if node >= 0 and len(nodes[node]) >= count:
            return sorted(nodes[node].values())[:count]
    raise RuntimeError(f'Need {count} allocated physical cores on one NUMA node; refusing cross-node or SMT fallback')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cpu-list', action='store_true')
    p.add_argument('--binary', type=Path); p.add_argument('--gindex', type=Path)
    p.add_argument('--inputs', type=Path); p.add_argument('--output', type=Path)
    p.add_argument('--task', type=int, choices=range(6)); p.add_argument('--mode', choices=('SR', 'HiC', 'mixed'), default='mixed')
    p.add_argument('--threads', type=int, default=64); p.add_argument('--budget', type=int, default=2140)
    p.add_argument('--sample', default='HG002')
    a = p.parse_args()
    if a.threads < 1 or a.budget < 0: p.error('Invalid compute/budget setting')
    if a.cpu_list: print(','.join(map(str, physical_cpus(a.threads)))); return
    if any(getattr(a, k) is None for k in ('binary','gindex','inputs','output','task')):
        p.error('binary, gindex, inputs, output and task are required')
    row = rows(a.inputs, a.mode)[a.task]; binary = a.binary.resolve(strict=True); index = a.gindex.resolve(strict=True)
    destination = a.output.resolve() / a.mode / row['depth']; destination.mkdir(parents=True, exist_ok=False)
    command = [str(binary), '--gindex', str(index), '--out', str(destination/'inference'),
               '-t', str(a.threads), '-R', str(a.budget), '--sample', a.sample,
               '--output-vcf', str(destination/'calls.vcf.gz')]
    if a.mode in ('SR', 'mixed'): command += ['--sr','-r1',row['sr_r1'],'-r2',row['sr_r2']]
    if a.mode in ('HiC', 'mixed'): command += ['--hic','-h1',row['hic_r1'],'-h2',row['hic_r2']]
    receipt = dict(depth=row['depth'], mode=a.mode, command=command, binary_sha256=digest(binary),
                   index_manifest_sha256=digest(index/'manifest.dgi'),
                   input_sha256={str(row[k]):digest(row[k]) for k in input_fields(a.mode)},
                   source_sha256=digest(__file__), cpus=sorted(os.sched_getaffinity(0)))
    (destination/'inputs.json').write_text(json.dumps(receipt, indent=2)+'\n')
    with (destination/'stdout.log').open('w') as out, (destination/'inference.log').open('w') as err:
        result = subprocess.run(['/usr/bin/time','-v','-o',str(destination/'resources.txt'),*command], stdout=out, stderr=err)
    (destination/'exit-code.txt').write_text(str(result.returncode)+'\n')
    if result.returncode: raise SystemExit(result.returncode)
    for name in ('inference/all-components.fa','inference/run.tsv','calls.vcf.gz','calls.vcf.gz.csi'):
        if not (destination/name).is_file() or not (destination/name).stat().st_size:
            raise RuntimeError(f'Missing/empty output {name}')
    verify_inputs(receipt, binary, index)
    (destination/'COMPLETE').write_text('Inference and variant calling complete; accuracy evaluation is separate.\n')

if __name__ == '__main__': main()
