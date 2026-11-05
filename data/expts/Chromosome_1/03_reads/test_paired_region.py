#!/usr/bin/env python3
"""Small deterministic checks before the expensive whole-reference alignment."""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile

import pysam
import paired_region as module


def record(name, mate, pos, cigar='10M', reverse=False, unmapped=False, extra=0):
    flag = 1 | (64 if mate == 1 else 128) | (16 if reverse else 0) | (4 if unmapped else 0) | extra
    return '\t'.join(map(str, (name, flag, '*' if unmapped else 'chr6', 0 if unmapped else pos,
        60, '*' if unmapped else cigar, '=', 1, 0, 'AACCGGTTAA', 'ABCDEFGHIJ')))


def run(root):
    sam = root/'toy.sam'
    groups = [
        [record('inside', 1, 101), record('inside', 2, 500)],
        [record('outside', 1, 50), record('outside', 2, 500)],
        [record('reverse', 2, 101, reverse=True), record('reverse', 1, 500)],
        [record('mate-unmapped', 1, 101), record('mate-unmapped', 2, 0, unmapped=True)],
        [record('supp-only', 1, 500), record('supp-only', 1, 101, extra=2048), record('supp-only', 2, 600)],
        [record('secondary-only', 1, 500), record('secondary-only', 1, 101, extra=256), record('secondary-only', 2, 600)],
        [record('deletion-only', 1, 86, '5M20D5M'), record('deletion-only', 2, 600)],
        [record('skip-only', 1, 86, '5M20N5M'), record('skip-only', 2, 600)],
        [record('left-touch', 1, 91), record('left-touch', 2, 600)],
        [record('right-touch', 1, 111), record('right-touch', 2, 600)],
        [record('one-base', 1, 92), record('one-base', 2, 600)],
        [record('softclip', 1, 101, '5S5M'), record('softclip', 2, 600)],
    ]
    sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr6\tLN:1000\n'+'\n'.join(r for g in groups for r in g)+'\n')
    selected = root/'selected'
    module.extract(argparse.Namespace(sam=str(sam), output=str(selected), region='chr6:101-110', accession='toy'))
    pairs = list(module.pairs(selected))
    assert len(pairs) == 5
    names = {p[0][0].removesuffix('/1') for p in pairs}
    assert names == {'@inside','@reverse','@mate-unmapped','@one-base','@softclip'}, names
    reverse = next(p for p in pairs if p[0][0] == '@reverse/1')[1]
    assert reverse[1] == 'TTAACCGGTT' and reverse[2] == 'JIHGFEDCBA'
    assert json.loads((selected/'selection.json').read_text())['selected_unmapped_mates'] == 1
    with pysam.AlignmentFile(str(selected/'selected.primary.bam')) as bam:
        assert sum(1 for _ in bam) == 10
    for length, expected in ((10,['2x','4x','8x','Total']), (13,['2x','4x','Total']), (25,['2x','Total']), (50,['Total'])):
        sampled = root/f'depths-{length}'
        module.sample(argparse.Namespace(input=str(selected), output=str(sampled), length=length, seed=42,
            accession='toy', reference='toy.fa', source_identity='toy'))
        actual = [r.split('\t')[0] for r in (sampled/'available_depths.tsv').read_text().splitlines()[1:]]
        assert actual == expected, (length, actual)
        previous = set()
        for depth in actual:
            content = set(module.pairs(sampled/depth))
            assert previous <= content
            previous = content
        assert previous == set(pairs)
    bad = root/'bad.sam'
    bad.write_text('@SQ\tSN:chr6\tLN:1000\n'+record('missing-mate', 1, 101)+'\n')
    try:
        module.extract(argparse.Namespace(sam=str(bad), output=str(root/'bad'), region='chr6:101-110', accession='toy'))
    except ValueError:
        pass
    else:
        raise AssertionError('Missing mate was silently accepted')
    print('PASS: primary-only selection, M/=X overlap, reverse orientation, unmapped mate, boundaries, soft clips, pair order, strict-depth eligibility and nested identity')


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='sr-wgs-toy-') as temp:
        run(pathlib.Path(temp))
