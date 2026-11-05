#!/usr/bin/env python3
"""Dipcall VCF contract: native normalization/CSI versus optional bcftools oracle.

The implementation always runs with empty PATH. bcftools is test-only.
"""
import argparse
import gzip
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('helper', type=Path)
p.add_argument('--bcftools', type=Path)
p.add_argument('--keep', type=Path)
a = p.parse_args()
helper = a.helper.resolve()
bcf = a.bcftools.resolve() if a.bcftools else None


def cmd(args, native=False, ok=True):
    r = subprocess.run(list(map(str, args)), env=dict(os.environ, PATH='') if native else None,
                       text=True, capture_output=True)
    assert (r.returncode == 0) == ok, (args, r.returncode, r.stderr)
    return r.stdout


def text(file):
    return gzip.open(file, 'rt').read() if str(file).endswith('.gz') else file.read_text()


def body(file):
    return [s for s in text(file).splitlines() if s and s[0] != '#']


def header(contigs, sample='HG002'):
    return ('##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID={n},length={len(s)}>\n' for n, s in contigs)
            + '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            + '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">\n'
            + '##FILTER=<ID=GAP1,Description="Gap">\n'
            + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n')


def run(root):
    rng = random.Random(8187)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    reference = dna(400)+'A'*1300+'AC'*800+dna(10000)
    contigs = [('chr2', reference), ('chr10', dna(5000)), ('empty', dna(100))]
    ref = root/'ref.fa'; ref.write_text(''.join('>'+n+'\n'+s+'\n' for n, s in contigs))
    records = []

    def add(name, pos, ref, alts, gt='0|1', filt='.', qual='30'):
        alts = list(dict.fromkeys(s for s in alts if s.upper() != ref.upper()))
        if not alts:
            return
        ad = ','.join(str(i+1) for i in range(len(alts)+1))
        records.append([name, str(pos), '.', ref, ','.join(alts), qual, filt, '.', 'GT:AD', gt+':'+ad])

    # Chromosome-start indels, repetitive shifts longer than the padding chunk,
    # substitutions with common context, multiallelic GT/AD, and missing calls.
    add('chr2', 1, reference[:2], [reference[:1]])
    add('chr2', 1, reference[0], [reference[0]*2])
    add('chr2', 1000, 'A', ['AA', 'AAA'], '1|2')
    add('chr2', 1100, 'AAA', ['A'])
    add('chr2', 2101, reference[2100:2110], [reference[2100:2104]])
    add('chr2', 3501, reference[3500:3505].lower(), ['a'+reference[3501:3505].lower()])
    # Dipcall emits '*' for sites overlapped by an upstream deletion. It is
    # not literal sequence and must not be trimmed or left-shifted. Include
    # mixed SNP/indel/star records, both phase orders and missing genotypes.
    add('chr2', 3700, reference[3699:3704], [reference[3699]])
    for i, gt in enumerate(('0|1', '1|0', '1|1', '.|1', '0|0')):
        add('chr2', 3701+i, reference[3700+i:3703+i], ['*'], gt)
    add('chr2', 3720, reference[3719:3722], [reference[3719], '*'], '1|2')
    add('chr2', 3740, reference[3739:3742], ['*', reference[3739]], '1|2')
    add('chr2', 3760, reference[3759], [next(c for c in 'ACGT' if c != reference[3759]), '*'], '2|1')
    add('chr2', 1001, 'AAA', ['*', 'A', 'AAAA'], '1|2')
    for name, seq in contigs[:2]:
        for j in range(700):
            pos = rng.randrange(1, len(seq)-100)
            n = rng.randint(1, 20); allele = seq[pos-1:pos+n-1]
            alts = [dna(rng.randint(1, 20)) for _ in range(rng.randint(1, 3))]
            alts = list(dict.fromkeys(s for s in alts if s != allele))
            if not alts:
                continue
            # Several short indels have common leading/trailing sequence.
            if j % 3 == 0:
                alts = [allele[0]+s+allele[-1] for s in alts]
            g = lambda: '.' if rng.random() < .04 else str(rng.randrange(len(alts)+1))
            add(name, pos, allele, alts, g()+'|'+g(), 'GAP1' if j % 17 == 0 else 'PASS' if j % 7 == 0 else '.')
    # Exact duplicates: retain first before the called/nonreference filter.
    add('chr2', 100, reference[99], ['C' if reference[99] != 'C' else 'G'], '0|0')
    add('chr2', 100, reference[99], ['C' if reference[99] != 'C' else 'G'], '1|1')
    order = {n: i for i, (n, _) in enumerate(contigs)}
    records.sort(key=lambda c: (order[c[0]], int(c[1])))
    source = root/'raw.vcf'
    source.write_text(header(contigs)+''.join('\t'.join(c)+'\n' for c in records))
    normalized = root/'normalized.vcf'
    cmd([helper, 'normalize', ref, source, normalized], native=True)
    star_rows = {(c[0], int(c[1]), c[3]): c[9]
                 for line in body(normalized) if (c := line.split('\t'))[4] == '*'}
    for pos, gt in ((3701, '0|1:1,2'), (3702, '1|0:1,2'), (3703, '1|1:1,2')):
        assert star_rows[('chr2', pos, reference[pos-1:pos+2])] == gt
    assert ('chr2', 3704, reference[3703:3706]) not in star_rows  # missing GT
    assert ('chr2', 3705, reference[3704:3707]) not in star_rows  # reference-only GT
    assert star_rows[('chr2', 3720, reference[3719:3722])] == '0|1:1,3'
    assert star_rows[('chr2', 3740, reference[3739:3742])] == '1|0:1,2'
    assert star_rows[('chr2', 3760, reference[3759])] == '1|0:1,3'
    assert star_rows[('chr2', 1001, 'AAA')] == '1|0:1,2'  # no repeat left shift
    compressed = root/'normalized.vcf.gz'
    cmd([helper, 'publish', normalized, compressed], native=True)
    assert body(normalized) == body(compressed)
    assert Path(str(compressed)+'.csi').is_file()
    cmd([helper, 'publish', normalized, compressed], native=True, ok=False)
    # Header-only and component dictionaries (including no-call contigs).
    empty = root/'empty.vcf'; empty.write_text(header([contigs[-1]]))
    cmd([helper, 'publish', empty, root/'empty.vcf.gz'], native=True)
    assert not body(root/'empty.vcf.gz')
    part = root/'part.vcf'; part.write_text(header([('another', 'A'*100)])+'another\t20\t.\tA\tC\t30\t.\t.\tGT:AD\t1|0:1,1\n')
    joined = root/'joined.vcf.gz'
    cmd([helper, 'combine', joined, empty, part, compressed], native=True)
    assert len(body(joined)) == len(body(compressed))+1
    wrong = root/'wrong.vcf'; wrong.write_text(header([contigs[-1]], 'other'))
    cmd([helper, 'combine', root/'wrong-joined.vcf', empty, wrong], native=True, ok=False)
    for label, row in (
        ('ref-mismatch', f'chr2\t1\t.\t{next(c for c in "ACGT" if c != reference[0])}\tN\t30\t.\t.\tGT:AD\t0|1:1,1'),
        ('unknown-annotation', f'chr2\t1\t.\t{reference[0]}\tN\t30\t.\tAC=1\tGT:AD\t0|1:1,1'),
        ('bad-genotype', f'chr2\t1\t.\t{reference[0]}\tN\t30\t.\t.\tGT:AD\t0|2:1,1'),
        ('star-ref-mismatch', f'chr2\t1\t.\t{next(c for c in "ACGT" if c != reference[0])}\t*\t30\t.\t.\tGT:AD\t0|1:1,1'),
        ('star-in-ref', 'chr2\t1\t.\t*\tA\t30\t.\t.\tGT:AD\t0|1:1,1'),
        ('star-in-sequence', f'chr2\t1\t.\t{reference[0]}\tA*\t30\t.\t.\tGT:AD\t0|1:1,1'),
        ('symbolic-alt', f'chr2\t1\t.\t{reference[0]}\t<DEL>\t30\t.\t.\tGT:AD\t0|1:1,1'),
    ):
        path = root/(label+'.vcf'); path.write_text(header(contigs)+row+'\n')
        cmd([helper, 'normalize', ref, path, root/(label+'.out.vcf')], native=True, ok=False)
    if bcf:
        # Same published normalization/filter order; uppercase before norm.
        upper = root/'upper.vcf'
        upper.write_text(header(contigs)+''.join('\t'.join(c[:3]+[c[3].upper(),c[4].upper()]+c[5:])+'\n' for c in records))
        cmd([bcf, 'view', '-f', '.,PASS', '-Ob', '-o', root/'pass.bcf', upper])
        cmd([bcf, 'norm', '-f', ref, '-m', '-any', '-Ob', '-o', root/'split.bcf', root/'pass.bcf'])
        cmd([bcf, 'norm', '-d', 'exact', '-Ob', '-o', root/'unique.bcf', root/'split.bcf'])
        cmd([bcf, 'view', '-i', 'GT="alt" && GT!="mis"', '-Ob', '-o', root/'called.bcf', root/'unique.bcf'])
        gold = root/'gold.vcf'; cmd([bcf, 'sort', '-Ov', '-o', gold, root/'called.bcf'])
        actual, expected = body(normalized), body(gold)
        assert actual == expected, ('normalization parity', len(actual), len(expected),
                                  next(((x,y) for x,y in zip(actual,expected) if x!=y), None))
        for file in (compressed, joined, root/'empty.vcf.gz'):
            assert cmd([bcf, 'query', '-l', file]).strip() == 'HG002'
            for name in ('chr2', 'chr10', 'empty', 'another'):
                observed = cmd([bcf, 'view', '-H', '-r', name, file]).splitlines()
                assert observed == [s for s in body(file) if s.startswith(name+'\t')], (file,name)
        # Partial interval queries must return overlapping deletions, not only
        # sites with POS in range. This independently checks CSI offsets/spans.
        for lo, hi in ((1,100),(500,505),(1000,1100),(2000,2010),(3500,3600)):
            observed = cmd([bcf,'view','-H','-r',f'chr2:{lo}-{hi}',compressed]).splitlines()
            expected = [s for s in body(compressed) if (c:=s.split('\t'))[0]=='chr2'
                        and int(c[1])<=hi and int(c[1])+len(c[3])-1>=lo]
            assert observed == expected, ('CSI overlap',lo,hi)
        print('PASS bcftools 1.23.1 differential normalization, multiallelic GT/AD, filters, duplicates and CSI queries')
    print('PASS empty-PATH normalization, compression, component joining, empty CSI and invalid-input guards')


if a.keep:
    a.keep.mkdir(parents=True, exist_ok=False);run(a.keep)
else:
    with tempfile.TemporaryDirectory(prefix='dipgenie-vcf-') as tmp:
        run(Path(tmp))
