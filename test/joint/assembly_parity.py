#!/usr/bin/env python3
"""Fresh assembly fixtures; optional differential check against frozen dipcall.

The native helper runs with PATH empty. --kit supplies the independent dipcall
0.3 kit (minimap2/HTSbox) oracle, not a production dependency. No biological
accuracy claim is made from these caller-contract fixtures.
"""
import argparse
import os
from pathlib import Path
import random
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('helper', type=Path)
p.add_argument('--kit', type=Path)
p.add_argument('--fasta-helper', type=Path)
p.add_argument('--keep', type=Path)
a = p.parse_args()
helper = a.helper.resolve()
kit = a.kit.resolve() if a.kit else None


def command(args, *, native=False, ok=True, output=None):
    env = dict(os.environ)
    if native:
        env['PATH'] = ''
    r = subprocess.run(list(map(str, args)), env=env, capture_output=True)
    assert (r.returncode == 0) == ok, (args, r.returncode, r.stderr.decode())
    if output:
        output.write_bytes(r.stdout)
    return r.stdout.decode()


def fasta(path, records):
    path.write_text(''.join('>'+name+'\n'+s+'\n' for name, s in records))


def run(root):
    rng = random.Random(9417)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    ref = dna(320000)
    alt = lambda x: next(c for c in 'ACGT' if c != x)
    rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

    def edit(changes):
        result = ''; pos = 0
        for start, end, value in sorted(changes):
            assert pos <= start <= end
            result += ref[pos:start]+value; pos = end
        return result+ref[pos:]

    h1 = edit([(11000, 11001, alt(ref[11000])), (40000, 40001, alt(ref[40000])),
               (90000, 90120, ''), (140000, 140000, dna(75)),
               (190000, 260000, rc(ref[190000:260000]))])
    h2 = edit([(23000, 23001, alt(ref[23000])), (40000, 40001, alt(ref[40000])),
               (120000, 120000, dna(81)), (220000, 220065, '')])
    ref_file = root/'reference.fa'
    fasta(ref_file, [('chr1', ref)])
    hap_files = [root/'h1.fa', root/'h2.fa']
    for h, seq in enumerate((h1, h2)):
        fasta(hap_files[h], [(f'hap{h+1}', seq), (f'short{h+1}', ref[70000:72000]),
                            (f'none{h+1}', 'N'*2000)])

    prefixes = [root/'one', root/'two']
    for threads, prefix in zip((1, 2), prefixes):
        stats = command([helper, 'align', ref_file, *hap_files, prefix, threads], native=True)
        assert all(int(row.split('\t')[1]) == 3 for row in stats.splitlines())
        command([helper, 'pileup', ref_file, *(str(prefix)+f'.hap{h}.bam' for h in (1, 2)),
                 str(prefix)+'.pair.vcf'], native=True)
        command([helper, 'pair', str(prefix)+'.pair.vcf', str(prefix)+'.dip.vcf', 'HG002'], native=True)
        command([helper, 'bed', *(str(prefix)+f'.hap{h}.paf' for h in (1, 2)), prefix], native=True)
    for ext in ('.hap1.paf', '.hap2.paf', '.hap1.bam', '.hap2.bam', '.hap1.bed', '.hap2.bed', '.dip.bed'):
        assert Path(str(prefixes[0])+ext).read_bytes() == Path(str(prefixes[1])+ext).read_bytes(), ext
    for ext in ('.pair.vcf', '.dip.vcf'):
        def body(prefix):
            return [s for s in Path(str(prefix)+ext).read_text().splitlines() if not s.startswith('#')]
        assert body(prefixes[0]) == body(prefixes[1]), ext
    called = {int(s.split('\t')[1]): s.split('\t') for s in body(prefixes[0])}
    assert called[11001][9].split(':')[0] == '1|0'
    assert called[23001][9].split(':')[0] == '0|1'
    assert called[40001][9].split(':')[0] == '1|1'
    assert any(abs(len(s[3])-len(s[4])) >= 50 for s in called.values())
    if a.fasta_helper:
        fasta(root/'diploid.fa', [('hap1',h1),('hap2',h2)])
        command([a.fasta_helper.resolve(),root/'diploid.fa',ref_file,root/'helper.vcf.gz'],native=True)
        expected=root/'expected.vcf'
        command([helper,'normalize',ref_file,str(prefixes[0])+'.dip.vcf',expected],native=True)
        import gzip
        with gzip.open(root/'helper.vcf.gz','rt') as f:
            observed=[s.rstrip('\n') for s in f if not s.startswith('#')]
        assert observed==[s for s in expected.read_text().splitlines() if not s.startswith('#')]
        assert (root/'helper.vcf.gz.csi').is_file()
        print('PASS VG FASTA helper uses the same embedded caller with empty PATH')
    command([helper, 'align', ref_file, *hap_files, prefixes[0], 1], native=True, ok=False)
    command([helper, 'align', ref_file, *hap_files, root/'bad-threads', 0], native=True, ok=False)
    fasta(root/'empty.fa', [])
    command([helper, 'align', ref_file, root/'empty.fa', hap_files[1], root/'empty', 1], native=True, ok=False)

    # Exhaustive diploid/no-PAR transformation on all allele/missing states,
    # including haploid assembly HET filters, missing mate, and multiallelic AD.
    variants = root/'pair-fixture.vcf'
    header = '##fileformat=VCFv4.2\n##source=fixture\n##reference=ref\n'
    header += '##contig=<ID=unit,length=10000>\n'
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    header += '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Depth">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thap1\thap2\n'
    rows = []
    states = ('0/0', '1/1', '2/2', '0/1', '1/2', './.', './1', '1/.')
    for g1 in states:
        for g2 in states:
            rows.append(f'unit\t{len(rows)+1}\t.\tA\tC,G\t99\t.\t.\tGT:AD\t{g1}:1,2,3\t{g2}:4,5,6\n')
    rows.append('unit\t999\t.\tA\tN\t99\t.\t.\tGT:AD\t0/0:1,0\t1/1:0,1\n')
    variants.write_text(header+''.join(rows))
    transformed = root/'pair-fixture.out.vcf'
    command([helper, 'pair', variants, transformed, 'HG002'], native=True)
    assert len([s for s in transformed.read_text().splitlines() if not s.startswith('#')]) == 64
    invalid = root/'invalid.vcf'; invalid.write_text(header+rows[0].replace('GT:AD', 'GT:DP'))
    command([helper, 'pair', invalid, root/'invalid-out.vcf', 'HG002'], native=True, ok=False)

    if kit:
        golden = command([kit/'k8', kit/'dipcall-aux.js', 'vcfpair', '-a', '-s', 'HG002', variants])
        assert transformed.read_text() == golden, 'vcfpair corner-case parity'
        golden_bams = []
        for h, query in enumerate(hap_files, 1):
            paf = command([kit/'minimap2', '-c', '--paf-no-hit', '-xasm5', '--cs', '-t2', ref_file, query])
            def no_rl(text):
                return ['\t'.join(c for c in line.split('\t') if not c.startswith('rl:i:'))
                        for line in text.splitlines() if line and not line.startswith('@')]
            assert no_rl(paf) == no_rl(Path(str(prefixes[0])+f'.hap{h}.paf').read_text()), f'PAF {h}'
            sam = root/f'gold{h}.sam'
            command([kit/'minimap2', '-a', '-xasm5', '--cs', '-t2', ref_file, query], output=sam)
            filtered = root/f'gold{h}.filtered.sam'
            command([kit/'k8', kit/'dipcall-aux.js', 'samflt', sam], output=filtered)
            bam = root/f'gold{h}.bam'
            command([kit/'samtools', 'sort', '-o', bam, filtered])
            golden_bams.append(bam)
            golden_sam = command([kit/'samtools', 'view', bam])
            native_sam = command([kit/'samtools', 'view', str(prefixes[0])+f'.hap{h}.bam'])
            assert no_rl(golden_sam) == no_rl(native_sam), f'sorted BAM {h}'
            sorted_paf = root/f'gold{h}.sorted.paf'
            lines = sorted(paf.splitlines(), key=lambda s: (s.split('\t')[5], int(s.split('\t')[7]), s))
            sorted_paf.write_text('\n'.join(lines)+'\n')
            calls = command([kit/'k8', kit/'paftools.js', 'call', sorted_paf])
            bed = ''.join(s[2:]+'\n' for s in calls.splitlines() if s.startswith('R\t'))
            assert bed == Path(str(prefixes[0])+f'.hap{h}.bed').read_text(), f'callable hap BED {h}'
        bed = command([kit/'bedtk', 'isec', '-m', *(str(prefixes[0])+f'.hap{h}.bed' for h in (1, 2))])
        assert bed == Path(str(prefixes[0])+'.dip.bed').read_text(), 'diploid callable BED'
        pair = root/'gold.pair.vcf'
        command([kit/'htsbox', 'pileup', '-q5', '-evcf', ref_file, *golden_bams], output=pair)
        golden_vcf = command([kit/'k8', kit/'dipcall-aux.js', 'vcfpair', '-a', '-s', 'HG002', pair])
        assert [s for s in golden_vcf.splitlines() if not s.startswith('#')] == body(prefixes[0]), 'raw diploid VCF parity'
        print('PASS frozen dipcall parity: all PAFs and sorted BAMs (excluding unused rl tag), raw diploid calls, 64 genotype/filter cases')
    print('PASS native assembly API: SNP/indel/inversion, unmapped/short contigs, 1/2-thread determinism, empty PATH and failure guards')


if a.keep:
    a.keep.mkdir(parents=True, exist_ok=False)
    run(a.keep)
else:
    with tempfile.TemporaryDirectory(prefix='dipgenie-assembly-') as tmp:
        run(Path(tmp))
