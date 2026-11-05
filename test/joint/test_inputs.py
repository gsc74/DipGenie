#!/usr/bin/env python3
"""Fresh paired inputs: GFA/VCF/index, component routing, invalid input guards."""
import csv
import gzip
from pathlib import Path
import random
import os
import shutil
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
def run(*args, ok=True):
    p = subprocess.run([binary, *map(str, args)], text=True, capture_output=True)
    if ok and p.returncode and '--out' in args:
        for log in Path(args[args.index('--out')+1]).rglob('caller.log'):
            print(log.read_text(), file=sys.stderr)
    assert (p.returncode == 0) == ok, p.stdout + p.stderr
    return p

with tempfile.TemporaryDirectory(prefix='dipgenie-inputs-') as tmp:
    root = Path(tmp)
    rng = random.Random(122)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    seq = [dna(3000), dna(3000)]
    (root/'reference.fa').write_text(''.join(f'>chr{i+1}\n{s}\n' for i,s in enumerate(seq)))
    vcf = '##fileformat=VCFv4.2\n' + ''.join(f'##contig=<ID=chr{i+1},length=3000>\n' for i in range(2))
    vcf += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    vcf += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpanel\n'
    for i,s in enumerate(seq):
        alt = next(b for b in 'ACGT' if b != s[1500])
        vcf += f'chr{i+1}\t1501\tv{i}\t{s[1500]}\t{alt}\t.\tPASS\t.\tGT\t0|1\n'
    (root/'panel.vcf').write_text(vcf)
    for side in (1,2):
        with (root/f'R{side}.fq').open('w') as f:
            for i,s in enumerate(seq):
                for n in range(100):
                    p = 100 + (n*13) % 2400 + (100 if side==2 else 0)
                    read = s[p:p+250]
                    f.write(f'@c{i}m{n}/{side}\n{read}\n+\n'+len(read)*'I'+'\n')
    idx = root/'panel.gindex'
    run('--vcf',root/'panel.vcf','--reference',root/'reference.fa','--write-index',idx)
    assert (idx/'manifest.dgi').is_file()
    assert (idx/'reference.fa').read_text().count('>') == 2
    common = ['--sr','-r1',root/'R1.fq','-r2',root/'R2.fq','-R','2','-t','4']
    for mode in ('index','vcf'):
        source = ['--gindex',idx] if mode=='index' else ['--vcf',root/'panel.vcf','--reference',root/'reference.fa']
        run(*source,*common,'--out',root/mode)
        with (root/mode/'read-bins.tsv').open() as f: bins=list(csv.DictReader(f,delimiter='\t'))
        assert [int(b['pairs']) for b in bins] == [100,100], bins
        assert (root/mode/'all-components.fa').read_text().count('>')==4
        assert len(list((root/mode/'components').glob('*/run.tsv'))) == 2
    assert (root/'index/all-components.fa').read_bytes()==(root/'vcf/all-components.fa').read_bytes()
    for c in ('0','1'):
        assert (root/'index/components'/c/'native-observations.tsv.gz').is_file()
    for c in ('0','1'):
        for f in ('all-components.fa','inputs.tsv','contacts.tsv','presence.tsv','pair.trace.tsv'):
            assert (root/'index/components'/c/f).read_bytes()==(root/'vcf/components'/c/f).read_bytes(),f
    # <50kb assemblies deliberately produce empty calls under the unchanged
    # dipcall gate. Both native binning and the entire caller run on empty PATH.
    before_path=os.environ['PATH'];os.environ['PATH']=''
    try:
        output=root/'combined.vcf.gz'
        run('--gindex',idx,*common,'--out',root/'called','--sample','HG002','--output-vcf',output)
        assert Path(str(output)+'.csi').exists()
        with gzip.open(output,'rt') as f: lines=f.read().splitlines()
        assert next(s for s in lines if s.startswith('#CHROM')).split('\t')[-1]=='HG002'
        assert not [s for s in lines if not s.startswith('#')]
        assert any('ID=chr1,' in s for s in lines) and any('ID=chr2,' in s for s in lines)
        assert len(list((root/'called').rglob('backend.tsv')))==2
        print('PASS empty-PATH raw native binning/caller: empty calls, sample, multi-contig header and CSI')
    finally:
        os.environ['PATH']=before_path
    # Fresh GFA index has the same oriented walks after an exact round trip.
    graph = root/'toy.gfa'
    graph.write_text('H\tVN:Z:1.0\nS\t0\t'+seq[0]+'\n'
        'P\tREF#0#chr1\t0+\t*\nP\tsample#1#chr1\t0+\t*\nP\tsample#2#chr1\t0+\t*\n')
    single=root/'single.gindex'
    run('--gfa',graph,'--REF','REF#0','--write-index',single)
    run('--gfa',graph,*common,'--out',root/'missing-ref',ok=False)
    run('--gfa',graph,'--REF','REF#0',*common,'--out',root/'gfa')
    run('--gindex',single,*common,'--out',root/'roundtrip')
    assert (root/'gfa/all-components.fa').read_bytes()==(root/'roundtrip/all-components.fa').read_bytes()
    # Existing paths cannot be overwritten, and unsupported flags fail.
    run('--gindex',single,*common,'--out',root/'roundtrip',ok=False)
    for flag in ('--ont','--hifi','--objective','--contacts','--evidence-mode','--libraries','--hic-contacts',
                 '--hic-reference','--phase-factor','--hic-split','--repair-charge','--panel-relations'):
        run('--gindex',single,*common,'--out',root/'invalid',flag,'x',ok=False)
    # Corrupt graph bytes must be detected before inference.
    shard=single/'component-0.gfa.gz';shard.write_bytes(shard.read_bytes()+b'bad')
    run('--gindex',single,*common,'--out',root/'corrupt',ok=False)
print('PASS VCF/GFA/index raw runs, full-reference paired bins, combined FASTA, exact round-trip, corruption and unsupported-option guards')
