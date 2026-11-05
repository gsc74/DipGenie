#!/usr/bin/env python3
"""Raw paired FASTQ + GFA integration test, with no inferred initializer and an empty PATH."""
import csv
import pathlib
import random
import subprocess
import tempfile
import argparse
import contextlib
import os
import gzip
import math

HERE=pathlib.Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--cyclic',action='store_true')
parser.add_argument('--dipcall',action='store_true')
parser.add_argument('--recombinant',action='store_true')
parser.add_argument('--unobserved-alleles',action='store_true',help='Extra panel-only alleles must be scored for absence, not added to phase markers')
parser.add_argument('--cross-locus',action='store_true',help='One copy in every intact donor, but two copies on one true recombinant haplotype')
parser.add_argument('--executable',type=pathlib.Path,default=HERE/'DipGenie')
parser.add_argument('--libraries',choices=('hic','sr','mixed'),default='hic')
parser.add_argument('--compare-executable',type=pathlib.Path,
                    help='Run the same generated inputs through another binary and require identical outputs')
options=parser.parse_args()
options.executable=options.executable.resolve()
if options.compare_executable:options.compare_executable=options.compare_executable.resolve()
rng=random.Random(42)
def dna(n): return ''.join(rng.choice('ACGT') for _ in range(n))
common=[dna(200) for _ in range(7)]
if options.dipcall:
    common[0]=dna(60000)+common[0]
    common[-1]+=dna(60000)
pair_count=18000 if options.dipcall else 1800
alleles=[]
for b in range(6):
    a=dna(200); alt=list(a)
    for p in range(8,200,17):alt[p]=rng.choice([c for c in 'ACGT' if c!=alt[p]])
    alleles.append([a,''.join(alt)])
donors=[[0]*6,[1]*6,[0,1,0,1,0,1],[1,0,1,0,1,0]]
if options.recombinant:donors=[[0,1,0,1,0,1],[1,0,1,0,1,0],[0,0,1,1,0,0],[1,1,0,0,1,1]]
if options.cross_locus:
    alleles[4][0]=alleles[0][0]
    for donor in donors:donor[4]=1-donor[0]
copies=lambda b:2 if options.cyclic and b==2 else 1
expected=[''.join(common[b]+alleles[b][h]*copies(b) for b in range(6))+common[6] for h in range(2)]
if options.unobserved_alleles:
    for words in alleles:words.append(dna(200))
    donors.append([2]*6)
# Preserve real-caller work and logs, including failed calls, for diagnosis.
workspace=contextlib.nullcontext(tempfile.mkdtemp(prefix='hic-raw-dipcall-')) if options.dipcall else tempfile.TemporaryDirectory(prefix='hic-raw-integration-')
with workspace as temp:
    d=pathlib.Path(temp)
    if options.dipcall:print('dipcall toy workspace:',d,flush=True)
    with open(d/'graph.gfa','w') as f:
        for b,s in enumerate(common):f.write(f'S\tc{b}\t{s}\n')
        for b,words in enumerate(alleles):
            for a,s in enumerate(words):
                f.write(f'S\tv{b}_{a}\t{s}\nL\tc{b}\t+\tv{b}_{a}\t+\t0M\nL\tv{b}_{a}\t+\tc{b+1}\t+\t0M\n')
                if copies(b)>1:f.write(f'L\tv{b}_{a}\t+\tv{b}_{a}\t+\t0M\n')
        for h,choices in enumerate(donors):
            walk=''.join(f'>c{b}'+f'>v{b}_{a}'*copies(b) for b,a in enumerate(choices))+'>c6'
            f.write(f'W\tpanel{h}\t0\ttoy\t0\t{len(expected[0])}\t{walk}\n')
    with open(d/'R1.fq','w') as a,open(d/'R2.fq','w') as b:
        for pair in range(pair_count):
            hap=rng.randrange(2);sites=rng.sample(range(len(expected[hap])-99),2)
            for side,(site,f) in enumerate(zip(sites,(a,b)),1):
                read=expected[hap][site:site+100]
                f.write(f'@m{pair}/{side}\n{read}\n+\n'+len(read)*'I'+'\n')
    if options.libraries in ('sr','mixed'):
        with open(d/'SR1.fq','w') as a,open(d/'SR2.fq','w') as b:
            for pair in range(pair_count):
                hap=rng.randrange(2);site=rng.randrange(len(expected[hap])-499)
                for side,(start,f) in enumerate(zip((site,site+400),(a,b)),1):
                    read=expected[hap][start:start+100]
                    f.write(f'@s{pair}/{side}\n{read}\n+\n'+100*'I'+'\n')
    cmd=[str(options.executable),'--gfa',str(d/'graph.gfa'),'--out',str(d/'out'),'-t','8','-R','18']
    if options.libraries in ('sr','mixed'):
        cmd+=['--sr','-r1',str(d/'SR1.fq'),'-r2',str(d/'SR2.fq')]
    if options.libraries in ('hic','mixed'):
        cmd+=['--hic','-h1',str(d/'R1.fq'),'-h2',str(d/'R2.fq')]
    cmd+=['--REF','panel0#0#toy']
    if options.dipcall:
        cmd+=['--output-vcf',str(d/'calls.vcf.gz'),'--sample','HG002']
    env=dict(os.environ,PATH='')
    subprocess.run(cmd,check=True,env=env)
    if options.compare_executable:
        reference=d/'out-reference'
        reference_cmd=list(cmd)
        reference_cmd[0]=str(options.compare_executable)
        reference_cmd[reference_cmd.index('--out')+1]=str(reference)
        if options.dipcall:reference_cmd[reference_cmd.index('--output-vcf')+1]=str(d/'reference.vcf.gz')
        subprocess.run(reference_cmd,check=True,env=env,capture_output=True)
        for name in ('all-components.fa','inputs.tsv','contacts.tsv','contacts.stage1.tsv','presence.tsv','pair.trace.tsv',
                     'library-counts.tsv','libraries.tsv','calibration-histogram.tsv','dosage-fit.tsv',
                     'native-factors.stage1.tsv','native-factors.uncalibrated.tsv','linkage-calibration.tsv',
                     'native-calibration-presence.tsv','orientation-search.tsv','global-phase-moves.tsv',
                     'unphased-segments.bed','copy-number-windows.tsv'):
            assert (reference/name).read_bytes()==(d/'out'/name).read_bytes(),name
        def stable_metrics(path):
            with path.open() as stream:
                return {r['metric']:r['value'] for r in csv.DictReader(stream,delimiter='\t') if r['metric']!='seconds'}
        assert stable_metrics(reference/'run.tsv')==stable_metrics(d/'out/run.tsv')
        if options.dipcall:
            def records(path):
                with gzip.open(path,'rt') as f:return [s for s in f.read().splitlines() if not s.startswith('##')]
            assert records(d/'reference.vcf.gz')==records(d/'calls.vcf.gz')
        print('PASS reference-binary parity: same inputs, calibration, factors, objective, trace, FASTAs, certificate'+(' and VCF' if options.dipcall else ''))
    if options.dipcall:
        assert (d/'calls.vcf.gz.csi').is_file()
        import gzip
        with gzip.open(d/'calls.vcf.gz','rt') as f:lines=f.read().splitlines()
        assert next(s for s in lines if s.startswith('#CHROM')).split('\t')[-1]=='HG002'
        records=[s.split('\t') for s in lines if not s.startswith('#')]
        assert records and all('|' in row[9].split(':')[0] for row in records)
        print('PASS empty-PATH raw reads to normalized indexed phased VCF with HG002 sample')
    with open(d/'out/run.tsv') as f: metrics={r['metric']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
    with (d/'out/stages.tsv').open() as f:ledger=list(csv.DictReader(f,delimiter='\t'))
    elapsed=0
    for row in ledger:
        begin,end,duration=(float(row[k]) for k in ('start_seconds','end_seconds','elapsed_seconds'))
        assert begin==elapsed and end>=begin and math.isclose(end-begin,duration,abs_tol=1e-12)
        elapsed=end
    assert elapsed>=float(metrics['seconds'])
    assert (d/'out/native-observations.tsv.gz').is_file()
    # End-to-end determinism includes placement/admission and calibrated
    # factors, not just equality of the final FASTAs.
    stable=('all-components.fa','inputs.tsv','contacts.tsv','presence.tsv','pair.trace.tsv',
            'library-counts.tsv','native-factors.stage1.tsv','native-factors.uncalibrated.tsv',
            'native-calibration-presence.tsv','linkage-calibration.tsv')
    inference_cmd=list(cmd)
    if '--output-vcf' in inference_cmd:
        at=inference_cmd.index('--output-vcf');del inference_cmd[at:at+2]
    for threads in (1,2,4):
        replica=d/f'out-native-{threads}';again=list(inference_cmd)
        again[again.index('-t')+1]=str(threads);again[again.index('--out')+1]=str(replica)
        subprocess.run(again,check=True,env=dict(os.environ,PATH=''),capture_output=True)
        for name in stable:assert (d/'out'/name).read_bytes()==(replica/name).read_bytes(),(threads,name)
        with gzip.open(d/'out/native-observations.tsv.gz','rb') as a,gzip.open(replica/'native-observations.tsv.gz','rb') as b:
            assert a.read()==b.read(),threads
    print('PASS end-to-end thread parity 1/2/4/8')
    if options.unobserved_alleles:
        with (d/'out/inputs.tsv').open() as f:targets=list(csv.DictReader(f,delimiter='\t'))
        zeros={int(r['target']) for r in targets if int(r['dosage_read_end_count'])==0}
        assert zeros and len(zeros)==int(metrics['zero_dosage_targets'])
        assert any(float(targets[x]['absent_probability'])>.99 for x in zeros)
        with (d/'out/contacts.tsv').open() as f:
            for r in csv.DictReader(f,delimiter='\t'):
                assert int(r['target1']) not in zeros and int(r['target2']) not in zeros
        print('PASS raw panel-only alleles retained as zero-count dosage evidence; excluded from phase-marker eligibility')
    if options.cross_locus:
        assert int(metrics['presence_only_targets'])>0
        with (d/'out/inputs.tsv').open() as source:
            cross=[r for r in csv.DictReader(source,delimiter='\t') if r['first_segment']!=r['last_segment']]
        assert cross and any(r['panel_copy_histogram']=='1:4' for r in cross)
        neutral=str(round(1024/3))
        assert all(r['reward_units']==r['other_reward_units']==neutral and float(r['hom_probability'])==.5 for r in cross)
    seqs=[]
    for line in (d/'out/all-components.fa').read_text().splitlines():
        if line.startswith('>'):seqs.append('')
        else:seqs[-1]+=line
    assert len(seqs)==2 and all(len(s)==len(expected[0]) for s in seqs)
    assert int(metrics['raw_pairs'])==pair_count*(2 if options.libraries=='mixed' else 1) and int(metrics['informative_contact_molecules'])>0
    assert int(metrics['raw_molecule_weight_units'])==1024*int(metrics['informative_contact_molecules'])
    confidence=int(metrics['linkage_confidence_scale'])
    assert 0<int(metrics['read_contact_weight_units'])<=confidence*int(metrics['raw_molecule_weight_units'])
    assert int(metrics['reward_integer'])<=int(metrics['segment_upper_integer'])<=int(metrics['full_graph_upper_integer'])
    assert sorted(seqs)==sorted(expected), 'Toy ground-truth diploid pair not recovered'
    # Replay from emitted factor tables, independently of the C++ scorer.
    import importlib.util
    audit_path=HERE/'audit_model.py'
    spec=importlib.util.spec_from_file_location('toy_factor_audit',audit_path)
    independent=importlib.util.module_from_spec(spec);spec.loader.exec_module(independent)
    assert independent.audit(d/'out')['status']=='PASS'
    print('PASS independent factor replay')
    print('PASS raw '+('cyclic' if options.cyclic else 'DAG')+f' GFA + {options.libraries} raw pairs -> exact expected diploid FASTAs; '+str(metrics))
