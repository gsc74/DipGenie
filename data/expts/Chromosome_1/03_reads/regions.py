#!/usr/bin/env python3
"""Hi-C regional extraction and pair-nested sampling. Never joins mates as one span."""
import argparse
import collections
import csv
import gzip
import hashlib
import itertools
import json
import pathlib
import subprocess
import sys
import pysam

REGIONS = {'MHC': ('chr6', 28381554, 33301857), 'chr1': ('chr1', 0, 248387328)}

def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)

def query_intervals(record, region):
    """Aligned M/=X query positions inside region, in ORIGINAL read orientation."""
    chrom, start, end = region
    if record.is_unmapped or record.reference_name != chrom:
        return []
    q, p, intervals = 0, record.reference_start, []
    length = record.infer_query_length()
    for op, n in record.cigartuples or []:
        if op in (0, 7, 8):
            lo, hi = max(p, start), min(p+n, end)
            if lo < hi:
                a, b = q+lo-p, q+hi-p
                intervals.append((length-b, length-a) if record.is_reverse else (a, b))
        if op in (0, 1, 4, 7, 8): q += n
        if op in (0, 2, 3, 7, 8): p += n
    return intervals

def union_length(intervals):
    total, end = 0, -1
    for a,b in sorted(intervals):
        total += max(0, b-max(a,end)); end=max(end,b)
    return total

def primary_pair(records):
    primary=[r for r in records if not r.is_secondary and not r.is_supplementary]
    if len(primary)<2: return None  # samtools view -L excluded an outside mate
    if len(primary)!=2 or not all(r.is_paired for r in primary):
        raise ValueError('Expected two primary paired records: '+records[0].query_name)
    primary.sort(key=lambda r: not r.is_read1)
    if not primary[0].is_read1 or not primary[1].is_read2:
        raise ValueError('Invalid mate flags')
    if any(op==5 for r in records for op,n in r.cigartuples or []):
        raise ValueError('Hard clipping is not permitted: align with -Y')
    return primary

def fq(record, mate):
    s=record.get_forward_sequence(); q=record.get_forward_qualities()
    if s is None or q is None or len(s)!=len(q): raise ValueError('Missing full read')
    return f'@{record.query_name}/{mate}\n{s}\n+\n'+''.join(chr(v+33) for v in q)+'\n'

def extract(args):
    out=pathlib.Path(args.output); out.mkdir(parents=True, exist_ok=False)
    sam=pysam.AlignmentFile(args.bam, 'rb' if args.bam.endswith('.bam') else 'r')
    stats={}; streams={}; seen={}
    for name,region in REGIONS.items():
        if sam.get_reference_length(region[0]) < region[2]: raise ValueError('Reference mismatch')
        d=out/name; d.mkdir(); seen[name]=set()
        streams[name]=(open(d/'R1.fastq','w'),open(d/'R2.fastq','w'),
                      open(d/'pairs.tsv','w'),pysam.AlignmentFile(str(d/'selected.bam'),'wb',header=sam.header),
                      open(d/'segments.tsv','w'))
        streams[name][2].write('readID\tinside_bases\traw_bases\trank\n')
        streams[name][4].write('readID\tmate\tchrom\tstart0\tcigar\tstrand\tmapq\tregional_query_intervals\n')
        stats[name]=collections.Counter()
    all_pairs=0
    for read_id, group in itertools.groupby(sam,key=lambda r:r.query_name):
        records=list(group); prim=primary_pair(records); all_pairs+=1
        if prim is None: continue
        if any(c.isspace() for c in read_id): raise ValueError('Whitespace in BAM QNAME')
        for name,region in REGIONS.items():
            # Explicit user requirement: BOTH primary alignment spans must
            # be wholly inside the interval. A one-base overlap is not enough.
            if not all(not r.is_unmapped and r.reference_name==region[0] and
                       r.reference_start>=region[1] and r.reference_end<=region[2]
                       and query_intervals(r,region) for r in prim): continue
            if read_id in seen[name]: raise ValueError('Repeated/non-grouped QNAME: '+read_id)
            seen[name].add(read_id)
            r1,r2,meta,bam,segments=streams[name]
            inside=0
            for mate,record,stream in zip((1,2),prim,(r1,r2)):
                stream.write(fq(record,mate))
                valid=[r for r in records if not r.is_secondary and int(r.is_read2)+1==mate]
                inside+=union_length([i for r in valid for i in query_intervals(r,region)])
            raw=sum(len(r.query_sequence) for r in prim)
            rank=hashlib.blake2b(f'{args.seed}:{read_id}'.encode(),digest_size=16).hexdigest()
            meta.write(f'{read_id}\t{inside}\t{raw}\t{rank}\n')
            for r in records:
                bam.write(r)
                if not r.is_secondary and not r.is_unmapped:
                    segments.write(f'{read_id}\t{2 if r.is_read2 else 1}\t{r.reference_name}\t{r.reference_start}\t{r.cigarstring}\t{"-" if r.is_reverse else "+"}\t{r.mapping_quality}\t{json.dumps(query_intervals(r,region))}\n')
            stats[name].update(pairs=1,inside_bases=inside,raw_bases=raw)
        if all_pairs%1000000==0: print(json.dumps({'all_pairs':all_pairs,'selected':stats}),file=sys.stderr,flush=True)
    for name,handles in streams.items():
        for h in handles: h.close()
        region=REGIONS[name]; stats[name].update(regional_bam_qname_groups=all_pairs,reference_bases=region[2]-region[1])
        (out/name/'extraction.json').write_text(json.dumps(dict(stats[name]),indent=2)+'\n')
        (out/name/'read_ids.txt').write_text(''.join(n+'\n' for n in sorted(seen[name])))
    sam.close()

def fastq(path):
    with open_text(path) as f:
        while header:=f.readline():
            seq=f.readline().rstrip(); plus=f.readline().rstrip(); qual=f.readline().rstrip()
            if not header.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):
                raise ValueError('Bad FASTQ: '+str(path))
            yield header.split()[0][1:],seq,qual

def pairs_rows(path):
    columns=None
    with open_text(path) as f:
        for line in f:
            if line.startswith('#columns:'): columns=line.strip().split()[1:]
            if line.startswith('#'): continue
            if columns is None: raise ValueError('Missing pairs columns')
            yield dict(zip(columns,line.rstrip('\n').split('\t')))

def nesting_members(ranked, length, depths):
    # Exact hash-prefix selection: target depth reached to within one pair's
    # LOCAL aligned bases. Both mates always retained. Never counts mate gap.
    total=sum(x[2] for x in ranked)
    available=[d for d in depths if d*length < total]
    thresholds=[d*length for d in available]
    buckets={}; cumulative=0; bucket=0
    for rank,name,inside,raw in sorted(ranked):
        while bucket<len(thresholds) and cumulative>=thresholds[bucket]: bucket+=1
        buckets[name]=bucket
        cumulative+=inside
    return available,buckets,total

def pair_rank(seed, name):
    """One reproducible pseudo-random priority per physical molecule."""
    if seed < 0:
        raise ValueError('Sampling seed must be nonnegative')
    return hashlib.blake2b(f'{seed}:{name}'.encode(),digest_size=16).hexdigest()

class Compressed:
    def __init__(self,path):
        self.file=open(path,'wb'); self.proc=subprocess.Popen(['pigz','-1','-p','2','-c'],stdin=subprocess.PIPE,stdout=self.file,text=True)
    def write(self,s): self.proc.stdin.write(s)
    def close(self):
        self.proc.stdin.close(); code=self.proc.wait(); self.file.close()
        if code: raise RuntimeError('pigz failed')

def sample(args):
    source=pathlib.Path(args.input); out=pathlib.Path(args.output); out.mkdir(parents=True,exist_ok=False)
    duplicates={r['readID'] for r in pairs_rows(source/'duplicates.pairs.gz')}
    ranked=[]; seen=set()
    with open(source/'pairs.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['readID'] in seen: raise ValueError('Duplicate input ID')
            seen.add(r['readID'])
            if r['readID'] not in duplicates:
                # Extraction's stored rank belongs to its original seed. Rehash
                # here so --seed really changes the sample, not just its receipt.
                ranked.append((pair_rank(args.seed,r['readID']),r['readID'],int(r['inside_bases']),int(r['raw_bases'])))
    length=REGIONS[args.region][2]-REGIONS[args.region][1]
    available,buckets,total=nesting_members(ranked,length,[2,4,8,16,32])
    if not ranked: raise ValueError('No regional pairs remain')
    labels=[f'{d}x' for d in available]+['Total']
    meta={n:(inside,raw) for rank,n,inside,raw in ranked}
    handles=[]; stats=[collections.Counter() for _ in labels]
    for label in labels:
        d=out/label; d.mkdir()
        handles.append([Compressed(d/'HG002.R1.fastq.gz'),Compressed(d/'HG002.R2.fastq.gz'),Compressed(d/'HG002.interleaved.fastq.gz')])
    consumed=0
    with open(out/'membership.tsv','w') as members:
        members.write('readID\tfirst_depth_index\tinside_bases\traw_bases\n')
        for a,b in itertools.zip_longest(fastq(source/'R1.fastq'),fastq(source/'R2.fastq')):
            if a is None or b is None or not a[0].endswith('/1') or not b[0].endswith('/2') or a[0][:-2]!=b[0][:-2]:
                raise ValueError('Broken FASTQ pairing')
            name=a[0][:-2]
            if name in duplicates: continue
            bucket=buckets.pop(name); inside,raw=meta.pop(name); consumed+=1
            records=[f'@{x[0]}\n{x[1]}\n+\n{x[2]}\n' for x in (a,b)]
            if raw!=len(a[1])+len(b[1]): raise ValueError('Sequence/metadata mismatch')
            members.write(f'{name}\t{bucket}\t{inside}\t{raw}\n')
            for i in range(bucket,len(labels)):
                for j in (0,1): handles[i][j].write(records[j])
                handles[i][2].write(''.join(records))
                stats[i].update(pairs=1,inside_bases=inside,raw_bases=raw)
    if buckets or meta: raise ValueError('Metadata IDs missing in FASTQ')
    for h in handles:
        for f in h: f.close()
    # Contact subsets use precisely the SAME pair membership, never independent
    # sampling and never use Total contacts for a lower-depth run.
    membership={}
    with open(out/'membership.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'): membership[r['readID']]=int(r['first_depth_index'])
    contacts=[Compressed(out/label/'contacts.pairs.gz') for label in labels]
    region=REGIONS[args.region]
    with open_text(source/'nodups.pairs.gz') as f:
        columns=None
        for line in f:
            if line.startswith('#'):
                if line.startswith('#columns:'): columns=line.strip().split()[1:]
                for h in contacts: h.write(line)
                continue
            if columns is None: raise ValueError('Missing contact header')
            r=dict(zip(columns,line.rstrip('\n').split('\t'))); name=r['readID']
            if name not in membership: continue
            if r['pair_type']!='UU': continue
            inside=[r[f'chrom{k}']==region[0] and region[1]<=int(r[f'pos{k}'])-1<region[2] for k in (1,2)]
            if not all(inside): continue
            for i in range(membership[name],len(labels)):
                contacts[i].write(line); stats[i]['internal_contacts']+=1
    for h in contacts: h.close()
    rows=[]
    for label,s in zip(labels,stats):
        rows.append(dict(region=args.region,depth_label=label,**s,depth=s['inside_bases']/length,
                         raw_sequence_depth=s['raw_bases']/length,reference_bases=length))
    with open(out/'coverage.csv','w') as f:
        fields=['region','depth_label','pairs','inside_bases','raw_bases','internal_contacts','depth','raw_sequence_depth','reference_bases']
        w=csv.DictWriter(f,fieldnames=fields,restval=0); w.writeheader(); w.writerows(rows)
    (out/'PASS.json').write_text(json.dumps({'nested_by_pair_hash_prefix':True,'paired_ids_verified':consumed,
        'excluded_duplicate_pairs':len(duplicates),'seed':args.seed,'coverage_basis':'union of aligned M/=X query bases inside region, divided by reference length',
        'unavailable_depths':[d for d in [2,4,8,16,32] if d not in available],
        'pair_selection':'both primary alignment spans fully contained in region; no outside-mate rescue',
        'contact_policy':'UU MAPQ>=30, deduplicated; both endpoints inside region',
        'rows':rows},indent=2)+'\n')

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='mode',required=True)
    e=sub.add_parser('extract'); e.add_argument('--bam',required=True); e.add_argument('--output',required=True); e.add_argument('--seed',type=int,default=42); e.set_defaults(func=extract)
    s=sub.add_parser('sample'); s.add_argument('--input',required=True); s.add_argument('--output',required=True); s.add_argument('--region',choices=REGIONS,required=True); s.add_argument('--seed',type=int,default=42); s.set_defaults(func=sample)
    args=p.parse_args(); args.func(args)
if __name__=='__main__': main()
