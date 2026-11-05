import argparse
import csv
import os
import gzip
import pathlib
import subprocess
import tempfile
import unittest
import pysam
import regions as m

HERE=pathlib.Path(__file__).resolve().parent
PT=os.environ.get('PAIRTOOLS', 'pairtools')

def record(header,name,mate,pos,chrom='chr6',reverse=False,cigar='100M',supp=False):
    r=pysam.AlignedSegment(header)
    r.query_name=name; r.flag=1+(64 if mate==1 else 128)+(16 if reverse else 0)+(2048 if supp else 0)
    r.reference_id=header.get_tid(chrom); r.reference_start=pos; r.mapping_quality=60
    r.cigarstring=cigar; r.query_sequence='ACGT'*25; r.query_qualities=pysam.qualitystring_to_array('I'*100)
    return r

class Tests(unittest.TestCase):
    def setUp(self):
        self.header=pysam.AlignmentHeader.from_references(['chr1','chr6'],[248387328,172126628])
    def test_cigar_boundaries_and_reverse(self):
        region=('chr6',100,200)
        self.assertEqual(m.query_intervals(record(self.header,'a',1,0),region),[])
        self.assertEqual(m.query_intervals(record(self.header,'a',1,1),region),[(99,100)])
        self.assertEqual(m.query_intervals(record(self.header,'a',1,199,reverse=True),region),[(99,100)])
        self.assertEqual(m.query_intervals(record(self.header,'a',1,90,cigar='10S40M10D50M'),region),[(20,50),(50,100)])
        self.assertEqual(m.union_length([(0,50),(30,80),(79,100)]),100)
    def test_prefix_nesting_variable_lengths(self):
        rows=[(f'{i:04}',str(i),i%13+4,100) for i in range(1000)]
        depths,buckets,total=m.nesting_members(rows,100,[2,4,8,16,32])
        prev=set()
        for k,d in enumerate(depths):
            ids={n for n in buckets if buckets[n]<=k}
            self.assertTrue(prev<=ids)
            bases=sum(r[2] for r in rows if r[1] in ids)
            self.assertGreaterEqual(bases,d*100)
            self.assertLess(bases-d*100,17)
            prev=ids
        self.assertEqual(m.nesting_members(list(reversed(rows)),100,[2,4,8,16,32]),
                         m.nesting_members(rows,100,[2,4,8,16,32]))
    def test_seed_rank_contract(self):
        import hashlib
        self.assertEqual(m.pair_rank(42,'read1'),
                         hashlib.blake2b(b'42:read1',digest_size=16).hexdigest())
        self.assertNotEqual(m.pair_rank(42,'read1'),m.pair_rank(43,'read1'))
        with self.assertRaises(ValueError):m.pair_rank(-1,'read1')
    def test_pairtools_and_extraction(self):
        with tempfile.TemporaryDirectory(prefix='hic-extract-test-') as td:
            td=pathlib.Path(td); bam=td/'toy.bam'; start=m.REGIONS['MHC'][1]
            with pysam.AlignmentFile(bam,'wb',header=self.header) as out:
                # Internal pair, exact PCR duplicate, long-range boundary pair,
                # one-base overlap, supplementary non-overlapping region hit,
                # and a chr1 pair. No FR or insert-size assumption is allowed.
                rows=[('internal',start+100,start+10000),('duplicate',start+100,start+10000),
                      ('boundary',start+500,m.REGIONS['MHC'][2]+100),('onebase',start-99,start-1000)]
                for name,a,b in rows:
                    out.write(record(self.header,name,1,a)); out.write(record(self.header,name,2,b,reverse=True))
                out.write(record(self.header,'split',1,start+1000,cigar='50M50S'))
                out.write(record(self.header,'split',1,start+1010,cigar='30S50M20S',supp=True))
                out.write(record(self.header,'split',2,start+20000))
                out.write(record(self.header,'chr1',1,100,chrom='chr1'))
                out.write(record(self.header,'chr1',2,100000,chrom='chr1'))
            (td/'regions.bed').write_text(f'chr6\t{start}\t{m.REGIONS["MHC"][2]}\nchr1\t0\t248387328\n')
            with open(td/'regional.sam','w') as stream:
                subprocess.run(['samtools','view','-h','-L',str(td/'regions.bed'),str(bam)],stdout=stream,check=True)
            m.extract(argparse.Namespace(bam=str(td/'regional.sam'),output=str(td/'extract'),seed=42))
            region=td/'extract/MHC'
            with open(region/'pairs.tsv') as stream: meta=list(csv.DictReader(stream,delimiter='\t'))
            self.assertEqual(len(meta),3)
            self.assertEqual(next(int(r['inside_bases']) for r in meta if r['readID']=='split'),180)
            self.assertFalse({'onebase','boundary'} & {r['readID'] for r in meta})
            # Restore all records from the full BAM by the retained pair IDs.
            subprocess.run(['samtools','view','-b','-N',str(region/'read_ids.txt'),'-o',str(region/'all-alignments.bam'),str(bam)],check=True)
            (td/'chrom.sizes').write_text('chr1\t248387328\nchr6\t172126628\n')
            subprocess.run([PT,'parse','--min-mapq','30','--walks-policy','5unique','--drop-sam',
                            '--add-columns','mapq,read_side,algn_idx,cigar','-c',str(td/'chrom.sizes'),
                            '-o',str(region/'raw.pairs'),str(region/'all-alignments.bam')],check=True,capture_output=True)
            subprocess.run([PT,'sort','--nproc','1','-o',str(region/'sorted.pairs'),str(region/'raw.pairs')],check=True,capture_output=True)
            subprocess.run([PT,'dedup','--max-mismatch','3','-o',str(region/'nodups.pairs.gz'),
                            '--output-dups',str(region/'duplicates.pairs.gz'),str(region/'sorted.pairs')],check=True,capture_output=True)
            dup=list(m.pairs_rows(region/'duplicates.pairs.gz'))
            self.assertEqual(len(dup),1)
            m.sample(argparse.Namespace(input=str(region),output=str(td/'samples'),region='MHC',seed=42))
            import json
            result=json.loads((td/'samples/PASS.json').read_text())
            self.assertEqual(result['paired_ids_verified'],2)
            self.assertEqual(result['unavailable_depths'],[2,4,8,16,32])
            reads=list(m.fastq(td/'samples/Total/HG002.interleaved.fastq.gz'))
            self.assertEqual(len(reads),4)
            # Original read orientations and pairing survive SAM reconstruction.
            for a,b in zip(reads[::2],reads[1::2]): self.assertEqual(a[0][:-2],b[0][:-2])

if __name__=='__main__': unittest.main()
