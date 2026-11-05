import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pysam
import prepare_chr1_sr as p
import audit_inputs


class Chr1SamplingTests(unittest.TestCase):
    def test_strict_pair_selection_and_nesting(self):
        with tempfile.TemporaryDirectory(prefix='chr1-pair-toy-') as tmp:
            root=Path(tmp);source=root/'source';d=source/'extracted';d.mkdir(parents=True)
            (d/'selection.json').write_text(json.dumps(dict(source_accession='SRR11321732',region='chr1:1-248387328')))
            (source/'extraction.command').write_text('bwa-mem2 mem chm13v2.0.fa /Input/SR/Total/HG002.SRR11321732.Total.R1.fastq.gz /Input/SR/Total/HG002.SRR11321732.Total.R2.fastq.gz')
            fq=[[],[]]
            with pysam.AlignmentFile(str(d/'selected.primary.bam'),'wb',header={'HD':{'VN':'1.6','SO':'queryname'},'SQ':[{'SN':'chr1','LN':100},{'SN':'chr2','LN':100}]}) as bam:
                for i in range(82):
                    name=f'pair{i:03}'
                    for mate in (0,1):
                        r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*30
                        r.query_qualities=pysam.qualitystring_to_array('I'*30);r.flag=1+(64 if mate==0 else 128)
                        r.reference_id=1 if i==80 and mate==1 else 0;r.reference_start=10+mate*30;r.cigarstring='30M';r.mapping_quality=60
                        if i==81 and mate==1:r.flag|=4;r.reference_id=-1;r.reference_start=-1;r.cigarstring=None
                        bam.write(r);fq[mate].append(f'@{name}/{mate+1}\n'+r.query_sequence+'\n+\n'+'I'*30+'\n')
            for mate in (0,1):(d/f'HG002.R{mate+1}.fastq').write_text(''.join(fq[mate]))
            with patch.dict(p.REGIONS,chr1=('chr1',0,100)):
                result=p.prepare(source,root/'out',42)
                self.assertEqual(result['paired_ids_verified'],80)
                self.assertEqual(result['excluded_noninternal_pairs'],2)
                self.assertEqual(audit_inputs.audit(root/'out','chr1')['status'],'PASS')
                with self.assertRaises(AssertionError):p.prepare(source,root/'out',42)


if __name__=='__main__':unittest.main()
