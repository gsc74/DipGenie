#!/usr/bin/env python3
"""Independent arbitrary-precision replay of soft dosage and phase reward."""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path

def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter='\t')

def audit(directory):
    d=Path(directory)
    m={x['metric']:x['value'] for x in rows(d/'run.tsv')}
    assert m['objective_model']=='zero_aware_dosage_all_library_phase'
    assert m['calibration_model']=='single_segment_panel_copy_NB_cross_segment_presence'
    assert int(m['objective_integer_scale'])==2
    scale=int(m['confidence_scale']);wg=int(m['genotype_weight_integer']);wc=int(m['contact_weight_integer'])
    targets=list(rows(d/'inputs.tsv'));presence=list(rows(d/'presence.tsv'))
    assert len(targets)==len(presence)
    libraries=list(rows(d/'libraries.tsv'))
    active=[x['library'] for x in libraries if int(x['dosage_evidence'])]
    has_sr=any(x['kind']=='SR' for x in libraries)
    assert active==[x['library'] for x in libraries if not has_sr or x['kind']=='SR']
    assert m['dosage_library_policy']==('SR_only' if has_sr else 'all_available')
    count_rows=list(rows(d/'library-counts.tsv'))
    assert len(count_rows)==len(targets)
    for target,counts in zip(targets,count_rows):
        assert target['target']==counts['target']
        assert int(target['read_end_count'])==sum(int(counts[x['library']]) for x in libraries)
        assert int(target['dosage_read_end_count'])==sum(int(counts[name]) for name in active)>=0
    fit={r['metric']:float(r['value']) for r in rows(d/'dosage-fit.tsv')}
    def components(n,copies):
        """Complete-count log likelihood for a total diploid copy number."""
        mean,shape=fit['mean'],fit['shape']
        def nb(n,mu):
            if math.isinf(shape):return n*math.log(mu)-mu-math.lgamma(n+1)
            return math.lgamma(n+shape)-math.lgamma(shape)-math.lgamma(n+1)+shape*math.log(shape/(shape+mu))+n*math.log(mu/(shape+mu))
        rate=fit['error_rate']
        if not copies:
            if not n:return math.log(fit['error_zero'])
            if not rate:return math.log1p(-fit['error_zero']) if n==1 else -math.inf
            lognorm=rate+math.log1p(-math.exp(-rate)) if rate>1 else math.log(math.expm1(rate))
            return math.log1p(-fit['error_zero'])+n*math.log(rate)-math.lgamma(n+1)-lognorm
        body=math.log1p(-fit['dropout'])+nb(n,mean*copies/2)
        return body if n else math.log(fit['dropout']+math.exp(body))
    def copy_level(n,classes,level):
        logs={};top=-math.inf
        for i,a in classes:
            for j,b in classes:
                v=math.log(a)+math.log(b)+components(n,i+j);logs[(i,j)]=v;top=max(top,v)
        p=[0.0,0.0,0.0];total=sum(math.exp(v-top) for v in logs.values())
        for (i,j),v in logs.items():p[int(i>=level)+int(j>=level)]+=math.exp(v-top)/total
        return p
    def posterior(n,carriers=None):
        mean,shape=fit['mean'],fit['shape']
        def nb(n,mu):
            if math.isinf(shape):return n*math.log(mu)-mu-math.lgamma(n+1)
            return math.lgamma(n+shape)-math.lgamma(shape)-math.lgamma(n+1)+shape*math.log(shape/(shape+mu))+n*math.log(mu/(shape+mu))
        rate=fit['error_rate']
        assert 0<=rate<=mean/2
        if not n:error=math.log(fit['error_zero'])
        elif not rate:error=math.log1p(-fit['error_zero']) if n==1 else -math.inf
        else:
            lognorm=rate+math.log1p(-math.exp(-rate)) if rate>1 else math.log(math.expm1(rate))
            error=math.log1p(-fit['error_zero'])+n*math.log(rate)-math.lgamma(n+1)-lognorm
        dropout=fit['dropout']
        def present(mu):
            body=math.log1p(-dropout)+nb(n,mu)
            return body if n else math.log(dropout+math.exp(body))
        logs=[error,present(mean/2),present(mean)]
        if m.get('dosage_prior','global')=='panel':
            f=(carriers+.5)/(int(m['panel_paths'])+1);prior=[(1-f)**2,2*f*(1-f),f*f]
        else:prior=[fit[f'prior{i}'] for i in range(3)]
        logs=[value+math.log(prior[i]) for i,value in enumerate(logs)]
        top=max(logs);p=[math.exp(value-top) for value in logs];total=sum(p)
        return [v/total for v in p]
    copies=0
    states=[];g=upper=0;ht=tt=at=0;covered=Counter();presence_only=0;error_aware=0
    for i,(x,y) in enumerate(zip(targets,presence)):
        assert int(x['target'])==int(y['target'])==i
        a,b=int(y['haplotype1']),int(y['haplotype2']);assert a in (0,1) and b in (0,1)
        h=int(x['hom']);assert h in (0,1)
        one,other=int(x['reward_units']),int(x['other_reward_units'])
        aw=int(x['absent_reward_units'])
        logodds=m.get('dosage_weights','probability')=='logodds'
        assert min(one,other,aw)>=0 and (one+other+aw<=scale if logodds else one+other+aw==scale)
        hw,tw=(one,other) if h else (other,one)
        if logodds:
            cap=float(m['log_odds_cap']);ph=float(x['hom_probability']);pa=float(x['absent_probability']);pt=max(0.0,1-ph-pa)
            ps=[pa,pt,ph];ls=[math.log(v) if v>0 else -math.inf for v in ps];top=max(ls)
            want=[round(max(0.0,cap-(top-l))*scale/(3*cap)) if l>-math.inf else 0 for l in ls]
            assert abs(aw-want[0])<=1 and abs(tw-want[1])<=1 and abs(hw-want[2])<=1,(want,aw,tw,hw)
        if x['first_segment']!=x['last_segment']:
            presence_only+=1
            assert (hw==tw==round(scale/3) if logodds else hw==tw==scale//2) and float(x['hom_probability'])==0.5
            assert aw==0
        if x.get('copy_level','0') not in ('0',''):
            copies+=1;level=int(x['copy_level'])
            classes=[c for c in (tuple(map(int,c.split(':'))) for c in x['copy_classes'].split(',')) if c[1]] # empty zero class omitted, as in C++
            assert x['graph_unique_locus']=='0' and x['first_segment']==x['last_segment'] and 1<=level<=max(c for c,_ in classes)
            p=copy_level(int(x['copy_occurrence_count']),classes,level)
            assert int(x['copy_occurrence_count'])>=int(x['dosage_read_end_count'])
            assert abs(p[0]-float(x['absent_probability']))<1e-6 and abs(p[2]-float(x['hom_probability']))<1e-6
        if x['graph_unique_locus']=='1' and x['first_segment']==x['last_segment']:
            error_aware+=1
            carriers=int(x['panel_carriers'])
            p=posterior(int(x['dosage_read_end_count']),carriers)
            assert abs(p[0]-float(x['absent_probability']))<1e-6
            assert abs(p[2]-float(x['hom_probability']))<1e-6
        if not logodds:
            assert abs(hw/scale-float(x['hom_probability']))<=1/scale+1e-14
            assert abs(aw/scale-float(x['absent_probability']))<=1/scale+1e-14
        ht+=hw;tt+=tw;at+=aw;g+=2*wg*(hw*a*b+tw*(a-b)**2+aw*(1-a)*(1-b))
        upper+=2*wg*max(hw,tw,aw);states.append(a+2*b)
        covered['hom' if h else 'het']+=int(a*b if h else (a-b)**2)
    sources=Counter();outcomes=Counter();c=0
    for e in rows(d/'contacts.tsv'):
        x,y=int(e['target1']),int(e['target2']);w=int(e['weight_units']);r=int(e['relation'])
        assert x!=y and 0<=x<len(states) and 0<=y<len(states) and w>0 and r in (-1,1)
        assert e['source'] in ('panel','molecule')
        a,b=states[x],states[y];sx=int(a==1)-int(a==2);sy=int(b==1)-int(b==2)
        kind=m.get('phase_factor_stage2',m.get('phase_factor','symmetric'))
        if kind=='presence':f=2 if ((a&b&3) if r>0 else ((a&1 and b&2) or (a&2 and b&1))) else 0
        elif kind=='neutral':f=0 if r*sx*sy<0 else 2
        else:f=1+r*sx*sy
        c+=wc*w*f;sources[e['source']]+=w;upper+=2*wc*w
        outcomes['nonheterozygous' if not sx*sy else 'agreeing' if r*sx*sy==1 else 'opposite']+=w
    W=sum(sources.values())
    assert presence_only==int(m['presence_only_targets'])
    if 'copy_strings' in m:assert copies==int(m['copy_strings'])+int(m['copy_level_targets'])
    assert ht==int(m['hom_reward_units']) and tt==int(m['het_reward_units'])
    assert at==int(m['absent_reward_units']) and error_aware==int(m['error_aware_targets'])==int(fit['observations'])
    zeros=sum(int(x['dosage_read_end_count'])==0 for x in targets)
    assert zeros==int(m['zero_dosage_targets'])==int(fit['zero_observations'])
    assert W==int(m['contact_weight_units'])
    calibrated=m.get('linkage_calibration','none')=='distance'
    native=m.get('evidence_mode')=='native'
    factor_den=int(m.get('linkage_factor_denominator','1'))
    assert factor_den in (1,1024)
    mass_scale=factor_den if native else (int(m['linkage_confidence_scale']) if calibrated else 1)
    assert sources['molecule']==int(m['read_contact_weight_units'])<=int(m['raw_molecule_weight_units'])*mass_scale
    if native and calibrated and m['recombination_charge']=='unordered_pair':
        # Replay the fixed pre-confidence classes independently of C++:
        # source/library/geometry/QC separation, eta, common Lmax, rounding,
        # compatible-factor pooling and the final normalization mass.
        eta=int(m['linkage_pseudocount_units']);cs=int(m['linkage_confidence_scale'])
        assert eta==512 and cs==1024
        reference={int(r['target']):int(r['presence']) for r in rows(d/'native-calibration-presence.tsv')}
        assert set(reference)==set(range(len(targets)))
        factors=list(rows(d/'native-factors.uncalibrated.tsv'))
        classes={r['class']:r for r in rows(d/'linkage-calibration.tsv')}
        expected={key:[0,0] for key in classes}
        for e in factors:
            key=e['class'];x,y=int(e['target1']),int(e['target2'])
            assert key in classes and (key=='panel')==(e['source']=='panel')
            if key!='panel':
                fields=key.split('/')
                assert fields[0] in ('SR','HIC') and fields[1].startswith('library=')
                assert fields[2] in ('unknown','trans') or fields[2].startswith('anchor_bin=')
                assert fields[3] in ('qc=1','qc=2')
            a,b=reference[x],reference[y];s=int(e['relation'])*(int(a==1)-int(a==2))*(int(b==1)-int(b==2))
            if s:expected[key][s<0]+=int(e['weight_units'])
        llrs={}
        for key,(a,o) in expected.items():
            r=classes[key]
            assert int(r['agree_units'])==a and int(r['oppose_units'])==o
            eps=(o+eta)/(a+o+2*eta);llr=max(0,math.log(a+eta)-math.log(o+eta))
            assert math.isclose(float(r['eps']),eps,rel_tol=1e-12,abs_tol=1e-15)
            assert math.isclose(float(r['llr']),llr,rel_tol=1e-12,abs_tol=1e-15)
            llrs[key]=llr
        maximum=max(llrs.values(),default=0)
        for key,llr in llrs.items():
            assert int(classes[key]['confidence'])==(math.floor(cs*llr/maximum+.5) if maximum else 0)
        pooled=Counter()
        for e in factors:
            product=int(e['weight_units'])*int(classes[e['class']]['confidence'])
            v=product if factor_den==cs else product//cs
            if v:pooled[(int(e['target1']),int(e['target2']),int(e['relation']),e['source'])]+=v
        actual=Counter()
        for e in rows(d/'contacts.tsv'):
            actual[(int(e['target1']),int(e['target2']),int(e['relation']),e['source'])]+=int(e['weight_units'])
        assert actual==pooled
    assert int(m['raw_molecule_weight_units'])==1024*int(m['informative_contact_molecules'])
    assert sources['panel']==int(m['panel_weight_units'])
    lam=int(m.get('phase_weight_diagnostic','1'))
    if W and tt:assert wc*W==lam*wg*tt
    else:assert wc==0
    for key in ('hom','het'):assert covered[key]==int(m[f'{key}_targets_satisfied'])
    for key in ('agreeing','opposite','nonheterozygous'):assert outcomes[key]==int(m[f'contact_{key}_units'])
    assert g==int(m['genotype_reward_integer']) and c==int(m['contact_reward_integer'])
    assert g+c==int(m['reward_integer'])<=int(m['segment_upper_integer'])<=upper==int(m['full_graph_upper_integer'])
    if 'genotype_upper_integer' in m:
        ug=int(m['genotype_upper_integer']);uc=int(m['phase_upper_integer'])
        assert ug+uc==upper and g<=ug and c<=uc
        def ceil_ratio(u,l):
            if l==0:return 'NULL' if u else '1.000000'
            q=(u*1000000+l-1)//l;return f'{q//1000000}.{q%1000000:06}'
        assert ceil_ratio(ug,g)==m['genotype_approximation_ratio'] and ceil_ratio(uc,c)==m['phase_approximation_ratio']
        # Mediant inequality: U/L <= max(U_G/G, U_C/C) when both parts are positive.
        if g>0 and c>0:assert upper*g*c<=max(ug*c,uc*g)*(g+c)
    trace=list(rows(d/'pair.trace.tsv'))
    def charge(a,b):
        ordered=(a['haplotype1_donor']!=b['haplotype1_donor'])+(a['haplotype2_donor']!=b['haplotype2_donor'])
        if m.get('recombination_charge','ordered')!='unordered_pair':return ordered
        return min(ordered,(a['haplotype1_donor']!=b['haplotype2_donor'])+(a['haplotype2_donor']!=b['haplotype1_donor']))
    switches=sum(charge(a,b) for a,b in zip(trace,trace[1:]))
    assert switches==int(m['recombinations'])
    if g+c:
        ratio=(upper*1000000+g+c-1)//(g+c);rho=f'{ratio//1000000}.{ratio%1000000:06}'
    else:rho='NULL' if upper else '1.000000'
    assert rho==m['approximation_ratio']
    return dict(status='PASS',genotype_reward=g,phase_reward=c,upper=upper,rho=rho,
                source_weight_units=dict(sources),phase_outcomes=dict(outcomes),recombinations=switches,
                scope='Exact output-factor replay; not evidence of better biological accuracy.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path)
    args=p.parse_args();result=audit(args.directory)
    if args.output:
        with args.output.open('x') as stream:json.dump(result,stream,indent=2)
    print(json.dumps(result,indent=2))
