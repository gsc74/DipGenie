#!/usr/bin/env python3
"""GIAB v5.0q evaluation: hap.py/vcfeval (SNPs) and Truvari refine (SVs). Phasing: phasing.sh (yak).

Use the matching reference, truth VCF and benchmark BED for each class.
The query's callable BED is deliberately not an evaluation mask.
The hap.py query omits records with an allele longer than --max-small-allele
(default 10 kb): such records are outside the small-variant scope and crash
hap.py preprocessing. Truvari uses all query records.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PROTOCOL = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def command(argv, log):
    start = time.monotonic()
    with Path(log).open('a') as f:
        f.write('$ ' + shlex.join(map(str, argv)) + '\n')
        f.flush()
        subprocess.run(list(map(str, argv)), check=True, stdout=f, stderr=f)
    return time.monotonic() - start


def shell(script, log):
    return command(['/bin/bash', '-euo', 'pipefail', '-c', script], log)


def canonical(source, dest, sample, exclude_star=False):
    import pysam
    with pysam.VariantFile(str(source)) as v:
        if list(v.header.samples) != [sample]:
            raise ValueError(f'{source}: expected only sample {sample}')
        with pysam.VariantFile(str(dest), 'wz', header=v.header) as out:
            for record in v:
                if exclude_star and '*' in (record.alts or ()):
                    continue
                record.alleles = tuple(a.upper() for a in record.alleles)
                out.write(record)
    pysam.tabix_index(str(dest), preset='vcf')


def small_only(source, dest, max_allele):
    import pysam
    kept = dropped = 0
    with pysam.VariantFile(str(source)) as v, pysam.VariantFile(str(dest), 'wz', header=v.header) as out:
        for record in v:
            if max(len(a) for a in record.alleles) > max_allele:
                dropped += 1
                continue
            kept += 1
            out.write(record)
    pysam.tabix_index(str(dest), preset='vcf', force=True)
    Path(dest).with_name('query.small.filter.txt').write_text(
        f'max_allele_length={max_allele} kept={kept} dropped={dropped}\n')


def bench_commands(truth, query, ref, bed, out, threads, truvari='truvari'):
    # Keep all SNP/indel context. --refine is split into a separately logged
    # command; no --subset or --write-phab changes the full-domain denominators.
    return [
        [truvari, 'bench', '-b', truth, '-c', query, '-f', ref, '-o', out,
         '--includebed', bed, '--pick', 'ac', '--passonly', '-r', '2000', '-C', '5000'],
        [truvari, 'refine', '--align', 'mafft', '--use-original-vcfs',
         '--threads', str(threads), out],
    ]


def refined_summary(out):
    before = json.loads((out/'summary.json').read_text())
    log = (out/'refine.log.txt').read_text()
    no_regions = '0 regions to be refined' in log and 'No regions to be refined. Exiting' in log
    if '[ERROR]' in log or (not no_regions and ('[CRITICAL]' in log or 'Finished refine' not in log)):
        raise RuntimeError('Truvari refinement did not finish successfully')
    after = before.copy() if no_regions else json.loads((out/'refine.variant_summary.json').read_text())
    for key in ('base cnt', 'comp cnt'):
        if after[key] != before[key]:
            raise RuntimeError('Refinement changed full-domain denominators')
    for total, positive, negative in [('base cnt', 'TP-base', 'FN'), ('comp cnt', 'TP-comp', 'FP')]:
        if after[total] != after[positive] + after[negative]:
            raise RuntimeError('Inconsistent variant counts')
    return dict(before=before, refined=after, no_candidate_regions=no_regions)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for field in ('query', 'reference', 'truth-small', 'truth-sv', 'bed-small', 'bed-sv', 'sdf', 'stratification', 'out'):
        p.add_argument('--'+field, type=Path, required=True)
    p.add_argument('--threads', type=int, default=64)
    p.add_argument('--sample', default='HG002')
    p.add_argument('--gender', choices=('male', 'female'), required=True)
    p.add_argument('--rtg', default='rtg')
    p.add_argument('--truvari', default='truvari')
    p.add_argument('--max-small-allele', type=int, default=10_000)
    a = p.parse_args()
    if a.threads < 1:
        p.error('threads must be positive')
    paths = [a.query, a.reference, a.truth_small, a.truth_sv, a.bed_small, a.bed_sv, a.stratification]
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            p.error(f'missing/empty input: {path}')
    if not a.sdf.is_dir():
        p.error('sdf must be the RTG template directory for the same reference')
    a.out.mkdir(parents=True, exist_ok=False)
    receipt = dict(protocol=PROTOCOL, input_sha256={str(x.resolve()): digest(x) for x in paths},
                   query_callable_mask=False, sample=a.sample, gender=a.gender,
                   source_sha256=digest(__file__), threads=a.threads, max_small_allele=a.max_small_allele,
                   SV_metric='Truvari refined variant-level F1, not genotype-aware Aardvark F1')
    (a.out/'inputs.json').write_text(json.dumps(receipt, indent=2)+'\n')
    for tool, args in [('bcftools', ['--version']), (a.truvari, ['version']),
                       ('mafft', ['--version']), (a.rtg, ['version'])]:
        command([tool, *args], a.out/'versions.log')
    for name, source, exclude in [('truth.small', a.truth_small, False), ('truth.sv', a.truth_sv, True),
                                 ('query', a.query, False), ('query.sv', a.query, True)]:
        canonical(source, a.out/(name+'.vcf.gz'), a.sample, exclude)
    small_only(a.out/'query.vcf.gz', a.out/'query.small.vcf.gz', a.max_small_allele)
    timings = {}
    timings['hap.py'] = command(['hap.py', a.out/'truth.small.vcf.gz', a.out/'query.small.vcf.gz',
        '-r', a.reference, '-f', a.bed_small, '-o', a.out/'happy', '--engine', 'vcfeval',
        '--engine-vcfeval-path', a.rtg, '--engine-vcfeval-template', a.sdf,
        '--gender', a.gender, '--threads', a.threads, '--stratification', a.stratification], a.out/'happy.log')
    for label, cmd in zip(('truvari bench', 'truvari refine'), bench_commands(
            a.out/'truth.sv.vcf.gz', a.out/'query.sv.vcf.gz', a.reference, a.bed_sv, a.out/'truvari', a.threads, a.truvari)):
        timings[label] = command(cmd, a.out/(label.replace(' ', '-')+'.log'))
    sv = refined_summary(a.out/'truvari')
    with (a.out/'happy.summary.csv').open() as f:
        snps = [row for row in csv.DictReader(f) if row['Type']=='SNP' and row['Filter']=='PASS']
    if len(snps) != 1:
        raise RuntimeError('Expected one PASS SNP summary')
    result = dict(SNP_F1=float(snps[0]['METRIC.F1_Score']), SNP_precision=float(snps[0]['METRIC.Precision']),
                  SNP_recall=float(snps[0]['METRIC.Recall']), SNP_counts=snps[0],
                  SV_F1=sv['refined']['f1'], SV_precision=sv['refined']['precision'],
                  SV_recall=sv['refined']['recall'], SV=sv, seconds=timings, protocol=receipt)
    for path, value in receipt['input_sha256'].items():
        if digest(path) != value:
            raise RuntimeError('Input changed during evaluation: '+path)
    (a.out/'metrics.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (a.out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', metrics_sha256=digest(a.out/'metrics.json')))+'\n')
    print(json.dumps({k: result[k] for k in ('SNP_F1', 'SV_F1')}, indent=2))


if __name__ == '__main__':
    main()
