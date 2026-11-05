#!/usr/bin/env python3
"""Select complete primary-alignment pairs, then sample nested whole pairs.

Input SAM must be query-name grouped (the direct paired bwa-mem2 output).
No read trimming, MAPQ threshold, duplicate removal or secondary rescue occurs.
"""
import argparse
import csv
import hashlib
import itertools
import json
import pathlib
import sys

import pysam


def overlap(record, contig, start, end):
    if record.is_unmapped or record.reference_name != contig:
        return 0
    position = record.reference_start
    bases = 0
    for op, length in record.cigartuples or ():
        if op in (0, 7, 8):
            bases += max(0, min(position + length, end) - max(position, start))
        if op in (0, 2, 3, 7, 8):
            position += length
    return bases


def write_fastq(out, record, mate):
    sequence = record.get_forward_sequence()
    qualities = record.get_forward_qualities()
    if not sequence or qualities is None or len(sequence) != len(qualities):
        raise ValueError(f'Missing sequence/qualities: {record.query_name}')
    if any(op == 5 for op, length in record.cigartuples or ()):
        raise ValueError(f'Hard-clipped primary cannot reconstruct full read: {record.query_name}')
    quality = ''.join(chr(q + 33) for q in qualities)
    out.write(f'@{record.query_name}/{mate}\n{sequence}\n+\n{quality}\n')
    return len(sequence)


def extract(args):
    out = pathlib.Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    contig, interval = args.region.split(':')
    first, last = map(int, interval.split('-'))
    start, end = first - 1, last
    stats = dict(all_pairs=0, selected_pairs=0, selected_bases=0,
                 aligned_inside_bases=0, selected_unmapped_mates=0,
                 selected_reverse_mates=0, ignored_nonprimary=0)
    with pysam.AlignmentFile(args.sam, 'r') as sam:
        if contig not in sam.references or end > sam.get_reference_length(contig):
            raise ValueError('MHC region is not contained in reference')
        with pysam.AlignmentFile(str(out/'selected.primary.bam'), 'wb', header=sam.header) as bam, \
                (out/'HG002.R1.fastq').open('w') as r1, (out/'HG002.R2.fastq').open('w') as r2:
            selected_names = set()
            for name, records in itertools.groupby(sam, key=lambda r: r.query_name):
                primary = []
                for record in records:
                    if record.is_secondary or record.is_supplementary:
                        stats['ignored_nonprimary'] += 1
                    else:
                        primary.append(record)
                if not primary:
                    continue
                if len(primary) != 2 or not all(r.is_paired for r in primary):
                    raise ValueError(f'Expected exactly two grouped primary mates: {name}')
                primary.sort(key=lambda r: not r.is_read1)
                if not primary[0].is_read1 or not primary[1].is_read2:
                    raise ValueError(f'Duplicate/missing mate flag: {name}')
                stats['all_pairs'] += 1
                inside = sum(overlap(r, contig, start, end) for r in primary)
                if not inside:
                    continue
                if name in selected_names:
                    raise ValueError(f'Duplicate selected read-pair identity: {name}')
                selected_names.add(name)
                stats['selected_pairs'] += 1
                stats['aligned_inside_bases'] += inside
                for mate, record, stream in zip((1, 2), primary, (r1, r2)):
                    stats['selected_bases'] += write_fastq(stream, record, mate)
                    stats['selected_unmapped_mates'] += int(record.is_unmapped)
                    stats['selected_reverse_mates'] += int(record.is_reverse)
                    bam.write(record)
                if stats['selected_pairs'] % 100000 == 0:
                    print(json.dumps(stats), file=sys.stderr, flush=True)
    if not stats['selected_pairs']:
        raise ValueError('No primary-aligned MHC read pairs selected')
    stats.update(region=args.region, reference_bases=end-start,
                 depth=stats['selected_bases']/(end-start),
                 aligned_depth=stats['aligned_inside_bases']/(end-start),
                 selection='Either primary mate has >=1 aligned M/=X query base inside MHC; retain both original untrimmed mates',
                 mapq_filter='none', source_accession=args.accession)
    (out/'selection.json').write_text(json.dumps(stats, indent=2)+'\n')


def fastq(path):
    with open(path) as stream:
        while name := stream.readline():
            sequence, plus, quality = (stream.readline().rstrip('\r\n') for _ in range(3))
            if not name.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise ValueError(f'Malformed FASTQ: {path}')
            yield name.split()[0], sequence, quality


def pairs(folder):
    for pair in itertools.zip_longest(fastq(folder/'HG002.R1.fastq'), fastq(folder/'HG002.R2.fastq')):
        if None in pair or pair[0][0].removesuffix('/1') != pair[1][0].removesuffix('/2'):
            raise ValueError('FASTQs are not synchronized read pairs')
        yield pair


def sample(args):
    source, out = pathlib.Path(args.input), pathlib.Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    ranked, total_bases = [], 0
    for index, pair in enumerate(pairs(source)):
        name = pair[0][0].removesuffix('/1')
        key = hashlib.blake2b(f'{args.seed}:{name}'.encode(), digest_size=16).digest()
        bases = sum(len(record[1]) for record in pair)
        ranked.append((key, index, bases))
        total_bases += bases
    if not ranked:
        raise ValueError('Empty Total read set')
    ranked.sort()
    depths = [d for d in (2, 4, 8) if d * args.length < total_bases]
    chosen, rows, counters = {}, [], {}
    for depth in depths:
        running, selected = 0, set()
        for key, index, bases in ranked:
            # Nearest total number of bases in one fixed random ordering.
            if abs(running + bases - depth*args.length) > abs(running - depth*args.length):
                break
            selected.add(index)
            running += bases
        if not selected:
            raise ValueError('Target depth selected no pairs')
        chosen[f'{depth}x'] = selected
    chosen['Total'] = None
    streams = {}
    for depth in chosen:
        folder = out/depth
        folder.mkdir()
        streams[depth] = [(folder/name).open('w') for name in ('HG002.R1.fastq','HG002.R2.fastq','HG002.interleaved.fastq')]
        counters[depth] = dict(reads=0, bases=0, identities=set())
    try:
        for index, pair in enumerate(pairs(source)):
            digest = hashlib.sha256(repr(pair).encode()).hexdigest()
            for depth, selected in chosen.items():
                if selected is not None and index not in selected:
                    continue
                for record, stream in zip(pair, streams[depth][:2]):
                    name, sequence, quality = record
                    text = f'{name}\n{sequence}\n+\n{quality}\n'
                    stream.write(text)
                    streams[depth][2].write(text)
                    counters[depth]['reads'] += 1
                    counters[depth]['bases'] += len(sequence)
                counters[depth]['identities'].add(digest)
    finally:
        for group in streams.values():
            for stream in group:
                stream.close()
    previous = set()
    for depth, data in counters.items():
        if not previous <= data['identities']:
            raise ValueError('Nested subset audit failed')
        previous = data['identities']
        if len(data['identities'])*2 != data['reads']:
            raise ValueError('Duplicate pair content')
        row = dict(assay='SR', depth=depth, reads=data['reads'], bases=data['bases'],
                   mean_read_length=data['bases']/data['reads'], coverage=data['bases']/args.length,
                   SRR_accessions=args.accession, reference=args.reference,
                   reference_source='https://github.com/marbl/CHM13',
                   source_identity=args.source_identity, nesting='PASS')
        rows.append(row)
        with (out/depth/'pair-identities.sha256').open('w') as stream:
            stream.write('\n'.join(sorted(data['identities']))+'\n')
    with (out/'coverage.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader(); writer.writerows(rows)
    with (out/'available_depths.tsv').open('w') as stream:
        stream.write('depth\tcoverage\treads\tbases\n')
        for row in rows:
            stream.write(f"{row['depth']}\t{row['coverage']:.9f}\t{row['reads']}\t{row['bases']}\n")
    (out/'sampling.json').write_text(json.dumps(dict(seed=args.seed, method='One fixed hash ordering; nested nearest-base-count prefixes, whole pairs only',
        depth_rule='Include nominal depth only if strictly less than measured Total', reference_bases=args.length,
        source_accession=args.accession, nested='PASS', depths=rows), indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('extract'); p.add_argument('--sam', default='-'); p.add_argument('--region', required=True)
    p.add_argument('--output', required=True); p.add_argument('--accession', required=True); p.set_defaults(func=extract)
    p = sub.add_parser('sample'); p.add_argument('--input', required=True); p.add_argument('--output', required=True)
    p.add_argument('--length', type=int, required=True); p.add_argument('--seed', type=int, default=42)
    p.add_argument('--accession', required=True); p.add_argument('--reference', required=True)
    p.add_argument('--source-identity', required=True); p.set_defaults(func=sample)
    args = parser.parse_args(); args.func(args)


if __name__ == '__main__':
    main()
