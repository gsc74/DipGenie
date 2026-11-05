#!/usr/bin/env python3
"""Re-select complete chr1 SR pairs and make nested, aligned-base-depth sets."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import pysam
from regions import REGIONS, primary_pair, query_intervals, union_length, fastq, nesting_members, pair_rank, Compressed



def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''): h.update(block)
    return h.hexdigest()


def prepare(source, out, seed):
    pair_rank(seed, 'validate')
    assert not out.exists(), 'Refuse to overwrite inputs'
    selection = json.loads((source / 'extracted/selection.json').read_text())
    assert selection['source_accession'] == 'SRR11321732'
    assert selection['region'] == 'chr1:1-248387328'
    command = (source / 'extraction.command').read_text()
    assert 'bwa-mem2 mem' in command and 'SRR11321732' in command and 'chm13v2.0.fa' in command
    region = REGIONS['chr1']; length = region[2] - region[1]
    bam = source / 'extracted/selected.primary.bam'
    accepted = {}; ranked = []; seen = set(); rejected = 0
    with pysam.AlignmentFile(str(bam), 'rb') as stream:
        assert stream.get_reference_length('chr1') == length
        for name, group in itertools.groupby(stream, key=lambda r: r.query_name):
            assert name not in seen, 'Non-name-grouped or duplicate primary records'
            seen.add(name); pair = primary_pair(list(group))
            if pair is None or not all(not r.is_unmapped and r.reference_name == 'chr1' and
                                      r.reference_start >= 0 and r.reference_end <= length and query_intervals(r, region) for r in pair):
                rejected += 1; continue
            inside = sum(union_length(query_intervals(r, region)) for r in pair)
            raw = sum(len(r.query_sequence) for r in pair)
            accepted[name] = (inside, raw)
            ranked.append((pair_rank(seed, name), name, inside, raw))
            if len(seen) % 1000000 == 0: print('primary_pairs', len(seen), 'accepted', len(accepted), flush=True)
    del seen
    depths, buckets, total = nesting_members(ranked, length, [2, 4, 8, 16, 32])
    del ranked
    assert depths == [2, 4, 8, 16, 32], 'Insufficient aligned depth for six-set study'
    labels = [f'{d}x' for d in depths] + ['Total']
    out.mkdir(parents=True)
    handles = {}; stats = {label: dict(pairs=0, inside_bases=0, raw_bases=0, internal_contacts=0) for label in labels}
    for label in labels:
        folder = out / label; folder.mkdir()
        handles[label] = [Compressed(folder / f'HG002.R{i}.fastq.gz') for i in (1, 2)] + [Compressed(folder / 'HG002.interleaved.fastq.gz')]
        # SR mates supply linkage directly; no Hi-C contact annotation.
        contacts = Compressed(folder / 'contacts.pairs.gz')
        contacts.write('## pairs format v1.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2 read_side1 read_side2\n')
        contacts.close()
    written = 0
    with (out / 'membership.tsv').open('x') as members:
        members.write('readID\tfirst_depth_index\tinside_bases\traw_bases\n')
        for a, b in itertools.zip_longest(fastq(source / 'extracted/HG002.R1.fastq'), fastq(source / 'extracted/HG002.R2.fastq')):
            assert a is not None and b is not None
            assert a[0].endswith('/1') and b[0].endswith('/2') and a[0][:-2] == b[0][:-2]
            name = a[0][:-2]
            if name not in accepted: continue
            inside, raw = accepted.pop(name); bucket = buckets.pop(name)
            assert raw == len(a[1]) + len(b[1])
            records = [f'@{r[0]}\n{r[1]}\n+\n{r[2]}\n' for r in (a, b)]
            members.write(f'{name}\t{bucket}\t{inside}\t{raw}\n')
            for label in labels[bucket:]:
                for mate in (0, 1): handles[label][mate].write(records[mate])
                handles[label][2].write(''.join(records))
                stats[label]['pairs'] += 1; stats[label]['inside_bases'] += inside; stats[label]['raw_bases'] += raw
            written += 1
            if written % 1000000 == 0: print('written_pairs', written, flush=True)
    for writers in handles.values():
        for writer in writers: writer.close()
    assert not accepted and not buckets and stats['Total']['inside_bases'] == total
    rows = [dict(region='chr1', depth_label=label, **stats[label], reference_bases=length,
                 depth=stats[label]['inside_bases']/length, raw_sequence_depth=stats[label]['raw_bases']/length) for label in labels]
    with (out / 'coverage.csv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    receipt = dict(status='PASS', library='SR', accession='SRR11321732', seed=seed,
                   nested_by_pair_hash_prefix=True, paired_ids_verified=written, excluded_noninternal_pairs=rejected,
                   coverage_basis='sum of aligned primary M/=X query bases, divided by chr1 length; mate gap never counted',
                   pair_selection='both primary mate spans fully contained in chr1; no outside/unmapped-mate rescue',
                   contact_policy='SR linkage comes from paired FASTQs; empty contacts files are not passed to inference',
                   unavailable_depths=[], rows=rows,
                   source_files={str(p): sha(p) for p in (source/'extracted/selection.json', source/'extraction.command', bam)},
                   source_FASTQs=[str(source/f'extracted/HG002.R{i}.fastq') for i in (1, 2)])
    with (out / 'PASS.json').open('x') as stream: json.dump(receipt, stream, indent=2)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='output directory of align_sr.sh')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output, args.seed)['rows'], indent=2))
