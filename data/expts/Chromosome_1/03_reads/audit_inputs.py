#!/usr/bin/env python3
"""Independently recount nested FASTQs and their same-molecule Hi-C contacts."""
import argparse
import collections
import csv
import hashlib
import itertools
import json
from pathlib import Path
import regions as r


def audit(folder, region_name):
    folder = Path(folder)
    manifest = json.loads((folder / 'PASS.json').read_text())
    seed = manifest['seed']
    assert isinstance(seed, int) and seed >= 0
    rows = manifest['rows']
    memberships = {}
    with (folder / 'membership.tsv').open() as stream:
        for x in csv.DictReader(stream, delimiter='\t'):
            name = x['readID']
            assert name not in memberships
            memberships[name] = (int(x['first_depth_index']), int(x['inside_bases']), int(x['raw_bases']))
    # Independently reconstruct the seeded prefix assignment from the stored
    # per-molecule local aligned-base counts, not from depth folder names.
    length = r.REGIONS[region_name][2] - r.REGIONS[region_name][1]
    ranked = sorted((hashlib.blake2b(f'{seed}:{name}'.encode(), digest_size=16).hexdigest(), name, inside)
                    for name, (_, inside, _) in memberships.items())
    total = sum(x[2] for x in ranked)
    depths = [d for d in (2, 4, 8, 16, 32) if d * length < total]
    assert [x['depth_label'] for x in rows] == [f'{d}x' for d in depths] + ['Total']
    cumulative = 0
    bucket = 0
    for _, name, inside in ranked:
        while bucket < len(depths) and cumulative >= depths[bucket] * length:
            bucket += 1
        assert memberships[name][0] == bucket
        cumulative += inside
    del ranked
    payloads = {}
    for one, two in itertools.zip_longest(r.fastq(folder/'Total/HG002.R1.fastq.gz'),
                                          r.fastq(folder/'Total/HG002.R2.fastq.gz')):
        assert one is not None and two is not None
        assert one[0].endswith('/1') and two[0].endswith('/2') and one[0][:-2] == two[0][:-2]
        name = one[0][:-2]
        assert name in memberships and name not in payloads
        payloads[name] = hashlib.sha256('\n'.join(one + two).encode()).digest()
    assert payloads.keys() == memberships.keys()
    outputs = []
    previous = set()
    previous_contacts = set()
    chrom, lo, hi = r.REGIONS[region_name]
    for i, row in enumerate(rows):
        depth = row['depth_label']
        directory = folder / depth
        expected = {name for name, value in memberships.items() if value[0] <= i}
        seen = set()
        raw_bases = inside_bases = 0
        interleaved = iter(r.fastq(directory / 'HG002.interleaved.fastq.gz'))
        for one, two in itertools.zip_longest(r.fastq(directory / 'HG002.R1.fastq.gz'),
                                              r.fastq(directory / 'HG002.R2.fastq.gz')):
            assert one is not None and two is not None
            assert one[0].endswith('/1') and two[0].endswith('/2') and one[0][:-2] == two[0][:-2]
            name = one[0][:-2]
            assert name in expected and name not in seen
            assert hashlib.sha256('\n'.join(one + two).encode()).digest() == payloads[name]
            seen.add(name)
            assert next(interleaved, None) == one and next(interleaved, None) == two
            bases = len(one[1]) + len(two[1])
            assert bases == memberships[name][2]
            raw_bases += bases
            inside_bases += memberships[name][1]
        assert next(interleaved, None) is None
        assert seen == expected and previous <= seen
        assert len(seen) == row['pairs'] and raw_bases == row['raw_bases'] and inside_bases == row['inside_bases']
        contact_ids = set()
        distances = collections.Counter()
        orientations = collections.Counter()
        side_pairs = collections.Counter()
        for contact in r.pairs_rows(directory / 'contacts.pairs.gz'):
            name = contact['readID']
            assert name in seen and name not in contact_ids
            contact_ids.add(name)
            assert contact['pair_type'] == 'UU'
            positions = []
            for end in (1, 2):
                assert contact[f'chrom{end}'] == chrom and int(contact[f'mapq{end}']) >= 30
                position = int(contact[f'pos{end}']) - 1
                assert lo <= position < hi
                positions.append(position)
            distance = abs(positions[1] - positions[0])
            key = ('<1kb' if distance < 1000 else '1-10kb' if distance < 10000 else
                   '10-100kb' if distance < 100000 else '100kb-1Mb' if distance < 1000000 else '>=1Mb')
            distances[key] += 1
            orientations[contact['strand1'] + contact['strand2']] += 1
            side_pairs[contact['read_side1'] + ',' + contact['read_side2']] += 1
        assert previous_contacts <= contact_ids
        assert len(contact_ids) == row.get('internal_contacts', 0)
        outputs.append(dict(depth=depth, pairs=len(seen), depth_aligned_bases=inside_bases / length,
                            raw_sequence_depth=raw_bases / length, contacts=len(contact_ids),
                            contact_distance_bins=dict(distances), strand_pairs=dict(orientations),
                            parsed_read_sides=dict(side_pairs)))
        previous, previous_contacts = seen, contact_ids
    return dict(status='PASS', region=region_name, seed=seed,
                input_folder=str(folder.resolve()),
                metadata_sha256=hashlib.sha256((folder/'PASS.json').read_bytes()).hexdigest(),
                actual_fastq_nesting=True,
                identical_pair_payload_across_depths=True,
                interleaved_equals_separate_mates=True, seed_prefix_verified=True,
                contact_nesting_and_same_molecule_membership=True, rows=outputs)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('folder', type=Path)
    p.add_argument('region', choices=r.REGIONS)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = audit(args.folder, args.region)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps(result, indent=2))
