#!/usr/bin/env python3
"""Recombination budget R for one chromosome (DipGenie -R).

Expected number of switches under the Li and Stephens copying model, doubled:
    E[R] = (L - 1) * (H - 1)/H * (1 - exp(-4 * Ne * r / H)),   R = 2 * round(E[R])

Example (human chromosome 1, CHM13):
    python3 scripts/compute_R.py --length 248387328 --rate 1.12e-8    # prints 2140
"""
import argparse
import math

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--length', type=int, required=True, help='chromosome length L in bp')
parser.add_argument('--rate', type=float, required=True, help='recombination rate r per bp (e.g. 1.12e-8)')
parser.add_argument('--haplotypes', type=int, default=51,
                    help='panel haplotypes H used by DipGenie (default 51: the reference and 25 samples)')
parser.add_argument('--ne', type=float, default=5000, help='effective population size Ne (default 5000)')
a = parser.parse_args()
if a.length < 2 or a.rate <= 0 or a.haplotypes < 2 or a.ne <= 0:
    parser.error('length >= 2, rate > 0, haplotypes >= 2 and ne > 0 are required')
expected = (a.length - 1) * (a.haplotypes - 1) / a.haplotypes * (1 - math.exp(-4 * a.ne * a.rate / a.haplotypes))
print(2 * math.floor(expected + 0.5))
