#!/usr/bin/env python3

import sys
import os
import argparse

# read prefix and n as arguments
parser = argparse.ArgumentParser(description='Join haplotypes')
parser.add_argument('prefix', type=str, help='prefix of the file')
parser.add_argument('n', type=int, help='number of haplotypes')
args = parser.parse_args()

print(args.prefix)
print(args.n)

# create a GFA with S lines L lines and W lines
# S lines are the sequences with s1 s2 s3 s4 s5 and so on
# prefix_1_1.fa will be read as s1 and sequence will be the node sequence for all _1_1.fa files write s1 ... sN
# prefix_1_2.fa will be read as sN+1 and sequence will be the node sequence for all _1_2.fa files write sN+1 ... s2N
# W lines are walks with w1 and w2 such that w1 is s1 s2 s3 s4 s5 ... SN and w2 is sN+1 sN+2 sN+3 sN+4 sN+5 ... s2N

# create a GFA file
gfa_file = open(args.prefix + ".gfa", "w")
gfa_file.write("H\tVN:Z:1.0\n")


Vertices_1 = []
Vertices_2 = []
# write S lines
# read sequences from all _1_1.fa files and write them as S lines
for i in range(1, args.n + 1):
    filename = args.prefix + "_" + str(i) + "_1.fa"
    with open(filename, "r") as file:
        # skip the first header line
        if file.readline().startswith(">"):
            pass
        sequence = file.read().replace('\n', '')
        gfa_file.write("S\t" + str(i) + "\t" + sequence + "\n")
        Vertices_1.append(i)
        
# read sequences from all _1_2.fa files and write them as S lines
for i in range(1, args.n + 1):
    filename = args.prefix + "_" + str(i) + "_2.fa"
    with open(filename, "r") as file:
        # skip the first header line
        if file.readline().startswith(">"):
            pass
        sequence = file.read().replace('\n', '')
        gfa_file.write("S\t" + str(args.n + i) + "\t" + sequence + "\n")
        Vertices_2.append(args.n + i)
        
# write L lines such that L lines are the links between the sequences
# It looks like L\Vertices_1[i]\t+\Vertices_1[i+1]\t+\t0M for all i from 0 to n-1 
# for Vertices_2 as well

# write L lines
for i in range(args.n - 1):
    gfa_file.write("L\t" + str(Vertices_1[i]) + "\t+\t" + str(Vertices_1[i + 1]) + "\t+\t0M\n")
    gfa_file.write("L\t" + str(Vertices_2[i]) + "\t+\t" + str(Vertices_2[i + 1]) + "\t+\t0M\n")
    
# write W lines
# W lines look like W	hap_1	1	chromosome_id	0	3	>s1>s2>s4>s7>s8 for hap_1
# W hap_2 1 chromosome_id 0 3 >sN+1>sN+2>sN+4>sN+7>sN+8 for hap_2

# write W lines
gfa_file.write("W\thap_1\t1\tchromosome_id\t0\t" + str(args.n) + "\t")
for i in range(args.n):
    gfa_file.write(">" + str(Vertices_1[i]))
gfa_file.write("\n")

gfa_file.write("W\thap_2\t1\tchromosome_id\t0\t" + str(args.n) + "\t")
for i in range(args.n):
    gfa_file.write(">" + str(Vertices_2[i]))
gfa_file.write("\n")