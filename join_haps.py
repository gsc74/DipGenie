#!/usr/bin/env python3

import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Join haplotypes')
parser.add_argument('prefix', type=str, help='prefix of the file')
parser.add_argument('n', type=int, help='number of haplotypes')
args = parser.parse_args()

print(args.prefix)
print(args.n)

gfa_file = open(args.prefix + ".gfa", "w")

def chunk_sequence(seq, length=30):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

Vertices_1 = []
Vertices_2 = []
current_node_id = 1

# For haplotype 1 side
for i in range(1, args.n + 1):
    filename = args.prefix + "_" + str(i) + "_1.fa"
    with open(filename, "r") as file:
        file.readline()  # skip header
        sequence = file.read().replace('\n', '')
        chunks = chunk_sequence(sequence, 30)
        node_ids = []
        for c in chunks:
            gfa_file.write("S\t" + str(current_node_id) + "\t" + c + "\n")
            node_ids.append(current_node_id)
            current_node_id += 1
        Vertices_1.append(node_ids)

# For haplotype 2 side
for i in range(1, args.n + 1):
    filename = args.prefix + "_" + str(i) + "_2.fa"
    with open(filename, "r") as file:
        file.readline()  # skip header
        sequence = file.read().replace('\n', '')
        chunks = chunk_sequence(sequence, 30)
        node_ids = []
        for c in chunks:
            gfa_file.write("S\t" + str(current_node_id) + "\t" + c + "\n")
            node_ids.append(current_node_id)
            current_node_id += 1
        Vertices_2.append(node_ids)

# Merge vertices
merged_vertices_1 = []
for v_list in Vertices_1:
    for v in v_list:
        merged_vertices_1.append(v)
merged_vertices_2 = []
for v_list in Vertices_2:
    for v in v_list:
        merged_vertices_2.append(v)


# Write L lines linking consecutive nodes within each haplotype path
for v in range(len(merged_vertices_1)-1):
    gfa_file.write("L\t" + str(merged_vertices_1[v]) + "\t+\t" + str(merged_vertices_1[v+1]) + "\t+\t0M\n")
    
for v in range(len(merged_vertices_2)-1):
    gfa_file.write("L\t" + str(merged_vertices_2[v]) + "\t+\t" + str(merged_vertices_2[v+1]) + "\t+\t0M\n")

# Write W lines for the two haplotypes
# Now these will include all the chunked nodes
# hap_1
gfa_file.write("W\thap_1\t1\tchromosome_id\t0\t")
gfa_file.write(str(5E6) + "\t")
for v in merged_vertices_1:
    gfa_file.write(">" + str(v))
gfa_file.write("\n")

# hap_2
gfa_file.write("W\thap_2\t1\tchromosome_id\t0\t")
gfa_file.write(str(5E6) + "\t")
for v in merged_vertices_2:
    gfa_file.write(">" + str(v))
gfa_file.write("\n")

gfa_file.close()