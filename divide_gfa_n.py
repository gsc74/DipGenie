#!/usr/bin/env python3

import sys
import networkx as nx

def read_gfa(filename):
    s_lines, l_links, w_lines = [], [], []
    haps, len_v_hap, all_haps = {}, {}, {}
    adj, in_degree, out_degree, node_length = {}, {}, {}, {}
    s_map = {}

    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            if line.startswith('S'):
                s_lines.append(parts)
            elif line.startswith('L'):
                l_links.append(parts)
            elif line.startswith('W'):
                w_lines.append(parts)

    for s_line in s_lines:
        s_map[s_line[1]] = s_line[2]
        node_length[s_line[1]] = len(s_line[2])

    for l_line in l_links:
        v1, v2 = l_line[1], l_line[3]
        adj.setdefault(v1, set()).add(v2)
        adj.setdefault(v2, set())

    # Compute topological order
    G = nx.DiGraph()
    for v in adj:
        for u in adj[v]:
            G.add_edge(v, u)
    vertices_ = list(nx.topological_sort(G))

    for v in adj:
        in_degree.setdefault(v, 0)
        out_degree.setdefault(v, 0)
        for u in adj[v]:
            in_degree[u] = in_degree.get(u, 0) + 1
            out_degree[v] = out_degree.get(v, 0) + 1

    haplotypes = []
    for w_line in w_lines:
        haplotype = w_line[1] + '.' + w_line[2]
        haplotypes.append(haplotype)
        path = [x for x in w_line[6].split('>') if x]
        for x in path:
            haps.setdefault(x, set()).add(haplotype)

    # Compute lengths per haplotype at each vertex
    for w_line in w_lines:
        haplotype = w_line[1] + '.' + w_line[2]
        path = [x for x in w_line[6].split('>') if x]
        len_v_hap[haplotype] = {}
        for v in vertices_:
            len_v_hap[haplotype][v] = 0
        curr_len = 0
        for v in path:
            len_v_hap[haplotype][v] = curr_len
            curr_len += node_length[v]

    for v in haps:
        all_haps[v] = (len(haps[v]) == len(w_lines) and in_degree[v] != 0 and out_degree[v] != 0)

    return s_lines, l_links, w_lines, haps, len_v_hap, all_haps, haplotypes, vertices_, s_map, node_length

def find_break_vertices(vertices, hap, len_v_hap, all_haps, node_length, n):
    # Find total length of haplotype path
    max_len = 0
    for v in vertices:
        val = len_v_hap[hap][v] + node_length[v]
        if val > max_len:
            max_len = val

    # We want n parts, so (n-1) breakpoints at multiples of max_len/n
    segment_length = max_len / float(n)
    targets = [segment_length * i for i in range(1, n)]  # from 1/n to (n-1)/n

    break_vertices = []
    current_target_index = 0

    # Iterate over vertices in topological order to find suitable break vertices
    for v in vertices:
        start_pos = len_v_hap[hap][v]
        end_pos = start_pos + node_length[v]

        # Check if we cross the next target(s)
        while current_target_index < len(targets) and end_pos >= targets[current_target_index]:
            if all_haps.get(v, False):
                break_vertices.append(v)
                current_target_index += 1
            else:
                # If this vertex doesn't satisfy all_haps[v], we do not increment target.
                # Just break from the while loop, move to the next vertex and hope to find
                # a suitable vertex for this same target later.
                break

        if current_target_index == len(targets):
            # Found all required break vertices
            break

    return break_vertices, max_len

def extract_segment(w_lines, start_vertex, end_vertex):
    # Extract segment of W lines between start_vertex and end_vertex
    # If start_vertex is None, we start from beginning of path.
    # If end_vertex is None, we go till the end of the path.
    new_w_lines = []
    for w_line in w_lines:
        parts = w_line[:6]
        path = w_line[6].split('>')[1:]
        subpath = []
        recording = (start_vertex is None)
        for v in path:
            if v == start_vertex:
                recording = True
                continue
            if recording:
                subpath.append(v)
            if v == end_vertex:
                break
        if subpath:
            new_w_lines.append(parts + [">" + ">".join(subpath)])
    return new_w_lines

def build_gfa(s_map, segment_w_lines):
    # Given truncated W lines, build the set of nodes and links
    vtx_set = set()
    new_l_links = set()
    for w_line in segment_w_lines:
        path = w_line[6].split('>')[1:]
        vtx_set.update(path)
        for i in range(len(path)-1):
            line_l = f"L\t{path[i]}\t+\t{path[i+1]}\t+\t0M"
            new_l_links.add(line_l)
    s_lines_out = []
    for vtx in vtx_set:
        s_lines_out.append(f"S\t{vtx}\t{s_map[vtx]}")
    return s_lines_out, new_l_links, segment_w_lines

def write_gfa_file(output_file, s_lines_out, l_links_out, w_lines_out):
    with open(output_file, 'w') as f:
        for line in s_lines_out:
            f.write(line + "\n")
        for ll in l_links_out:
            f.write(ll + "\n")
        for wl in w_lines_out:
            f.write("\t".join(wl) + "\n")

if __name__ == "__main__":
    # Usage: python script.py input.gfa n output_prefix
    # Divides the graph into n parts based on haplotype 0.
    if len(sys.argv) != 4:
        print("Usage: python script.py input.gfa n output_prefix")
        sys.exit(1)

    input_gfa = sys.argv[1]
    n = int(sys.argv[2])
    output_prefix = sys.argv[3]

    print(f"Input GFA file: {input_gfa}")
    print(f"Chop into {n} equal parts based on haplotype 0.")

    (s_lines, l_links, w_lines,
     haps, len_v_hap, all_haps, haplotypes, vertices, s_map, node_length) = read_gfa(input_gfa)

    # We'll pick the first haplotype
    h = haplotypes[0]

    break_vertices, total_length = find_break_vertices(vertices, h, len_v_hap, all_haps, node_length, n)
    print("Break vertices:", break_vertices)

    # Segments: from start(None) to break_vertices[0], break_vertices[0] to break_vertices[1], etc.
    segment_boundaries = [None] + break_vertices + [None]

    num_segments = len(segment_boundaries)-1
    print(f"Will produce {num_segments} GFA files (segments).")

    for i in range(num_segments):
        start_v = segment_boundaries[i]
        end_v = segment_boundaries[i+1]
        segment_w_lines = extract_segment(w_lines, start_v, end_v)
        s_lines_out, l_links_out, w_lines_out = build_gfa(s_map, segment_w_lines)
        output_file = f"{output_prefix}_{i+1}.gfa"
        write_gfa_file(output_file, s_lines_out, l_links_out, w_lines_out)

    print(f"Generated {num_segments} GFA files.")