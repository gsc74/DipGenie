#include "solver.h"
//#include "edlib.h"
#include "misc.h"
#include "Classifier.hpp"
#include "Fitter.hpp"


// Constructor
Solver::Solver(gfa_t *g) {
    this->g = g;
}

KSEQ_INIT(gzFile, gzread)


uint64_t hash128_to_64_(const std::string &str) {
    uint64_t hash_output[2]; // To store the 128-bit output

    // Hash the string using MurmurHash3_x64_128
    MurmurHash3_x64_128(str.data(), str.length(), 0, hash_output);

    // Combine the two 64-bit parts of the 128-bit hash into a single 64-bit hash
    return hash_output[0] ^ hash_output[1]; // Simple XOR combination
}


void::Solver::read_gfa() 
{
    uint v;
    n_vtx = gfa_n_vtx(g);
    // Resize node_len
    node_len.resize(n_vtx/2, 0);
    node_seq.resize(n_vtx/2, "");
    std::vector<std::vector<uint32_t>> adj_;
    adj_.resize(n_vtx); 
    /* Node_len */
    for (int v = 0; v < n_vtx/2; v++)
    {
        gfa_seg_t segment = (g)->seg[v];
        int len =  segment.len;
        node_len[v] = len;
        node_seq[v] = segment.seq;
    }
    // look for all the edges , if the sum of all the 
    // edges are zero then, that's a linear reference
    u_int num_edges = 0;
    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        int v_ = av->v_lv >> 32;
        // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::std::endl; 
        for (int i = 0; i < n_edges; i++)
        {
            num_edges++;
        }
    }


    for (v = 0; v < n_vtx; v++)
    {
        gfa_arc_t *av = gfa_arc_a(g, v);
        int n_edges = gfa_arc_n(g, v);
        if (num_edges == 0) // Linear Reference
        {
            lin_ref = 1; // Mark as a linear reference
        }else
        {
            int v_ = av->v_lv >> 32;
            // std::cerr << " node : " << v << " node_len : " << node_len[v] << std::std::endl; 
            for (int i = 0; i < n_edges; i++)
            {
                uint w = av[i].w;
                adj_[v_].push_back(w);
            }
        }
    }

    // Copy the adjacency list only for forward strand
    n_vtx /= 2;
    adj_list.resize(n_vtx);
    for (int i = 0; i < 2 * n_vtx; i++)
    {
        if (i % 2 == 0) // Only forward strand
        {
            for (auto v : adj_[i])
            {
                adj_list[i/2].push_back(v/2);
            }
        }
    }

    adj_.clear();

    // Read walks
    if (g->n_walk > 0) is_gfa_v12 = true;
    
    num_walks = g->n_walk; // count number of walks
    
    //num_walks = 10;


    haps.resize(n_vtx);
    paths.resize(num_walks);
    in_paths.resize(num_walks, std::vector<int32_t>(n_vtx, 0));

    // Fill the paths
    for (size_t w = 0; w < num_walks; w++)
    {
        std::string walk_name = std::string(g->walk[w].sample) + "." + std::to_string(g->walk[w].hap);
        hap_id2name.push_back(walk_name);
        int32_t idx = 0;
        for (size_t n = 0; n < g->walk[w].n_v; n++)
        {
            int v = g->walk[w].v[n];
            if (v%2 != 0) {
                //fprintf(stderr, "Error: Walk %d has reverse strand vertices %d\n", w, v);
                exit(1);
            }
			v /= 2; // for forward strand
            haps[v].push_back(w); // Map forward walk to haplotype
            paths[w].push_back(v); // Map forward walk to path
            in_paths[w][v] = 1;
        }
    }

    const int32_t V = n_vtx; // forward-strand vertex count (already halved above)

    // 1) Seed with earliest position on any walk.
    //    Use int64_t to be safe on very long walks.
    const int64_t INF = std::numeric_limits<int64_t>::max()/4;
    std::vector<int64_t> pos(V, INF);

    for (size_t w = 0; w < paths.size(); ++w) {
        const auto &pw = paths[w];
        for (int64_t t = 0; t < (int64_t)pw.size(); ++t) {
            int32_t vtx = pw[t];
            if (vtx < 0 || vtx >= V) continue;
            if (t < pos[vtx]) pos[vtx] = t; // earliest column where it appears
        }
    }

    // Any vertex that never appears in any walk: give it a stable default.
    // We’ll drop them at the end or park them after all walked columns.
    // Here we set to a large bucket that will be densified later.
    int64_t fallback_start = 0;
    {
        // If we have any seeded positions, set fallback to max seeded + 1; else 0
        int64_t max_seed = -1;
        for (int32_t v = 0; v < V; ++v) if (pos[v] != INF) max_seed = std::max(max_seed, pos[v]);
        fallback_start = (max_seed >= 0 ? max_seed + 1 : 0);
        for (int32_t v = 0; v < V; ++v) if (pos[v] == INF) pos[v] = fallback_start; // park them after seen columns
    }

    // 2) Enforce monotonicity along each walk:
    //    For each consecutive pair (u -> v) on a walk, require pos[v] >= pos[u] + 1.
    //    Iterate until convergence (cap iterations to V to avoid pathological cycles).
    bool changed = true;
    int iter = 0, iter_cap = std::max(10, V); // generous cap; walks are linear so this converges fast
    while (changed && iter++ < iter_cap) {
        changed = false;
        for (size_t w = 0; w < paths.size(); ++w) {
            const auto &pw = paths[w];
            for (size_t t = 1; t < pw.size(); ++t) {
                int32_t u = pw[t-1], v = pw[t];
                if (u < 0 || u >= V || v < 0 || v >= V) continue;
                int64_t need = pos[u] + 1;
                if (pos[v] < need) { pos[v] = need; changed = true; }
            }
        }
    }

    // 3) Densify columns: map (possibly gappy) pos[] to a compact 0..K-1.
    std::vector<std::pair<int64_t,int32_t>> by_pos;
    by_pos.reserve(V);
    for (int32_t v = 0; v < V; ++v) by_pos.emplace_back(pos[v], v);
    std::sort(by_pos.begin(), by_pos.end()); // primary by column, then vertex id

    // Assign dense ranks
    std::vector<int32_t> dense_pos(V, -1);
    int32_t cur_rank = -1;
    int64_t prev_col = std::numeric_limits<int64_t>::min();
    for (auto &p : by_pos) {
        int64_t col = p.first;
        int32_t v = p.second;
        if (col != prev_col) { ++cur_rank; prev_col = col; }
        dense_pos[v] = cur_rank;
    }
    int32_t K = cur_rank + 1; // total number of MSA-like columns that contain at least one vertex

    // Build the final "top_order" replacement using these dense positions.
    top_order.clear();
    top_order.reserve(V);
    for (auto &p : by_pos) top_order.push_back(p.second);

    // Build map: vertex -> global order index (0..V-1), i.e., rank in top_order
    top_order_map.assign(V, -1);
    for (int32_t i = 0; i < (int32_t)top_order.size(); ++i)
        top_order_map[top_order[i]] = i;


    struct SkipEdge { int32_t v_from, v_to; int32_t col_from, col_to; int32_t skip_len; };
    std::vector<SkipEdge> skip_edges; skip_edges.reserve(V);

    for (size_t w = 0; w < paths.size(); ++w) {
        const auto &pw = paths[w];
        for (size_t t = 1; t < pw.size(); ++t) {
            int32_t u = pw[t-1], v = pw[t];
            int32_t cu = dense_pos[u], cv = dense_pos[v];
            if (cu >= 0 && cv >= 0 && cv > cu + 1) {
                skip_edges.push_back({u, v, cu, cv, cv - cu - 1});
            }
        }
    }

    for (int32_t u = 0; u < V; ++u) {
        auto &nei = adj_list[u];
        std::sort(nei.begin(), nei.end(), [&](int32_t a, int32_t b){
            int32_t ca = dense_pos[a], cb = dense_pos[b];
            if (ca != cb) return ca < cb;
            return a < b;
        });
    }

    // From this point on, treat "top_order" and "top_order_map" as the MSA-like order.
    // top_order is sorted by global column, and top_order_map[v] gives its order index.
}

// Read the reads
void Solver::read_ip_reads(std::vector<std::pair<std::string, std::string>> &ip_reads, std::string ip_reads_file)
{
    // Read the reads with gzip
    gzFile fp;
    kseq_t *seq;
    int32_t l;

    fp = gzopen(ip_reads_file.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        ip_reads.push_back(std::make_pair(seq->name.s, seq->seq.s));
    }

    kseq_destroy(seq);
    gzclose(fp);
}


std::string reverse_strand(std::string seq)
{
    std::string rev_seq = "";
    for (int i = seq.size() - 1; i >= 0; i--)
    {
        if (seq[i] == 'A' || seq[i] == 'a')
        {
            rev_seq += 'T';
        }
        else if (seq[i] == 'T' || seq[i] == 't')
        {
            rev_seq += 'A';
        }
        else if (seq[i] == 'C' || seq[i] == 'c')
        {
            rev_seq += 'G';
        }
        else if (seq[i] == 'G' || seq[i] == 'g')
        {
            rev_seq += 'C';
        }else
        {
            rev_seq += seq[i];
        }
        
    }
    return rev_seq;
}

std::vector<std::pair<uint64_t, Anchor>> Solver::index_kmers(int32_t hap) {

    std::vector<std::pair<uint64_t, Anchor>> kmer_index;
    std::string haplotype;

    // Concatenate sequences in the path to form the haplotype
    for (size_t i = 0; i < paths[hap].size(); i++) {
        haplotype += node_seq[paths[hap][i]];
    }

    // Convert to uppercase
    std::transform(haplotype.begin(), haplotype.end(), haplotype.begin(), ::toupper);
    size_t hap_size = haplotype.size();

    if (hap_size < window + k_mer - 1) return kmer_index;

    // Precompute index to vertex map
    std::vector<int32_t> idx_vtx_map(hap_size, -1);
    int32_t count_idx = 0;
    for (size_t i = 0; i < paths[hap].size(); i++) {
        for (size_t j = 0; j < node_seq[paths[hap][i]].size(); j++) {
            idx_vtx_map[count_idx++] = paths[hap][i]; // base_idx -> vertex_idx
        }
    }

    uint64_t prev_hash = UINT64_MAX; // Initialize to the maximum possible value

    // Deque to store k-mers and their positions
    std::deque<std::pair<std::string, int32_t>> window_deque;

    for (int32_t i = 0; i <= hap_size - k_mer; ++i) {
        // Extract the forward and reverse k-mers
        std::string fwd_kmer = haplotype.substr(i, k_mer);

        std::string rev_kmer = reverse_strand_(fwd_kmer);
        // Choose the lexicographically smallest k-mer
        std::string min_kmer = std::min(fwd_kmer, rev_kmer);

        // Remove elements from the back of the deque if they are larger than the current k-mer
        while (!window_deque.empty() && window_deque.back().first >= min_kmer) {
            window_deque.pop_back();
        }

        // Add current element at the back of the deque
        window_deque.emplace_back(min_kmer, i);

        // Remove elements from the front if they are out of the window
        if (!window_deque.empty() && window_deque.front().second <= i - window) {
            window_deque.pop_front();
        }

        // Process the minimizer once the window is fully populated
        if (i >= window - 1) {  
            std::string best_kmer = window_deque.front().first;
            //different seed from earlier
            uint64_t best_hash = hash128_to_64_(best_kmer.c_str());
            if (best_hash != prev_hash) {
                prev_hash = best_hash;

                // Create anchor
                Anchor anchor;
                anchor.h = hap;
                std::unordered_set<int32_t> unique_vtx_set;
                std::vector<int32_t> unique_vtxs_vec;

                // Map k-mer positions to vertex indices
                int32_t best_start_idx = window_deque.front().second;
                for (int32_t j = best_start_idx; j < best_start_idx + k_mer; j++) {
                    int32_t vtx_idx = idx_vtx_map[j];
                    if (unique_vtx_set.find(vtx_idx) == unique_vtx_set.end()) {
                        unique_vtx_set.insert(vtx_idx);
                        unique_vtxs_vec.push_back(vtx_idx);
                    }
                }
                // sort this by topological order
                std::sort(unique_vtxs_vec.begin(), unique_vtxs_vec.end(), [&](int32_t a, int32_t b) {
                    return top_order_map[a] < top_order_map[b];
                });

                // Assign unique vertices to the anchor
                anchor.k_mers = std::move(unique_vtxs_vec);
                kmer_index.emplace_back(best_hash, anchor);
            }
        }
    }
    return kmer_index;
}


std::set<uint64_t> Solver::compute_hashes(std::string &read_seq) {
    // Convert to uppercase
    std::transform(read_seq.begin(), read_seq.end(), read_seq.begin(), ::toupper);

    std::set<uint64_t> read_hashes;
    int32_t seq_size = read_seq.size();
    if (seq_size < window + k_mer - 1) return read_hashes;

    uint64_t prev_hash = UINT64_MAX; // Initialize to maximum possible value

    // Deque to store k-mers and their positions
    std::deque<std::pair<std::string, int32_t>> window_deque;

    for (int32_t i = 0; i <= seq_size - k_mer; ++i) {
        // Extract the forward and reverse k-mers
        std::string fwd_kmer = read_seq.substr(i, k_mer);
        std::string rev_kmer = reverse_strand_(fwd_kmer);

        // Choose the lexicographically smallest k-mer
        std::string min_kmer = std::min(fwd_kmer, rev_kmer);

        // Remove elements from the back of the deque if they are larger than the current k-mer
        while (!window_deque.empty() && window_deque.back().first >= min_kmer) {
            window_deque.pop_back();
        }

        // Add current element at the back of the deque
        window_deque.emplace_back(min_kmer, i);

        // Remove elements from the front if they are out of the window
        if (!window_deque.empty() && window_deque.front().second <= i - window) {
            window_deque.pop_front();
        }

        // Only insert the hash if it is different from the previous hash
        if (i >= window - 1) {  // Once we reach the window size
            std::string best_kmer = window_deque.front().first;
            uint64_t best_hash = hash128_to_64_(best_kmer.c_str());
            if (best_hash != prev_hash) {
                read_hashes.insert(best_hash);
                prev_hash = best_hash;
            }
        }
    }

    return read_hashes;
}


std::vector<std::vector<std::vector<int32_t>>> Solver::compute_anchors(std::vector<std::pair<uint64_t, Anchor>> &minimizers, std::map<uint64_t, int32_t> &read_hashes)
{
    std::vector<std::vector<std::vector<int32_t>>> anchors;
    std::vector<std::vector<std::pair<int32_t, std::vector<int32_t>>>> local_anchors(num_threads);
    #pragma omp parallel for num_threads(num_threads)
    for (int64_t i = 0; i < minimizers.size(); i++)
    {
        int32_t tid = omp_get_thread_num();
        auto minimizer  = minimizers[i];
        auto hash = minimizer.first;
        if (read_hashes.find(hash) != read_hashes.end()) // Found a match
        {
            std::vector<int32_t> anchor;
            for (size_t j = 0; j < minimizer.second.k_mers.size(); j++)
            {
                anchor.push_back(minimizer.second.k_mers[j]);
            }
            local_anchors[tid].push_back(std::make_pair(read_hashes[hash], anchor)); // id, anchor
        }
    }

    anchors.resize(read_hashes.size());
    for (int32_t i = 0; i < num_threads; i++)
    {
        for (auto anchor: local_anchors[i])
        {
            anchors[anchor.first].push_back(anchor.second);  // (id -> anchor_path)
        }
    }
    local_anchors.clear();
    return anchors;
}


void Solver::compute_and_classify_anchors(std::vector<std::pair<std::string, std::string>> &full_ip_reads)
{

    std::vector<int32_t> hap_sizes(num_walks);

    #pragma omp parallel for num_threads(num_threads)
    for (size_t h = 0; h < num_walks; h++)
    {
        std::string haplotype = "";
        for (size_t i = 0; i < paths[h].size(); i++) haplotype += node_seq[paths[h][i]];
        hap_sizes[h] = haplotype.size();
    }

    
    auto ip_reads = full_ip_reads;
    int32_t num_reads = ip_reads.size();

    // Index the kmers
    std::cerr << "Number of Minimizers" << std::endl;
    std::vector<std::vector<std::pair<uint64_t, Anchor>>> kmer_index(num_walks);
    std::vector<std::unordered_map<uint64_t, std::vector<int32_t>>> kmer_hist_per_hap(num_walks);
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t h = 0; h < num_walks; h++) {
        kmer_index[h] = index_kmers(h);
        std::string hap_name = hap_id2name[h];
        fprintf(stderr, "%s : %d\n", hap_name.c_str(), kmer_index[h].size());

        auto& kmer_hist = kmer_hist_per_hap[h];  // Get the map for the current haplotype

        for (const auto& minimizer : kmer_index[h]) {
            // Check if the kmer is present in the current haplotype's kmer_hist
            if (kmer_hist.find(minimizer.first) == kmer_hist.end() || (std::find(kmer_hist[minimizer.first].begin(), kmer_hist[minimizer.first].end(), h) == kmer_hist[minimizer.first].end())) {
                kmer_hist[minimizer.first].push_back(h);
            }
        }
    }

    // Merge kmer_hist_per_hap into a single map
    std::unordered_map<uint64_t, std::vector<int32_t>> kmer_hist;
    std::unordered_map<uint64_t, int32_t> uniqe_kmers;
    for (const auto& hap_hist : kmer_hist_per_hap) {
        for (const auto& entry : hap_hist) {
            auto& global_vector = kmer_hist[entry.first];
            global_vector.insert(global_vector.end(), entry.second.begin(), entry.second.end());
            std::sort(global_vector.begin(), global_vector.end());
            global_vector.erase(std::unique(global_vector.begin(), global_vector.end()), global_vector.end());
            uniqe_kmers[entry.first]++;
        }
    }

    std::vector<int32_t> kmer_hist_count(num_walks + 1, 0);
    for (const auto& kmer : kmer_hist) {
        kmer_hist_count[kmer.second.size()]++;
    }

    if (debug)
    {
        // plot the histogram [x, y] -> [number of haplotypes, number of kmers] with **** on the screen with two for loops
        // print shared fraction of unique kmers by haplotypes
        //fprintf(stderr, "Shared fraction of unique kmers by haplotypes\n");
        int32_t count_unique_kmers = uniqe_kmers.size();
        int32_t sum_unique_kmers = 0;
        for (int32_t i = 1; i <= num_walks; i++) {
            // fprintf(stderr, "[%d, %.5f]\n", i, (float)kmer_hist_count[i]/(float)count_unique_kmers);
            //fprintf(stderr, "[Haplotypes: %d, fraction of unique shared kmers: %.5f]\n", i, (float)kmer_hist_count[i]/(float)count_unique_kmers);
            sum_unique_kmers += kmer_hist_count[i];
        }
        assert(sum_unique_kmers == count_unique_kmers);
    }
    // free unnecessary memory
    kmer_hist.clear();
    kmer_hist_per_hap.clear();
    kmer_hist_count.clear();
    //fprintf(stderr, "[M::%s::%.3f*%.2f] Haplotypes sketched\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0));

    // Compute the anchors
    int64_t num_kmers = 0;
    std::vector<std::set<uint64_t>> Read_hashes(num_reads);
    std::map<uint64_t, int32_t> Sp_R;
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t r = 0; r < num_reads; r++)
    {
        Read_hashes[r] = compute_hashes(ip_reads[r].second);
    }
    // push all the unique read hashes to a set
    for (int32_t r = 0; r < num_reads; r++)
    {
        for (auto hash: Read_hashes[r])
        {
            Sp_R[hash]++; // duplicate hashes are not allowed
        }
    }
    // reset Sp_R values from 0 to max_Sp_R
    count_sp_r = 0;
    for (auto &hash: Sp_R)
    {
        hash.second = count_sp_r++; // unique hash -> map_id
    }
    assert(Sp_R.size() == count_sp_r);
    // clear the read hashes
    //Read_hashes.clear();

    // get a map id -> hash
    std::map<int32_t, uint64_t> Sp_R_hash;
    for (const auto& hash_pair : Sp_R) {
        Sp_R_hash[hash_pair.second] = hash_pair.first; // map_id -> hash
    }

    // print Indexed reads with spectrum size: Sp_R.size()
    fprintf(stderr, "[M::%s::%.3f*%.2f] Indexed reads with spectrum size: %d\n", __func__, realtime() - mg_realtime0, cputime() / (realtime() - mg_realtime0), Sp_R.size());

    Anchor_hits.assign(Sp_R.size(), std::vector<std::vector<std::vector<int32_t>>>(num_walks));
    std::vector<int32_t> count_anchors_hap(num_walks, 0);
    // compute the anchors
    for (int32_t h = 0; h < num_walks; h++)
    {
        // if (paths[h].size() == 0) continue;
        std::vector<std::vector<std::vector<int32_t>>> loc_match = compute_anchors(kmer_index[h], Sp_R); // parallel execution
        for (int32_t r = 0; r < Sp_R.size(); r++)
        {
            for (auto anchor: loc_match[r])
            {
                Anchor_hits[r][h].push_back(anchor);
            }
            if (loc_match[r].size() > 0) count_anchors_hap[h]++; // count the number of unique anchors
        } 
    }
    Sp_R.clear();

    // find number of kmers
    int32_t num_kmers_tot = 0;
    for (int32_t h = 0; h < num_walks; h++)
    {
        int32_t loc_count = 0;
        for (int32_t r = 0; r < count_sp_r; r++)
        {
            loc_count += Anchor_hits[r][h].size();
        }
        num_kmers_tot += loc_count;
    }

    std::vector<std::vector<std::vector<std::vector<int32_t>>>> Anchor_hits_1(count_sp_r, std::vector<std::vector<std::vector<int32_t>>>(num_walks));

    std::vector<int64_t> filtered_kmers_vec(num_threads, 0);
    #pragma omp parallel for num_threads(num_threads)
    for (int32_t r = 0; r < count_sp_r; r++) {
        std::map<std::string, std::pair<int32_t, std::vector<std::pair<int32_t, std::vector<int32_t>>>>> Anchor_hits_map;

        for (int32_t h = 0; h < num_walks; h++) {
            if (paths[h].size() == 0) continue;
            for (int32_t k = 0; k < Anchor_hits[r][h].size(); k++) {
                std::string anchor_str;
                for (auto v : Anchor_hits[r][h][k]) {
                    anchor_str += std::to_string(v) + "_";
                }

                if (Anchor_hits_map.find(anchor_str) == Anchor_hits_map.end()) { // does not exist
                    Anchor_hits_map[anchor_str].first = 1;
                    Anchor_hits_map[anchor_str].second.push_back(std::make_pair(h, Anchor_hits[r][h][k]));
                } else {
                    Anchor_hits_map[anchor_str].first++;
                    Anchor_hits_map[anchor_str].second.push_back(std::make_pair(h, Anchor_hits[r][h][k]));
                }
            }
        }

        bool all_haps = false;
        for (const auto &anchor : Anchor_hits_map)
        {
            if (anchor.second.first >= threshold * num_walks) {
                all_haps = true;
                break;
            }
        }

        if (!all_haps) { // only keep the minimizer if any of its anchor is not shared by all the haplotypes
            for (const auto &anchor : Anchor_hits_map) {
                for (const auto &hap_anchor : anchor.second.second) {
                    Anchor_hits_1[r][hap_anchor.first].push_back(hap_anchor.second);
                }
            }
        }else {
            filtered_kmers_vec[omp_get_thread_num()] += 1;
        }
    }


    // Replace Anchor_hits with Anchor_hits_1
    Anchor_hits = std::move(Anchor_hits_1);
    Anchor_hits_1.clear();

    // reorder occurrences of hits within each haplotype
    for (size_t r = 0; r < Anchor_hits.size(); ++r) {
        for (size_t h = 0; h < Anchor_hits[r].size(); ++h) {
            std::sort(
                Anchor_hits[r][h].begin(),
                Anchor_hits[r][h].end(),
                [](const std::vector<int32_t>& a,
                const std::vector<int32_t>& b)
                {
                    /* put empty vectors last */
                    if (a.empty()) return false;   // b precedes a
                    if (b.empty()) return true;    // a precedes b

                    /* primary key: left-most vertex */
                    if (a[0] != b[0])
                        return a[0] < b[0];

                    /* secondary key: earlier *ending* vertex
                    (i.e. shorter / earlier-ending path first) */
                    return a.back() < b.back();
                }
            );
        }
    }


    // kmer classifier *****************************************************


    int64_t filtered_kmers = 0;
    for (int32_t i = 0; i < num_threads; i++) filtered_kmers += filtered_kmers_vec[i];
    int64_t retained_kmers = count_sp_r - filtered_kmers;
    filtered_kmers_vec.clear();
    // find number of kmers
    std::cerr << "Number of Anchors" << std::endl;
    for (int32_t h = 0; h < num_walks; h++)
    {
        int32_t loc_count = 0;
        for (int32_t r = 0; r < count_sp_r; r++)
        {
            loc_count += Anchor_hits[r][h].size();
        }
        num_kmers += loc_count;
        std::string hap_name = hap_id2name[h];
        fprintf(stderr, "%s : %d\n", hap_name.c_str(), loc_count);
    }

    // pritn num_kmers/num_kmers_tot * 100 % are part of the haplotype
    //fprintf(stderr, "[M::%s::%.3f*%.2f] Filtered/Retained Minimizers: %.2f/%.2f%\n",
        __func__, 
        realtime() - mg_realtime0, 
        cputime() / (realtime() - mg_realtime0),
        (float)filtered_kmers/(float)count_sp_r * 100,
        (float)retained_kmers/(float)count_sp_r * 100;

    std::vector<std::vector<int32_t>> in_nodes(n_vtx);
    for (int32_t i = 0; i < n_vtx; i++)
    {
        for (auto v : adj_list[i])
        {
            in_nodes[v].push_back(i);
        }
    }

    // Compute Heterozygous and Homozygous Kmers
    
    // Iterate over all the reads and do +=1 if kmer is seen more than once if kmer is not there in kmer_count then set it to 0
    // std::vector<std::set<uint64_t>> kmer_reads(num_reads);
    // std::map<uint64_t, int32_t> Sp_R;

    // Shared map
    std::map<uint64_t, int32_t> kmer_count;

    // Mutex for merging
    std::mutex map_mutex;

    #pragma omp parallel num_threads(num_threads)
    {
        std::map<uint64_t, int32_t> local_map;

        #pragma omp for nowait
        for (int32_t r = 0; r < num_reads; r++) {
            for (auto hash : Read_hashes[r]) {
                local_map[hash] += 1;
            }
        }

        // Merge local_map into global map (serial part)
        std::lock_guard<std::mutex> lock(map_mutex);
        for (const auto& kv : local_map) {
            kmer_count[kv.first] += kv.second;
        }
    }

    // // print first 10 elements from kmer_count
    // for (auto it = kmer_count.begin(); it != std::next(kmer_count.begin(), 1000) && it != kmer_count.end(); ++it) {
    //     fprintf(stderr, "[%lu, %d]\n", it->first, it->second);
    // }

    // exit(0);


    Read_hashes.clear();

    // build {multiplicity, freq} table of kmers
    std::map<int32_t, int32_t> kmer_freq;
    std::map<int32_t, std::vector<uint64_t>> rev_kmer_freq;
    for (const auto& kv : kmer_count) {
        kmer_freq[kv.second] += 1;
        rev_kmer_freq[kv.second].push_back(kv.first);
    }

    std::vector<HistBin> Hist_kmer; // {multiplicity, freq}
    for (const auto& kv : kmer_freq) {
        Hist_kmer.push_back({(int)kv.first, (double)kv.second});
    }

    // // print first 40 elements from Hist_kmer
    // for (size_t i = 0; i < std::min((size_t)40, Hist_kmer.size()); i++) {
    //     fprintf(stderr, "[%d, %.2f]\n", Hist_kmer[i].multiplicity, Hist_kmer[i].freq);
    // }

    // find max_multiplicity and max_freq
    int max_multiplicity = 0;
    double max_freq = 0.0;
    for (const auto& bin : Hist_kmer) {
        if (bin.multiplicity > max_multiplicity) {
            max_multiplicity = bin.multiplicity;
        }
        if (bin.freq > max_freq) {
            max_freq = bin.freq;
        }
    }

    // print max_multiplicity and max_freq
    //fprintf(stderr, "Max Multiplicity: %d, Max Frequency: %.2f\n", max_multiplicity, max_freq);

    KGFitOptions opt;
    opt.max_copy = 10;
    opt.max_x_use = max_multiplicity;
    opt.u_hi = max_multiplicity;
    opt.fit_error = true;
    opt.fit_varw  = true;

    std::cout << "Classifying kmers..." << std::endl;
    auto res = KGFitterBO::fit(Hist_kmer, opt);
    KmerGenieDiploidLike model(res.P);

    // std::cerr << "best NLL="<<res.nll<<"\n"
    //     << "u_v="<<res.P.u_v<<" (hom SD) sd_v="<<res.P.sd_v
    //     << "\n(het VAR) var="<<res.P.var_w
    //     << "\n(proportion het)p_d="<<res.P.p_d
    //     << "\nzp_copy="<<res.P.zp_copy<<" zp_copy_het="<<res.P.zp_copy_het
    //     << "\n(zipf exponent) s="<<res.P.err_shape
    //     << "\nmax_copy="<<res.P.max_copy<<"\n";

    fprintf(stderr,
    "[M::%s] Fitted model: best NLL=%.2f, u_v=%.2f (hom mean), sd_v=%.2f (hom SD), "
    "var_w=%.2f, p_d=%.2f, zp_copy=%.2f, zp_copy_het=%.2f, err_shape=%.2f, max_copy=%d\n",
    __func__,
    res.nll,
    res.P.u_v,
    res.P.sd_v,
    res.P.var_w,
    res.P.p_d,
    res.P.zp_copy,        // %.2f (was %d)
    res.P.zp_copy_het,    // %.2f (was %d)
    res.P.err_shape,
    res.P.max_copy        // %d is correct for the int
    );
        
    // 2) build two lookup sets
    std::unordered_set<uint64_t> S_homo, S_hetero, S_ambiguous;
    S_homo.reserve(kmer_count.size());
    S_hetero.reserve(kmer_count.size());
    S_ambiguous.reserve(kmer_count.size());
    // int32_t total_count = 0;
    
    KmerGenieDiploidLike M(res.P);
    // classify a few multiplicities for output
    // for (int m = 1; m < max_multiplicity; m++) {
    //     auto r = M.classify(m);
    //     const char* L = (r.label==KGPosterior::HET?"HET":(r.label==KGPosterior::HOM?"HOM":"AMB"));
    //     if(m <= 20) std::cout << "m="<<m<<"  post[het,hom]=("<<r.p_het<<","<<r.p_hom<<")  "<<L<<"\n";
    // }

    // before the loop, use plain ints
    int64_t count_homo  = 0;
    int64_t count_hetero= 0;

    homo_bv.assign(count_sp_r,0);
    Anchor_hits_homo.assign(count_sp_r, std::vector<std::vector<std::vector<int32_t>>>(num_walks));
    Anchor_hits_hetero.assign(count_sp_r, std::vector<std::vector<std::vector<int32_t>>>(num_walks));
    Anchor_hits_ambiguous.assign(count_sp_r, std::vector<std::vector<std::vector<int32_t>>>(num_walks));

    //#pragma omp parallel for reduction(+:total_kmers,count_homo,count_hetero) num_threads(num_threads)
    for (int32_t kmer_idx = 0; kmer_idx < Anchor_hits.size(); ++kmer_idx) {
        auto it = Sp_R_hash.find(kmer_idx);
        if (it == Sp_R_hash.end()) continue; // shouldn’t happen

        uint64_t hash = it->second;
        int multiplicity = 0;
        auto kc = kmer_count.find(hash);
        if (kc != kmer_count.end()) multiplicity = kc->second;

        if (multiplicity == 0) continue; // unseen or filtered; skip

        auto classification = M.classify(multiplicity);
        bool is_homo = (classification.label == KGPosterior::HOM);
        bool is_hetero = (classification.label == KGPosterior::HET);

        if (is_homo + is_hetero == 0) {
            fprintf(stderr, "[M::%s] Error: k-mer %lu==  is not classified\n",
                __func__, kmer_idx);
            exit(1);
        }

        // check only one of is_homo, is_hetero, is_ambiguous is true
        if (is_homo + is_hetero != 1) { // Sanity check
            fprintf(stderr, "[M::%s] Error: k-mer %lu== classified as multiple types\n",
                __func__, kmer_idx);
            exit(1);
        }

        if(is_homo){
            homo_bv[kmer_idx] = 1;
            for(int h = 0; h < num_walks; h++){
                auto &out = Anchor_hits_homo[kmer_idx][h];
                out.insert(out.end(), Anchor_hits[kmer_idx][h].begin(), Anchor_hits[kmer_idx][h].end());
            }
            count_homo++;
        }else{
            for(int h = 0; h < num_walks; h++){
                auto &out = Anchor_hits_hetero[kmer_idx][h];
                out.insert(out.end(), Anchor_hits[kmer_idx][h].begin(), Anchor_hits[kmer_idx][h].end());
            }
            count_hetero++;
        }

    }

    float frac_h  = float(count_homo)   / std::max<int64_t>(1, count_homo+count_hetero);
    float frac_he = float(count_hetero) / std::max<int64_t>(1, count_homo+count_hetero);
    fprintf(stderr,
        "[M::%s] Phasing done. Homozygous: %.2f%%, Heterozygous: %.2f%%, Total kmers: %lld\n",
        __func__, frac_h*100, frac_he*100, count_homo+count_hetero);

}

