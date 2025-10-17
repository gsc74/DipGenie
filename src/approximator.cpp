#include "approximator.h"
#include "graph.hpp"
#include <cstdint>
#include <functional>
#include <numeric>
#include <random>
#include <cmath>
#include <queue>
#include <unordered_set>
//#include "edlib.h"
#include "misc.h"
#include "Classifier.hpp"
#include "Fitter.hpp"
#include <iomanip>
#include <chrono>
#include <sstream>

// Constructor
Approximator::Approximator(gfa_t *g) {
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


void::Approximator::read_gfa() 
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

    // Optionally expose the per-vertex MSA column if you want it directly:
    // (You can keep this around as a member if helpful.)
    // this->msa_column = std::move(dense_pos);

    // --- Optional: build "skip edges" per walk (jumps over >= 1 columns) ---
    // For each consecutive pair on a walk, if dense_pos jumps by > 1, record a skip.
    // This can help you explicitly model gaps in your haplotype paths.
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

    // If you want to replace/augment adj_list with edges in MSA order:
    //  - You can keep your existing graph edges as-is.
    //  - Or you can build a new adjacency that uses "columns" (dense_pos) as the sort key
    //    and includes skip_edges for haplotype gap jumps when cv > cu + 1.
    //
    // Example: stable-sort each adjacency by the MSA order of targets.
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


void remove_duplicates(std::vector<int>& vec) {
    std::unordered_set<int> seen;
    std::vector<int> result;
    for (int num : vec) {
        if (seen.find(num) == seen.end()) {
            seen.insert(num);
            result.push_back(num);
        }
    }
    vec = std::move(result);
}



// Read the reads
void Approximator::read_ip_reads(std::vector<std::pair<std::string, std::string>> &ip_reads, std::string ip_reads_file)
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


std::vector<std::pair<uint64_t, Anchor>> Approximator::index_kmers(int32_t hap) {

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


std::set<uint64_t> Approximator::compute_hashes(std::string &read_seq) {
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


std::vector<std::vector<std::vector<int32_t>>> Approximator::compute_anchors(std::vector<std::pair<uint64_t, Anchor>> &minimizers, std::map<uint64_t, int32_t> &read_hashes)
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



// treat all colors as distinct
std::vector<int> Approximator::dp_approximation_solver(Graph g, int R) {

    int n = g.adj_list.size();

    const std::size_t N = static_cast<std::size_t>(n) * (R + 1);

    std::vector<int> dp(N, 0);
    std::vector<int> back_vtx(N, -1);
    std::vector<int> back_r  (N, -1);

    auto idx = [&](int v, int r) -> std::size_t { return std::size_t(v)*(R+1) + r; };
    for(int u = 0; u < n; u++){
        for(int r = 0; r <= R; r++){
            //std::cout << "u: " << u << " r: " << r << " dp[u][r]: " << dp[u][r] << std::endl;
            for(auto& [v, w_uv] : g.adj_list[u]){
                
                if(r+w_uv <= R && dp[idx(u,r)] + g.color[v].size() > dp[idx(v,r+w_uv)]){
                    dp[idx(v, r+w_uv)] = dp[idx(u,r)] + g.color[v].size();
                    back_vtx[idx(v, r+w_uv)] = u;
                    back_r[idx(v, r+w_uv)] = r;
                }
            }
        }
    }


    std::vector<int> colors_by_r;
    //std::vector<int> edit_dist_by_r;
    //for(int r = 0; r <= R; r++){
    std::vector<std::map<int,int>> occ_count_by_r;
    for(int r = 0; r <= R; r++){
        //std::vector<int>path_temp;
        std::unordered_set<int> true_colours;
        std::map<int,int> occ_count;
        // backtrack to obtain path
        int cur_vtx = n-1;
        int cur_r = r;
        while(cur_vtx != -1){
            //path_temp.push_back(cur_vtx); // comment
            for(auto c : g.color[cur_vtx]){
                true_colours.insert(c);
                if(occ_count.find(c) == occ_count.end()){
                    occ_count[c] = 1;
                }else{
                    occ_count[c] = occ_count[c]+1;
                }
            }
            int temp_vtx = cur_vtx;
            
            //cur_vtx = backptr_vtx[cur_vtx][cur_r];
            //cur_r = backptr_r[temp_vtx][cur_r];

            cur_vtx = back_vtx[idx(cur_vtx, cur_r)];
            cur_r = back_r[idx(temp_vtx, cur_r)];
        }
        colors_by_r.push_back(static_cast<int>(true_colours.size()));
        occ_count_by_r.push_back(occ_count);

    }
    //compute approximation certificate
    for (size_t i = 0; i + 1 < occ_count_by_r.size(); ++i) {
        float avg = 0;
        const auto& occ_count = occ_count_by_r[i];
        for(const auto& c: occ_count){
            avg += c.second;
            //std::cout << c.first << ".second = " << c.second << std::endl;
        }
        avg = avg/occ_count.size();
        std::cout << "Approximation ratio certificate: " << avg << std::endl;
    }

    // compute best r
    int best_r = 0;
    double max_delta = 0;
    for (size_t i = 0; i + 1 < colors_by_r.size(); ++i) {
        std::cout << "r: " << i << " true score: " << colors_by_r[i] << std::endl;
        int delta = colors_by_r[i + 1] - colors_by_r[i];
        if (std::abs(delta) > max_delta) max_delta = std::abs(delta);
    }

    for (size_t r = 0; r + 1 < colors_by_r.size(); ++r) {
        int delta = colors_by_r[r + 1] - colors_by_r[r];
        double angle_rad = std::atan(static_cast<double>(delta) / max_delta);
        double angle_deg = angle_rad * 180.0 / M_PI;
        std::cout << "r: " << r << " -> " << r + 1 << ", Δcolors: " << delta
                << ", angle: " << angle_deg << "°" //<<", ed: " << edit_dist_by_r.at(i)
                << std::endl;
        if(angle_deg < 10){
            best_r = r;
            break;
        }

    }

    std::cerr << "Recombination count: " << best_r << std::endl;

    // recover corresponding path
    std::vector<int>path;

    int cur_vtx = n-1;
    int cur_r = best_r;
    while(cur_vtx != -1){
        path.push_back(cur_vtx);
        int temp_vtx = cur_vtx;

        cur_vtx = back_vtx[idx(cur_vtx, cur_r)];
        cur_r = back_r[idx(temp_vtx, cur_r)];
    }

    std::reverse(path.begin(), path.end());
    std::vector<int> dp_path_in_original_graph;
    for(auto u : path){
        for(auto u_original : g.original_vertex[u]){
            dp_path_in_original_graph.push_back(u_original);
        }
    }

    auto uniq_consecutive = [](std::vector<int>& v){
        v.erase(std::unique(v.begin(), v.end()), v.end());
    };
    remove_duplicates(dp_path_in_original_graph);

    return dp_path_in_original_graph;
}


template <class T>
std::size_t intersection_size(const std::set<T>& a, const std::set<T>& b) {
    std::size_t c = 0;
    auto i = a.begin(), j = b.begin();
    while (i != a.end() && j != b.end()) {
        if (*i < *j) ++i;
        else if (*j < *i) ++j;
        else { ++c; ++i; ++j; }     // match
    }
    return c;                       // |A ∩ B|
}

template <class T>
std::size_t symdiff_size(const std::set<T>& a, const std::set<T>& b) {
    std::size_t inter = intersection_size(a, b);
    return a.size() + b.size() - 2 * inter;   // |A Δ B| = |A|+|B|−2|A∩B|
}


template <class T>
inline void dedup_inplace(std::vector<T>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

template <class T>
inline std::size_t intersection_size_sorted(const std::vector<T>& a,
                                            const std::vector<T>& b) {
    std::size_t i = 0, j = 0, c = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++c; ++i; ++j; }
    }
    return c;
}

template <class T>
inline std::size_t symdiff_size_sorted(const std::vector<T>& a,
                                       const std::vector<T>& b) {
    std::size_t i = 0, j = 0, c = 0;
    while (i < a.size() || j < b.size()) {
        if (j == b.size() || (i < a.size() && a[i] < b[j])) { ++c; ++i; }
        else if (i == a.size() || b[j] < a[i]) { ++c; ++j; }
        else { ++i; ++j; } // equal → not in symmetric diff
    }
    return c;
}

template <class T>
inline std::size_t union_size_sorted(const std::vector<T>& a,
                                     const std::vector<T>& b) {
    std::size_t i = 0, j = 0, c = 0;
    while (i < a.size() || j < b.size()) {
        if (j == b.size() || (i < a.size() && a[i] < b[j])) { ++c; ++i; }
        else if (i == a.size() || b[j] < a[i]) { ++c; ++j; }
        else { ++c; ++i; ++j; } // equal → count once
    }
    return c;
}

// Count |(A ∪ B) ∩ (C ∪ D)| without building unions
inline int inter_size_union2x2(const std::vector<int>& A, const std::vector<int>& B,
                               const std::vector<int>& C, const std::vector<int>& D) {
    size_t i=0,j=0,k=0,m=0; int cnt=0;
    while (i<A.size() || j<B.size() || k<C.size() || m<D.size()) {
        int x = std::numeric_limits<int>::max();
        if (i<A.size()) x = std::min(x, A[i]);
        if (j<B.size()) x = std::min(x, B[j]);
        if (k<C.size()) x = std::min(x, C[k]);
        if (m<D.size()) x = std::min(x, D[m]);

        bool inLeft  = false, inRight = false;
        while (i<A.size() && A[i]==x) { inLeft  = true; ++i; }
        while (j<B.size() && B[j]==x) { inLeft  = true; ++j; }
        while (k<C.size() && C[k]==x) { inRight = true; ++k; }
        while (m<D.size() && D[m]==x) { inRight = true; ++m; }

        if (inLeft && inRight) ++cnt;
    }
    return cnt;
}

// Count |(E ∪ F) △ (G ∪ H)| without building unions
inline int symdiff_size_union2x2(const std::vector<int>& E, const std::vector<int>& F,
                                 const std::vector<int>& G, const std::vector<int>& H) {
    size_t i=0,j=0,k=0,m=0; int cnt=0;
    while (i<E.size() || j<F.size() || k<G.size() || m<H.size()) {
        int x = std::numeric_limits<int>::max();
        if (i<E.size()) x = std::min(x, E[i]);
        if (j<F.size()) x = std::min(x, F[j]);
        if (k<G.size()) x = std::min(x, G[k]);
        if (m<H.size()) x = std::min(x, H[m]);

        bool inLeft  = false, inRight = false;
        while (i<E.size() && E[i]==x) { inLeft  = true; ++i; }
        while (j<F.size() && F[j]==x) { inLeft  = true; ++j; }
        while (k<G.size() && G[k]==x) { inRight = true; ++k; }
        while (m<H.size() && H[m]==x) { inRight = true; ++m; }

        if (inLeft ^ inRight) ++cnt;   // XOR → symmetric diff
    }
    return cnt;
}

// Count |(A ∪ B ∪ C ∪ D)| without building unions
inline int union_size_union2x2(const std::vector<int>& A, const std::vector<int>& B,
                               const std::vector<int>& C, const std::vector<int>& D) {
    size_t i=0,j=0,k=0,m=0; 
    int cnt=0;
    while (i<A.size() || j<B.size() || k<C.size() || m<D.size()) {
        int x = std::numeric_limits<int>::max();
        if (i<A.size()) x = std::min(x, A[i]);
        if (j<B.size()) x = std::min(x, B[j]);
        if (k<C.size()) x = std::min(x, C[k]);
        if (m<D.size()) x = std::min(x, D[m]);

        bool present = false;
        while (i<A.size() && A[i]==x) { present = true; ++i; }
        while (j<B.size() && B[j]==x) { present = true; ++j; }
        while (k<C.size() && C[k]==x) { present = true; ++k; }
        while (m<D.size() && D[m]==x) { present = true; ++m; }

        if (present) ++cnt;
    }
    return cnt;
}


#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>

std::vector<int> gather_colors(Graph g, int start, int end) {
    int last_level = g.level.at(end);
    std::vector<int> colours;

    if (g.level.at(start) > last_level) return colours;

    std::unordered_set<int> visited;
    std::queue<int> q;

    // color -> seen at this level?
    std::unordered_map<int, std::unordered_set<int>> seen_at_level;

    visited.insert(start);
    q.push(start);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        int lu = g.level[u];

        // collect colors at u (once per level)
        auto& seen_here = seen_at_level[lu];
        const auto& cols = g.color[u];
        for (int c : cols) {
            if (seen_here.insert(c).second) {
                colours.push_back(c);
            }
        }

        // expand along weight-0 edges to allowed levels
        for (const auto& [v, w_uv] : g.adj_list[u]) {
            if (w_uv != 0) continue;
            if (g.level[v] <= last_level && visited.insert(v).second) {
                q.push(v);
            }
        }
    }

    return colours;
}


static std::string format_hms(double sec) {
    int s = (int)std::round(sec);
    int h = s / 3600; s %= 3600;
    int m = s / 60;   s %= 60;
    std::ostringstream os;
    if (h) os << h << "h";
    if (h || m) os << m << "m";
    os << s << "s";
    return os.str();
}

static void progress_bar(std::size_t current, std::size_t total,
                         std::chrono::steady_clock::time_point start,
                         std::size_t width = 40)
{
    using namespace std::chrono;
    double frac   = total ? (double)current / (double)total : 1.0;
    std::size_t n = (std::size_t)std::floor(frac * width);

    double elapsed = duration<double>(steady_clock::now() - start).count();
    double rate    = elapsed > 0 ? current / elapsed : 0.0;                // it/s
    double eta     = rate > 0 ? (total > current ? (total - current)/rate : 0.0) : 0.0;

    std::ostringstream line;
    line << '\r' << '[';
    for (std::size_t i = 0; i < width; ++i)
        line << (i < n ? '=' : (i == n ? '>' : ' '));
    line << "] " << std::setw(3) << int(frac * 100) << "%  "
         << current << '/' << total
         << "  | " << std::fixed << std::setprecision(1) << rate << " it/s"
         << "  | ETA " << format_hms(eta)
         << "         ";

    std::cout << line.str() << std::flush;
    if (current == total) std::cout << '\n';
}





std::vector<std::tuple<int, int, std::string, std::string>> Approximator::diploid_dp_approximation_solver(Graph g, int R, std::vector<bool> color_homo_bv) {
    
    
    // Build vertex -> position-in-its-level map once
    const int L = (int)g.vertices_in_level.size();
    std::vector<int> pos_in_level(g.adj_list.size(), -1);
    for (int l = 0; l < L; ++l) {
        const auto& Lv = g.vertices_in_level[l];
        for (int i = 0; i < (int)Lv.size(); ++i) pos_in_level[ Lv[i] ] = i;
    }

    // Backpointers: allocate per level (only from l to l+1, so store at level l+1)
    //std::vector<std::vector<std::pair<std::uint16_t,std::uint16_t>>> back_vtx(L);
    //std::vector<std::vector<std::uint8_t>>  back_r(L);
    //const std::uint16_t NO_VTX = std::numeric_limits<std::uint16_t>::max();
    //const std::uint8_t  NO_R   = std::numeric_limits<std::uint8_t>::max();

    struct dp_entry {
        int value;
        int s_het;  // for approximation cert later

        std::vector<std::pair<int, int>> weighted_p1_edges; 
        std::vector<std::pair<int, int>> weighted_p2_edges; 

        dp_entry(int v = 0, int s = 0) : value(v), s_het(0) {}
    };

    // Rolling DP buffers
    const int32_t NEG_INF = std::numeric_limits<int32_t>::min() / 4;
    std::vector<dp_entry> dp_cur, dp_next;

    // helper to index flattened (R+1)×k×k
    auto at3 = [](std::vector<dp_entry>& buf, int k, int i, int j, int r) -> dp_entry& {
        return buf[ ((std::size_t)r * k + i) * k + j ];
    };

    // Build per-vertex sorted, deduped lists for homo/hetero colors
    std::cout << "Creating hetro/hom-zygous colors per vertex lists" << std::endl;
    std::vector<std::vector<int>> homo_sorted(g.adj_list.size());
    std::vector<std::vector<int>> hetero_sorted(g.adj_list.size());

    auto sort_dedup = [](std::vector<int>& v){
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    };

    for (int v = 0; v < (int)g.adj_list.size(); ++v) {
        auto& H = homo_sorted[v];
        auto& T = hetero_sorted[v];
        H.reserve(g.color[v].size());
        T.reserve(g.color[v].size());
        for (int c : g.color[v]) {
            if(color_homo_bv[c]){ 
                H.push_back(c);
            }else{
                T.push_back(c);
            }
        }
        sort_dedup(H);
        sort_dedup(T);
    }

    // ---- DP profiling (lightweight) ---------------------------------
    #define ENABLE_DP_PROF 1
    using clock_ = std::chrono::steady_clock;

    double prof_progress = 0.0;
    double prof_init_cur = 0.0;
    double prof_alloc_next = 0.0;
    double prof_scored = 0.0;
    double prof_relax = 0.0;
    double prof_swap = 0.0;

    long long prof_score_pairs = 0;   // # of (i,j,u2,v2) pairs per level (score_deltas.size())
    long long prof_relax_checks = 0;  // # of relax candidates considered
    long long prof_relax_wins = 0;    // # of times we improved dst

    auto dp_wall_start = clock_::now();



    std::cout << "Running DP" << std::endl;
    auto t0 = std::chrono::steady_clock::now();


    for (int l = 0; l < L; ++l) {

        // progress bar (timed)
        {
            auto _t = clock_::now();
            if (l % 10 == 0 || l == L-1) {
                progress_bar((std::size_t)l + 1, (std::size_t)L, t0);
            }
            prof_progress += std::chrono::duration<double>(clock_::now() - _t).count();
        }

        const auto& Lnow  = g.vertices_in_level[l];
        const int   k     = (int)Lnow.size();

        // allocate/resize current level buffer if first iteration (timed)
        if (l == 0) {
            auto _t = clock_::now();
            const std::size_t sz0 = (std::size_t)(R + 1) * k * k;
            dp_cur.assign(sz0, dp_entry(0));
            prof_init_cur += std::chrono::duration<double>(clock_::now() - _t).count();
        }

        // last level has no outgoing edges; nothing to push to next
        if (l + 1 >= L) break;

        const auto& Lnext = g.vertices_in_level[l + 1];
        const int   k2    = (int)Lnext.size();

        // dp_next init (timed)
        const std::size_t szN = (std::size_t)(R + 1) * k2 * k2;
        {
            auto _t = clock_::now();
            //dp_next.assign(szN, dp_entry(NEG_INF));
            dp_next.resize(szN);  // grows without realloc if capacity is enough; shrinks but keeps capacity

            // reset in place: set value to NEG_INF and clear edge lists
            for (auto &e : dp_next) {
                e.value = NEG_INF;
                e.weighted_p1_edges.clear();
                e.weighted_p2_edges.clear();
            }
            prof_alloc_next += std::chrono::duration<double>(clock_::now() - _t).count();
        }

        // ---------------- precompute score_deltas (timed) ----------------
        std::vector<int> score_deltas;
        std::vector<int> s_hets;
        {
            auto _t = clock_::now();

            score_deltas.reserve((size_t)k * k * 4); // cheap guess to cut reallocs
            s_hets.reserve((size_t)k * k * 4); // cheap guess to cut reallocs
            for (int i = 0; i < k; ++i) {
                for (int j = 0; j < k; ++j) {
                    int u1 = Lnow[i];
                    int v1 = Lnow[j];

                    for (const auto& [u2, _wu] : g.adj_list[u1]) {
                        int iu2 = pos_in_level[u2];
                        if ((unsigned)iu2 >= (unsigned)k2) continue;
                        for (const auto& [v2, _wv] : g.adj_list[v1]) {
                            int jv2 = pos_in_level[v2];
                            if ((unsigned)jv2 >= (unsigned)k2) continue;

                            int inter = inter_size_union2x2( homo_sorted[u1], homo_sorted[v1],
                                                            homo_sorted[u2], homo_sorted[v2] );
                            int symd  = symdiff_size_union2x2( hetero_sorted[u1], hetero_sorted[v1],
                                                            hetero_sorted[u2], hetero_sorted[v2] );
                            s_hets.push_back(symd);
                            score_deltas.push_back(inter + symd);
                        }
                    }
                }
            }
            prof_score_pairs += (long long)score_deltas.size();
            prof_scored += std::chrono::duration<double>(clock_::now() - _t).count();
        }

        // ---------------- relaxation over r (timed) ---------------------
        {
            auto _t = clock_::now();

            for (int r = 0; r <= R; r ++) {
                std::size_t idx = 0;  // IMPORTANT: this must advance once per (u2,v2) pair
                for (int i = 0; i < k; ++i) {
                    for (int j = 0; j < k; ++j) {

                        dp_entry& src = at3(dp_cur, k, i, j, r);
                        if (src.value == NEG_INF) {
                            // still need to skip the same number of (u2,v2) entries:
                            // count how many (u2,v2) exist for this (i,j) to advance idx
                            int u1 = Lnow[i], v1 = Lnow[j];
                            for (const auto& [u2, _wu] : g.adj_list[u1]) {
                                int iu2 = pos_in_level[u2];
                                for (const auto& [v2, _wv] : g.adj_list[v1]) {
                                    int jv2 = pos_in_level[v2];
                                    ++idx; // consume one precomputed score
                                }
                            }
                            continue;
                        }

                        int u1 = Lnow[i];
                        int v1 = Lnow[j];

                        for (const auto& [u2, wu] : g.adj_list[u1]) {
                            int iu2 = pos_in_level[u2];

                            for (const auto& [v2, wv] : g.adj_list[v1]) {
                                int jv2 = pos_in_level[v2];

                                ++prof_relax_checks;

                                int r2 = r + wu + wv;
                                if (r2 > R) { ++idx; continue; }

                                auto& dst = at3(dp_next, k2, iu2, jv2, r2);
                                int cand = src.value + score_deltas[idx];
                                if (cand > dst.value) {
                                    ++prof_relax_wins;

                                    dst.value = cand;
                                    dst.s_het = s_hets[idx];
                                    dst.weighted_p1_edges = src.weighted_p1_edges;
                                    dst.weighted_p2_edges = src.weighted_p2_edges;

                                    if (wu > 0) dst.weighted_p1_edges.emplace_back(u1, u2);
                                    if (wv > 0) dst.weighted_p2_edges.emplace_back(v1, v2);

                                    if (l + 1 == L - 1) {
                                        dst.weighted_p1_edges.emplace_back(u1, u2);
                                        dst.weighted_p2_edges.emplace_back(v1, v2);
                                    }
                                }
                                ++idx; // advance exactly once per (u2,v2)
                            } // v2
                        } // u2
                    } // j
                } // i
            } // r

            prof_relax += std::chrono::duration<double>(clock_::now() - _t).count();
        }

        // ---------------- roll buffers (timed) --------------------------
        {
            auto _t = clock_::now();
            //std::vector<dp_entry>().swap(dp_cur);
            dp_cur.swap(dp_next);
            //std::vector<dp_entry>().swap(dp_next);
            prof_swap += std::chrono::duration<double>(clock_::now() - _t).count();
        }
    } // l





    // ---- print profile summary ---------------------------------------
    double prof_total = std::chrono::duration<double>(clock_::now() - dp_wall_start).count();

    auto pct = [&](double t){ return (prof_total > 0) ? (100.0 * t / prof_total) : 0.0; };

    std::cout << "\n=== DP profile (seconds) ===\n"
            << " progress UI:   " << prof_progress   << "s (" << pct(prof_progress)   << "%)\n"
            << " init dp_cur:   " << prof_init_cur   << "s (" << pct(prof_init_cur)   << "%)\n"
            << " alloc dp_next: " << prof_alloc_next << "s (" << pct(prof_alloc_next) << "%)\n"
            << " score_deltas:  " << prof_scored     << "s (" << pct(prof_scored)     << "%)\n"
            << " relax/update:  " << prof_relax      << "s (" << pct(prof_relax)      << "%)\n"
            << " roll buffers:  " << prof_swap       << "s (" << pct(prof_swap)       << "%)\n"
            << " total (wall):  " << prof_total      << "s (100%)\n"
            << " pairs precomp: " << prof_score_pairs << "\n"
            << " relax checks:  " << prof_relax_checks << "\n"
            << " relax wins:    " << prof_relax_wins << "\n";



    std::cout << "\nDone with diploid DP, reporting solutions..." << std::endl;

    // indexer for back tables at a given level l (r-major, then i, then j)
    auto idx3_level = [&](int l, int i, int j, int r) -> std::size_t {
        int k = (int)g.vertices_in_level[l].size();
    #ifndef NDEBUG
        if ((unsigned)i >= (unsigned)k || (unsigned)j >= (unsigned)k) throw std::out_of_range("i/j");
    #endif
        return ((std::size_t)r * k + i) * k + j;
    };

    int k_sink = (int)g.vertices_in_level.back().size(); // should be 1

    // // for debug
    // Outputting weighted edges, use these and comment out DP above for evaluating different scoring
    for( int r = 0; r <= R; r++ ) {
        int start = g.vertices_in_level.at(0).at(0);
        std::cout << "r = " << r << std::endl;
        std::cout << "P1 weighted edges: " << std::endl;
        auto& sink_dp = at3(dp_cur, k_sink, 0, 0, r); // corresponds to single sink vertex
        for( auto edge : sink_dp.weighted_p1_edges ) {
            std::cout << "(" << edge.first << ", " << edge.second << ")" << std::endl;
        }

        std::cout << "P2 weighted edges: " << std::endl;
        for( auto edge : sink_dp.weighted_p2_edges ) {
            std::cout << "(" << edge.first << ", " << edge.second << ")" << std::endl;
        }
    }


    // Find the nearest vertex (by number of 0-weight hops) from src whose hap == target_hap.
    // Returns -1 if none reachable via only 0-weight edges.
    auto find_next_zero_hap = [&](int src, int target_hap) -> int {
        // If zero hops are allowed and src already matches, return it.
        if (g.haplotype.at(src) == target_hap && g.original_vertex.at(src).size() > 0) return src;

        std::queue<int> q;
        std::unordered_set<int> visited;
        q.push(src);
        visited.insert(src);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (const auto& [v, w] : g.adj_list[u]) {
                if (w != 0) continue;                 // only follow weight-0 edges
                if (!visited.insert(v).second) continue;

                if (g.haplotype.at(v) == target_hap && g.original_vertex.at(v).size() > 0){  // first hit = closest by BFS
                    return v;
                }

                q.push(v);
            }
        }
        return -1; // not found
    };
    std::vector<std::tuple<int, int, std::string, std::string>> solutions;

    for( int r = 0; r <= R; r++ ) {

        std::cout << "Computing haplotype paths for r = " << r << std::endl;
        auto& sink_dp = at3(dp_cur, k_sink, 0, 0, r); // corresponds to single sink vertex
        int r1 = sink_dp.weighted_p1_edges.size()-1;
        int r2 = sink_dp.weighted_p2_edges.size()-1;
        std::cerr << "Recombinations in P1: " << r1 << std::endl;
        std::cerr << "Recombinations in P2: " << r2 << std::endl;

        int start_exp = g.vertices_in_level.at(0).at(0);
        int h;

        std::string hap_1 = "";
        for(int i = 0; i < sink_dp.weighted_p1_edges.size(); i++){
            auto& edge = sink_dp.weighted_p1_edges.at(i);

            //std::cout << "(" << edge.first << ", " << edge.second << ")" << std::endl;
            if(g.original_vertex[edge.first].size() != 1){
                std::cout << "P1: Vertex " << edge.first << " in map back has " << g.original_vertex[edge.first].size() << " original vertices" << std::endl;
                exit(1);
            }
            int end_exp = edge.first;

            //std::cout << "P1 (string): getting h = " << h << std::endl;
            h = g.haplotype.at(end_exp);

            if(start_exp == g.vertices_in_level.at(0).at(0)){
                //std::cout << "P1 (string): obtaining first hap starting vertex" << std::endl;
                for(auto& v: g.vertices_in_level.at(1)){
                    if(g.haplotype.at(v) == h ){
                        start_exp = v;
                    }
                }
            }
            
            // go from vertex start to vertex end in paths[h], concatenating vertex labels
            //std::cout << "P1 (string): getting original start vertex: size: " << g.original_vertex.at(start_exp).size() << std::endl;
            int start_org = g.original_vertex.at(start_exp).at(0);
            //std::cout << "P1 (string): getting original end vertex" << std::endl;
            int end_org = g.original_vertex.at(end_exp).at(0);
            bool activated = false;
            for(int i = 0; i < paths[h].size(); i++){
                if(paths[h][i] == start_org){
                    activated = true;
                }
                if(activated){
                    hap_1 += node_seq[paths[h][i]];
                }
                if(paths[h][i] == end_org){
                    activated = false;
                    break;
                }
            }

            if(g.level.at(edge.second) == L-1){
                //std::cout << "P1 (path recovery): done" << std::endl;
                break;
            }
            auto& next_edge = sink_dp.weighted_p1_edges.at(i+1);
            int next_hap = g.haplotype.at(next_edge.first);
            bool found_next = false;
            //std::cout << "P1 (string): Finding next starting point" << std::endl;

            int next_start = find_next_zero_hap(edge.second, next_hap);
            if (next_start != -1) {
                start_exp = next_start;
                found_next = true;
            } else {
                std::cout << "P1 (path recovery) Could not find next_hap=" << next_hap
                        << " from " << edge.second << " via 0-weight edges\n";
            }
            //std::cout << "P1 (string): Done traversing path for h = " << h << std::endl;
        }



        std::string hap_2 = "";
        start_exp = g.vertices_in_level.at(0).at(0);

        for(int i = 0; i < sink_dp.weighted_p2_edges.size(); i++){

            auto& edge = sink_dp.weighted_p2_edges.at(i);
            //std::cout << "(" << edge.first << ", " << edge.second << ")" << std::endl;
            if(g.original_vertex[edge.first].size() != 1){
                std::cout << "P2: Vertex "<< edge.first <<" in map back has " << g.original_vertex[edge.first].size() << " original vertices" << std::endl;
                exit(1);
            }
            int end_exp = edge.first;
            h = g.haplotype.at(end_exp);

            if(start_exp == g.vertices_in_level.at(0).at(0)){
                for(auto& v: g.vertices_in_level.at(1)){
                    if(g.haplotype.at(v) == h){
                        start_exp = v;
                    }
                }
            }
            // go from vertex start to vertex end in paths[h], concatenating vertex labels
            //std::cout << "P2 (string): getting original start vertex" << std::endl;
            int start_org = g.original_vertex.at(start_exp).at(0);
            //std::cout << "P2 (string): getting original end vertex" << std::endl;
            int end_org = g.original_vertex.at(end_exp).at(0);
            bool activated = false;
            for(int i = 0; i < paths[h].size(); i++){
                if(paths[h][i] == start_org){
                    activated = true;
                }
                if(activated){
                    hap_2 += node_seq[paths[h][i]];
                }
                if(paths[h][i] == end_org){
                    activated = false;
                    break;
                }
            }

            if(g.level.at(edge.second) == L-1){
                //std::cout << "P2 (string): done" << std::endl;
                break;
            }
            auto& next_edge = sink_dp.weighted_p2_edges.at(i+1);
            int next_hap = g.haplotype.at(next_edge.first);
            bool found_next = false;
            //std::cout << "Finding next starting point" << std::endl;

            int next_start = find_next_zero_hap(edge.second, next_hap);
            if (next_start != -1) {
                start_exp = next_start;
                found_next = true;
            } else {
                std::cout << "P2 (path recovery) Could not find next_hap=" << next_hap
                        << " from " << edge.second << " via 0-weight edges\n";
            }
            //std::cout << "P2: Done traversing path for h = " << h << std::endl;            

        }

        solutions.emplace_back(r1, r2, hap_1, hap_2);
    }




    // Delta score approach
    // computes score and approx cert.
    std::vector<int> score_by_r;

    for( int r = 0; r <= R; r++ ) {

        std::cout << "Computing score for r = " << r << std::endl;
        auto& sink_dp = at3(dp_cur, k_sink, 0, 0, r); // corresponds to single sink vertex

        int s_het = sink_dp.s_het;

        std::vector<int> path_1_homo_colors;
        std::vector<int> path_1_hetero_colors;

        int start = g.vertices_in_level.at(0).at(0);

        for(int i = 0; i < sink_dp.weighted_p1_edges.size(); i++){
            
            //std::cout << "Getting edge" << std::endl;
            auto& edge = sink_dp.weighted_p1_edges.at(i);

            int end = edge.first; // head of weighted edge
            
            //std::cout << "Getting haplotype" << std::endl;
            int h = g.haplotype.at(end);

            //std::cout << "h = " << h << std::endl;
            if(start == g.vertices_in_level.at(0).at(0)){
                for(auto& v: g.vertices_in_level.at(1)){
                    if(g.haplotype.at(v) == h && g.original_vertex.at(v).size() > 0){
                        start = v;
                    }
                }
            }

            //std::cout << "getting colors" << std::endl;
            std::vector<int> reachable_colors = gather_colors(g, start, end);
            
            for(auto c : reachable_colors){

                if(color_homo_bv[c]){
                    path_1_homo_colors.push_back(c);
                }else{
                    path_1_hetero_colors.push_back(c);
                }
            }
            //std::cout << "Checking level of edge.second" << std::endl;

            if(g.level.at(edge.second) == L-1){
                //std::cout << "P1 (score): done" << std::endl;
                break;
            }
            //std::cout << "Getting next edge" << std::endl;

            auto& next_edge = sink_dp.weighted_p1_edges.at(i+1);

            //std::cout << "Getting next edge haplotype" << std::endl;
            int next_hap = g.haplotype.at(next_edge.first);
            bool found_next = false;
            //std::cout << "Finding next starting point" << std::endl;

            int next_start = find_next_zero_hap(edge.second, next_hap);
            if (next_start != -1) {
                start = next_start;
                found_next = true;
            } else {
                std::cout << "P1 (score) Could not find next_hap=" << next_hap
                        << " from " << edge.second << " via 0-weight edges\n";
            }
        }


        std::vector<int> path_2_homo_colors;
        std::vector<int> path_2_hetero_colors;

        start = g.vertices_in_level.at(0).at(0);

        for(int i = 0; i < sink_dp.weighted_p2_edges.size(); i++){

            auto& edge = sink_dp.weighted_p2_edges.at(i);
            
            int end = edge.first; // head of weighted edge
            
            int h = g.haplotype.at(end);

            if(start == g.vertices_in_level.at(0).at(0)){
                for(auto& v: g.vertices_in_level.at(1)){
                    if(g.haplotype.at(v) == h && g.original_vertex.at(v).size() > 0){
                        start = v;
                    }
                }
            }

            std::vector<int> reachable_colors = gather_colors(g, start, end);
            
            for(auto c : reachable_colors){

                if(color_homo_bv[c]){
                    path_2_homo_colors.push_back(c);
                }else{
                    path_2_hetero_colors.push_back(c);
                }
            }
            if(g.level.at(edge.second) == L-1){
                //std::cout << "P2 (score): done" << std::endl;
                break;
            }
            auto& next_edge = sink_dp.weighted_p2_edges.at(i+1);
            int next_hap = g.haplotype.at(next_edge.first);
            bool found_next = false;
            //std::cout << "Finding next starting point" << std::endl;

            int next_start = find_next_zero_hap(edge.second, next_hap);
            if (next_start != -1) {
                start = next_start;
                found_next = true;
            } else {
                std::cout << "P2 (score) Could not find next_hap=" << next_hap
                        << " from " << edge.second << " via 0-weight edges\n";
            }
            //std::cout << "P2: Done gathering colors for h = " << h << std::endl;
        }


        // compute multiplicity values for hom
        int m_G_hom = 0;
        std::unordered_map<int, std::array<int,2>> cnt_hom; // {color -> {count_in_p1, count_in_p2}}

        for (int c : path_1_homo_colors) cnt_hom[c][0]++;  // occ_c(path_1)
        for (int c : path_2_homo_colors) cnt_hom[c][1]++;  // occ_c(path_2)

        for (const auto& [color, ab] : cnt_hom) {
            m_G_hom += std::max(ab[0], ab[1]);
        }


        // compute multiplicity values for het
        int m_G_het = 0;
        std::unordered_map<int, std::array<int,2>> cnt_het;
        for (int c : path_1_hetero_colors) cnt_het[c][0]++;  // occ_c(path_1)
        for (int c : path_2_hetero_colors) cnt_het[c][1]++;  // occ_c(path_2)

        for (const auto& [color, ab] : cnt_het) {
            m_G_het += ab[0] + ab[1];
        }

        dedup_inplace(path_1_homo_colors);
        dedup_inplace(path_1_hetero_colors);
        dedup_inplace(path_2_homo_colors);
        dedup_inplace(path_2_hetero_colors);
    
    
        int intersection_count = intersection_size_sorted(path_1_homo_colors, path_2_homo_colors);
        int symdiff_count = symdiff_size_sorted(path_1_hetero_colors, path_2_hetero_colors);

        // approximation cert.
        float m_G_hom_avg = m_G_hom / (float) intersection_count;
        float m_G_het_avg = m_G_het / (float) symdiff_count;
        float m_bar = std::max(m_G_hom_avg, m_G_het_avg);

        int loss_het = s_het - m_G_het;

        std::cout << "Approximation certificate: multiplicative factor: " << m_bar 
            << " additive factor: " << -loss_het/(float)m_G_het_avg << std::endl;
    
        int score = intersection_count + symdiff_count;
        std::cout << "r: " << r << " score: " << score << std::endl;
        score_by_r.push_back(score);

    }

    // compute best r
    int best_r = 0;
    double max_delta = 0;
    for (size_t i = 0; i + 1 < score_by_r.size(); ++i) {
        int delta = score_by_r[i + 1] - score_by_r[i];
        if (std::abs(delta) > max_delta) max_delta = std::abs(delta);
    }

    for (size_t r = 0; r + 1 < score_by_r.size(); ++r) {
        int delta = score_by_r[r + 1] - score_by_r[r];
        double angle_rad = std::atan(static_cast<double>(delta) / max_delta);
        double angle_deg = angle_rad * 180.0 / M_PI;
        std::cout << "r: " << r << " -> " << r + 1 << ", Δscore: " << delta
                << ", angle: " << angle_deg << "°" //<<", ed: " << edit_dist_by_r.at(i)
                << std::endl;
        if(angle_deg < 10){
            best_r = r;
            auto best_sol = solutions.at(r);
            solutions.clear();
            solutions.push_back(best_sol);
            break;
        }
    }

    if(solutions.size() > 1){
        std::cout << "No clear solution based on marginal gain, returning all solutions" << std::endl;
    }else{
        std::cerr << "Best Recombination count: " << best_r << std::endl;
    }

    return solutions;

}



void Approximator::approximate(std::vector<std::pair<std::string, std::string>> &full_ip_reads, bool diploid)
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
    int32_t count_sp_r = 0;
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

    std::vector<std::vector<std::vector<std::vector<int32_t>>>> Anchor_hits(Sp_R.size(), std::vector<std::vector<std::vector<int32_t>>>(num_walks));
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

    // print first 40 elements from Hist_kmer
    for (size_t i = 0; i < std::min((size_t)40, Hist_kmer.size()); i++) {
        fprintf(stderr, "[%d, %.2f]\n", Hist_kmer[i].multiplicity, Hist_kmer[i].freq);
    }

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
    fprintf(stderr, "Max Multiplicity: %d, Max Frequency: %.2f\n", max_multiplicity, max_freq);

    KGFitOptions opt;
    opt.max_copy = 10;
    opt.max_x_use = max_multiplicity;
    opt.u_hi = max_multiplicity;
    opt.fit_error = true;
    opt.fit_varw  = true;

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
    int32_t total_count = 0;
    
    KmerGenieDiploidLike M(res.P);
    // classify a few multiplicities
    for (int f = 1; f < max_multiplicity; f++) {
        auto r = M.classify(f);
        const char* L = (r.label==KGPosterior::HET?"HET":(r.label==KGPosterior::HOM?"HOM":"AMB"));
        if(f <= 30) std::cout << "f="<<f<<"  post[het,hom]=("<<r.p_het<<","<<r.p_hom<<")  "<<L<<"\n";
        // fill S_homo, S_hetero, S_ambiguous with
        if (r.label == KGPosterior::HOM) {
            // rev_kmer_freq
            for (auto& kmer : rev_kmer_freq[f]) {
                S_homo.insert(kmer);
                total_count++;
            }
        } else if (r.label == KGPosterior::HET) {
            for (auto& kmer : rev_kmer_freq[f]) {
                S_hetero.insert(kmer);
                total_count++;
            }
        } else {
            for (auto& kmer : rev_kmer_freq[f]) {
                S_ambiguous.insert(kmer);
                total_count++;
            }
        }
    }

    // print the fraction of S_homo, S_hetero and S_ambiguous
    fprintf(stderr, "[M::%s] Classification done. Homozygous: %.2f%%, Heterozygous: %.2f%%, Total kmers: %d\n",
        __func__,
        100.0 * S_homo.size() / total_count,
        100.0 * S_hetero.size() / total_count,
        total_count);

    // before the loop, use plain ints
    int64_t total_kmers = 0;
    int64_t count_homo  = 0;
    int64_t count_hetero= 0;
    int64_t count_ambiguous= 0;


    std::vector<bool> homo_bv(count_sp_r,0);
    std::vector<bool> het_bv(count_sp_r,0);
    #pragma omp parallel for reduction(+:total_kmers,count_homo,count_hetero,count_ambiguous) num_threads(num_threads)
    for (int32_t r = 0; r < count_sp_r; ++r) {
        auto it = Sp_R_hash.find(r);
        if (it == Sp_R_hash.end()) continue;
        uint64_t hash = it->second;
        bool is_homo = (S_homo.count(hash) != 0);
        bool is_hetero = (S_hetero.count(hash) != 0);
        bool is_ambiguous = (S_ambiguous.count(hash) != 0);

        if (is_homo + is_hetero + is_ambiguous == 0 ) continue;

        // check only one of is_homo, is_hetero, is_ambiguous is true
        if (is_homo + is_hetero + is_ambiguous != 1) { // Sanity check
            fprintf(stderr, "[M::%s] Error: k-mer %lu classified as multiple types\n",
                __func__, hash);
            // print multiple classified types as well like, homo, hetero and ambigous
            fprintf(stderr, "[M::%s] k-mer %lu classified as homo: %d, hetero: %d, ambiguous: %d\n",
                __func__, hash, is_homo, is_hetero, is_ambiguous);
            exit(1);
        }

        total_kmers++;
        if(is_homo){
            homo_bv[r] = 1;
            count_homo++;
        }else if (is_hetero){
            het_bv[r] = 1;
            count_hetero++;
        }else{
            count_ambiguous++;
        }

    }

    // now you can compute your fractions just as before
    float frac_h  = float(count_homo)   / std::max<int64_t>(1, total_kmers);
    float frac_he = float(count_hetero) / std::max<int64_t>(1, total_kmers);
    fprintf(stderr,
        "[M::%s] Phasing done. Homozygous: %.2f%%, Heterozygous: %.2f%%, Total kmers: %lld\n",
        __func__, frac_h*100, frac_he*100, total_kmers);

    // clear Anchor_hits_homo[][h]



    // END kmer classifier *****************************************************
    

    // for(int r = 0; r < Anchor_hits.size(); r++){
    //     std::cout << "\n\nr (color): " << r << std::endl;
    //     for(int h = 0; h < Anchor_hits[r].size(); h++){
    //         std::cout << "\th (haplotype path): " << h << std::endl;
    //         for(int _ = 0; _ < Anchor_hits[r][h].size(); _++){
    //             std::cout << "\t\toccurrence: " << _ << std::endl;
    //             for(int u = 0; u < Anchor_hits[r][h][_].size(); u++){
    //                 std::cout << "\t\t\t" << Anchor_hits[r][h][_][u] << std::endl;
    //             }
    //         }
    //     }
    // }

    // build expanded graph representation - weighted edges
    int32_t number_of_vertices = 0;
    for(int h = 0; h < paths.size(); h++){
        number_of_vertices += paths[h].size();
    }
    std::vector<std::vector<std::pair<int32_t, int32_t>>> expanded_graph_adj_list(2+number_of_vertices);
    std::vector<std::vector<int32_t>> vertex_to_expanded_map(n_vtx, std::vector<int32_t>(paths.size(),-1));
    std::vector<std::vector<int32_t>> expanded_vertex_to_original_map(2+number_of_vertices, std::vector<int32_t>());
    std::vector<int32_t> vertex_to_haplotype_map(2+number_of_vertices);

    int sink = expanded_graph_adj_list.size()-1;
    int32_t current_vertex = 1;
    for(int h = 0; h < paths.size(); h++){

        // from source
        expanded_graph_adj_list[0].push_back({current_vertex,0});
        for(int i = 0; i < paths[h].size(); i++){

            vertex_to_expanded_map[paths[h][i]][h] = current_vertex;

            expanded_vertex_to_original_map[current_vertex].push_back(paths[h][i]);

            vertex_to_haplotype_map[current_vertex] = h;
            if(i < paths[h].size()-1){
                expanded_graph_adj_list[current_vertex].push_back({current_vertex+1, 0});
                current_vertex++;
            }else{
                //to sink
                expanded_graph_adj_list[current_vertex].push_back({expanded_graph_adj_list.size()-1, 0});
                current_vertex++;
            }
        }
    }

    // for(int u = 0; u < expanded_graph_adj_list.size(); u++){
    //     std::cout << "\n" << u << ": ";
    //     for(auto v_w : expanded_graph_adj_list[u]){
    //         std::cout << " (" << u << ", " << v_w.first << ", " << v_w.second << ") ";
    //     }
    // }
    // std::cout << "\nvertex map" <<std::endl;
    // for(int u; u < vertex_to_expanded_map.size(); u++){
    //     std::cout <<"\n" << u << ": ";
    //     for(auto v : vertex_to_expanded_map[u]){
    //         std::cout  << v << ", ";
    //     }
    // }

    // construct recombination edges
    std::vector<std::vector<int>> vertex_w_uv (n_vtx);
    for(int u = 0; u < adj_list.size(); u++){
        vertex_w_uv[u] = std::vector<int>(adj_list[u].size(), -1);
    }

    std::vector<std::vector<int32_t>> recombination_vertex(paths.size(), std::vector<int32_t>(0));

    int recombination_edge_count = 0;
    current_vertex = expanded_graph_adj_list.size();
    for(int h = 0; h < paths.size(); h++){
       for(auto i = 0; i < paths[h].size(); i++){
            int u = paths[h][i];
            for (int j = 0; j < adj_list[u].size(); j++){
                int v = adj_list[u][j];
                //std::cout << "h: " << h << " u: " << u << " v: " << v << std::endl;
                if (i == paths[h].size()-1 || v != paths[h][i+1]){
                    if(vertex_w_uv[u][j] == -1){
                        expanded_graph_adj_list.push_back(std::vector<std::pair<int32_t, int32_t>>());
                        expanded_vertex_to_original_map.push_back(std::vector<int32_t>());
                        vertex_to_haplotype_map.push_back(-1);
                        
                        vertex_w_uv[u][j] = current_vertex;
                        current_vertex++;
  
                    }
                    //std:: cout << "recombination edge needed from " << vertex_to_expanded_map[u][h] << " to ";
                    expanded_graph_adj_list[vertex_to_expanded_map[u][h]].push_back({vertex_w_uv[u][j],1});

                    recombination_vertex[h].push_back(vertex_to_expanded_map[u][h]);
                    recombination_edge_count++;
                    if(expanded_graph_adj_list[vertex_w_uv[u][j]].empty()){
                        for(auto v_e : vertex_to_expanded_map[v]){
                            if(v_e >= 0){
                                expanded_graph_adj_list[vertex_w_uv[u][j]].push_back({v_e,0});
                                recombination_vertex[vertex_to_haplotype_map[v_e]].push_back(v_e);
                                //std::cout << v_e << ", ";
                            }
                        }
                    }
                    //std::cout << "vertex_w_uv = " << vertex_w_uv[u][j] << std::endl;
                }
            }
       }
    }
    for (auto &vec : recombination_vertex)          // iterate every haplotype
    {
        std::sort(vec.begin(), vec.end());          // ❶ order
        vec.erase(std::unique(vec.begin(), vec.end()),
                vec.end());                       // ❷ shrink – remove dups
    }
    //std::cout << "Number of recombination edges: " << recombination_edge_count << std::endl;
    // std::cout << "Recombination vertex" << std::endl;
    // for(int h =0; h < recombination_vertex.size(); h++){
    //     std::cout << h << ": " << "count: " << recombination_vertex[h].size() << " / " << paths[h].size() << ": ";
    //     // for(auto v : recombination_vertex[h]){
    //     //     std::cout << " " << v;
    //     // }
    //     std::cout << std::endl;
    // }
    // for(int u = 0; u < expanded_graph_adj_list.size(); u++){
    //     std::cout << "\n" << u << ": ";
    //     for(auto v_w : expanded_graph_adj_list[u]){
    //         std::cout << " (" << u << ", " << v_w.first << ", " << v_w.second << ") ";
    //     }
    // }
    // std::cout << std::endl;

    //create expanded graph hit paths
    
    // // For DEBUG
    // std::vector<std::vector<std::vector<std::vector<int32_t>>>> Anchor_hits_test(5, std::vector<std::vector<std::vector<int32_t>>>(num_walks));
    // Anchor_hits_test[0][0].resize(2);
    // Anchor_hits_test[0][0][0] = {1,3,6};
    // Anchor_hits_test[0][0][1] = {3,6};

    // Anchor_hits_test[0][1].resize(1);
    // Anchor_hits_test[0][1][0] = {1,3,5};

    // Anchor_hits_test[1][2].resize(1);
    // Anchor_hits_test[1][2][0] = {0,2};

    // Anchor_hits_test[1][3].resize(1);
    // Anchor_hits_test[1][3][0] = {0,2};

    // Anchor_hits_test[2][2].resize(1);
    // Anchor_hits_test[2][2][0] = {6};

    // Anchor_hits_test[3][1].resize(1);
    // Anchor_hits_test[3][1][0] = {5,6};

    // Anchor_hits_test[3][4].resize(2);
    // Anchor_hits_test[3][4][0] = {3};
    // Anchor_hits_test[3][4][1] = {2};

    // Anchor_hits_test[4][4].resize(1);
    // Anchor_hits_test[4][4][0] = {4,6};

    //Anchor_hits = Anchor_hits_test;



        /* ----------------------------------------------------------------
    * helper: push one (u→v,0) edge
    * ---------------------------------------------------------------- */
    auto addZeroEdge = [&](int u, int v)
    {
        expanded_graph_adj_list[u].push_back({v, 0});
    };

    /* ----------------------------------------------------------------
    * one record per anchor
    * ---------------------------------------------------------------- */
    struct AnchorRec {
        int startExp;                  // left-most vertex   (expanded graph)
        int endExp;                    // right-most vertex
        std::vector<int> colours;      // colour set of *this* anchor
        int nodeID;                    // super-node ID that represents it
    };


        /* ----------------------------------------------------------------
    * 1.  create a vector<AnchorRec> *per haplotype*
    * ---------------------------------------------------------------- */ 
    std::vector<std::vector<int32_t>> color(expanded_graph_adj_list.size(),std::vector<int32_t>());
    std::vector<std::vector<AnchorRec>> anchorsByHap(paths.size());   // Hₕ
    std::vector<int32_t> color_to_anchor;                // color_id -> anchor_id (a)
    color_to_anchor.reserve(Anchor_hits.size());

    int nextID = static_cast<int>(expanded_graph_adj_list.size());    // first
                                                                    // free ID
    int colourID = 0;
    for (int a = 0; a < Anchor_hits.size(); ++a)
    {
        bool new_color_used = false;
        for (int h = 0; h < Anchor_hits[colourID].size(); ++h)
        {
            for (const auto& occ : Anchor_hits[a][h])
            {
                if (occ.empty()) continue;

                new_color_used = true;
                /* original vertex IDs */
                int startOrig = occ.front();
                int endOrig   = occ.back();

                /* map to expanded-graph vertices that already exist */
                int startExp = vertex_to_expanded_map[startOrig][h];
                int endExp   = vertex_to_expanded_map[endOrig  ][h];

                /* ----------------------------------------------------------------
                * create / reuse a super-node
                * ---------------------------------------------------------------- */
                int nodeID;
                if (startExp == endExp) {                // anchor spans 1 vertex
                    nodeID = startExp;                   // reuse existing vertex
                } else {
                    /* build a fresh super-node <nextID>  and connect it         */
                    addZeroEdge(startExp, nextID);       // startExp → nextID  (0)
                    expanded_graph_adj_list.push_back({{endExp, 0}}); // nextID → endExp

                    expanded_vertex_to_original_map.push_back(std::vector<int32_t>());
                    for(auto u_original : occ){
                        expanded_vertex_to_original_map.back().push_back(u_original);

                    }

                    color.push_back({});                 // reserve colour slot
                    vertex_to_haplotype_map.push_back(h);
                    nodeID = nextID++;
                }

                /* store the record */
                anchorsByHap[h].push_back({startExp,
                                        endExp,
                                        {colourID},   // initial colour set
                                        nodeID});
            }
        }
        if(new_color_used){ 
            color_to_anchor.push_back(a);   // map color_id -> anchor_id
            colourID++;
        }
    }
    
    std::unordered_set<int> colours_set;
    for (int h = 0; h < paths.size(); ++h)
    {
        for(auto anc : anchorsByHap[h]){
            for(auto c : anc.colours){
                colours_set.insert(c);
            }
        }
    }
    std::cout << "Total color count: " << colours_set.size() << std::endl;
    auto it = std::max_element(colours_set.begin(), colours_set.end());
    if (it != colours_set.end()) {
        int max_color_val = *it;
        // use max_val
        std::cout << "Max color value: " << max_color_val << std::endl;
    }

    /* ----------------------------------------------------------------
    * 2.  sweep each haplotype to
    *     – link touching / overlapping anchors
    *     – propagate colours through containment
    * ---------------------------------------------------------------- */
    for (int h = 0; h < paths.size(); ++h)
    {
        //std::cout << "haplotype: " << h << std::endl;
        auto& vec = anchorsByHap[h];
        if (vec.empty()) continue;

        /* sort by left end */
        std::sort(vec.begin(), vec.end(),
                [](const AnchorRec& a, const AnchorRec& b)
                {   
                    if(a.startExp != b.startExp){
                        return a.startExp < b.startExp;
                    }else{
                        return a.endExp < b.endExp;
                    } 
                });

        std::vector<AnchorRec*> stk;      // “active” intervals, ordered by end

        for (auto& anc : vec)
        {
            //std::cout << "anchor: " << anc.startExp << " - " << anc.endExp << " colors:";
            // for(auto c : anc.colours){
            //     std::cout << " " << c;
            // }
            // std::cout << std::endl;
            /* ---------- pop anchors that ended before ‘anc’ starts ---------- */
            while (!stk.empty() && stk.back()->endExp < anc.startExp)
                stk.pop_back();

            /* ---------- overlap / touch  ---------- */
            if (!stk.empty() && anc.startExp <= stk.back()->endExp && stk.back()->nodeID != anc.nodeID){
                //std::cout << "Added edge: (" << stk.back()->nodeID << ", " << anc.nodeID <<", 0)" << std::endl; 
                addZeroEdge(stk.back()->nodeID , anc.nodeID);       // safe edge
            }

            /* ---------- propagate colours through *containment* ------------- */
            for (int i = static_cast<int>(stk.size()) - 1; i >= 0; --i) {
                if (anc.endExp <= stk[i]->endExp) {          // fully inside
                    for (int c : anc.colours)
                        if (std::find(stk[i]->colours.begin(),
                                    stk[i]->colours.end(),
                                    c) == stk[i]->colours.end())
                            stk[i]->colours.push_back(c);
                } else break;               // earlier items end even later
            }

            stk.push_back(&anc);            // current anchor becomes active
        }

        /* ---------- write final colour sets back to `color` ---------------- */
        for (const auto& anc : vec) {
            auto& dst = color[anc.nodeID];
            dst.insert(dst.end(), anc.colours.begin(), anc.colours.end());
            std::sort(dst.begin(), dst.end());
            dst.erase(std::unique(dst.begin(), dst.end()), dst.end());
        }
    }

    // //remove shunted simple chains
    // for haplotype path finding this is not needed
    // for (int h = 0; h < paths.size(); ++h)
    // {
    //     //std::cout << "haplotype: " << h << std::endl;
    //     auto& vec = anchorsByHap[h];
    //     if (vec.empty()) continue;

    //     for (auto& anc : vec)
    //     {
    //         // remove chain if possible
    //         if(anc.startExp != anc.endExp){
    //             bool simple = true;
    //             for(int u = anc.startExp+1; u < anc.endExp; u++){
    //                 if(expanded_graph_adj_list[u].size() > 1){
    //                     simple = false;
    //                 }
    //             }
    //             if(simple){
    //                 auto& nbrs = expanded_graph_adj_list[anc.startExp];
    //                 int v = anc.startExp+1; 
    //                 nbrs.erase( std::remove_if(nbrs.begin(), nbrs.end(),
    //                         [v](const std::pair<int,int>& p) {
    //                             return p.first == v;   // keep only p.first ≠ v
    //                             }), nbrs.end() );
    //                 for(int u = anc.startExp+1; u < anc.endExp; u++){
    //                     expanded_graph_adj_list[u].clear();
    //                     color[u].clear();
    //                 }

    //             }
    //         }
    //     }
    // }


    // run our approximation algorithm
    Graph g;
    g.adj_list = expanded_graph_adj_list;
    g.color = color;
    g.original_vertex = expanded_vertex_to_original_map;  
    g.haplotype = vertex_to_haplotype_map;  
    

    /* compactify and keep the mapping */
    std::cout << "Compactifying..." << std::endl;
    sink = g.compactify(sink);

    std::cout << "Topologically ordering..." << std::endl;

    g.topologically_reorder(sink);
    //g.print();
    std::cout << "Topologically reordered" << std::endl;


    // DP approximation approach
    if(!diploid){
        // HAPLOID CASE

        std::vector<int> dp_path = dp_approximation_solver(g, recombination_limit);
        std::string dp_output;
        for(auto u_original : dp_path){
            dp_output += node_seq[u_original];
        }
        //std::cout << "DP output: " << dp_output << "\n" << std::endl;

        // write haplotype as to a file as fasta from the path
        std::ofstream dp_file_stream(hap_file, std::ios::out);
        dp_file_stream << ">" << "dp_sol" << " LN:" << dp_output.size() << std::endl;
        // write the path_str to the file 80 characters per line
        for (size_t i = 0; i < dp_output.size(); i += 80) {
            dp_file_stream << dp_output.substr(i, 80) << std::endl;
        }
        dp_file_stream.close();



    }else{
        // DIPLOID CASE

        std::vector<bool> color_homo_bv(colours_set.size(), 0);

        for(const auto& c : colours_set){
            int anchor = color_to_anchor[c];
            if(homo_bv[anchor]){
                color_homo_bv[c] = 1;
            }
        }


        // count number of weighted edges
        int weighted_edge_count = 0;
        for(int i = 0; i < g.adj_list.size(); i++){
            for(const auto& [u,w] : g.adj_list[i]){
                if(w > 0){
                    weighted_edge_count++;
                }
            }
        }


        std::cout << "Number of weighted (recombination) edges: " << weighted_edge_count << std::endl;

        // Exact (ILP) diploid solution under the same budget:
        //auto [opt_r1, opt_r2, opt_s1, opt_s2] = optimal_diploid_solver(g, recombination_limit, color_homo_bv, color_het_bv);

        std::cout << "Leveling" << std::endl;
        int width = g.strict_bfs_levelize_and_reorder();
        //g.print();
        std::cout << "number of levels = " << g.vertices_in_level.size() << std::endl;
        std::cout << "max level width = " << width << std::endl;

        int sum_of_level_widths = 0;
        for(auto l : g.vertices_in_level){
            sum_of_level_widths += l.size();
        }
        std::cout << "Average level width = " << sum_of_level_widths/(float)g.vertices_in_level.size() << std::endl;

        // on each level l, if a color appears on hap h, 
        // make sure it appears for all vertices on in hap h on that level
        std::cout << "Copying het colors across hap vertices on level..." << std::endl;
        
        int lev = 0;
        for(auto l : g.vertices_in_level){
            std::cout << "\rlevel: " << lev++;
            for(auto v : l){
                for(auto c : g.color[v]){
                    if(!homo_bv[c]){
                        for(auto u : l){
                            if(g.haplotype[u] == g.haplotype[v]){
                                g.color[v].push_back(c);
                            }
                        }
                    }
                }
            }
        }
        std::cout << std::endl;


        std::vector<std::tuple<int, int, std::string, std::string>> solutions = diploid_dp_approximation_solver(g, recombination_limit, color_homo_bv);
        
        for(const auto& [r1, r2, s1, s2] : solutions){
            std::cout << r1 << ", " << r2 << ", " << s1.length() << ", " << s2.length() << std::endl;
        }

        if(solutions.size() == 1){
            // int best_sol = evaluate_diploid_solution_based_on_kmers(solutions, validation_reads);
            // std::cout << "Best sol based validation reads:" 
            //     << " r1 = " << std::get<0>(solutions[best_sol]) 
            //     << " r2 = " << std::get<1>(solutions[best_sol]) 
            //     << std::endl;

            // const auto& [r1, r2, s1, s2] = solutions.at(best_sol);

            const auto& [r1, r2, s1, s2] = solutions.at(0);
            std::ofstream dp_file_stream(hap_file, std::ios::out);
            dp_file_stream << ">" << "sol_1" << " bp:" << s1.size() << std::endl;
            // write the path_str to the file 80 characters per line
            for (size_t i = 0; i < s1.size(); i += 80) {
                dp_file_stream << s1.substr(i, 80) << std::endl;
            }
            dp_file_stream << ">" << "sol_2" << " bp:" << s2.size() << std::endl;
            // write the path_str to the file 80 characters per line
            for (size_t i = 0; i < s2.size(); i += 80) {
                dp_file_stream << s2.substr(i, 80) << std::endl;
            }
            dp_file_stream.close();
        }else{
            std::cout << "No clear solution, output file not written." << std::endl;
        }
    }
}

