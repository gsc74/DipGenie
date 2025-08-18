#pragma once

#include <vector>
#include <utility>
#include <ostream>
#include <iostream>
#include <queue>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <cassert>

class Graph {
public:
    std::vector<int> level;
    std::vector<int> haplotype;
    std::vector<std::vector<int>> color;
    std::vector<std::vector<int>> original_vertex;
    std::vector<std::vector<int>> vertices_in_level;


    // adj_list[u] is a vector of (v, weight) pairs
    std::vector<std::vector<std::pair<int,int>>> adj_list;

    // ----- read_from_file -----
    // Format:
    //   n
    //   color[0]
    //   …
    //   color[n-1]
    //   then one line per vertex u:
    //     u v1 w1 v2 w2 …   (zero or more neighbor/weight pairs)
    //
    // e.g. “3 2 1 7 0” → from 3→2 weight=1, and 3→7 weight=0
    // static Graph read_from_file(const std::string& path) {
    //     std::ifstream in(path);
    //     if (!in) throw std::runtime_error("Unable to open file " + path);

    //     std::size_t n;
    //     if (!(in >> n))
    //         throw std::runtime_error("Missing vertex count");
    //     Graph g;
    //     g.color.resize(n);
    //     g.adj_list.resize(n);

    //     // read n colors
    //     for (std::size_t i = 0; i < n; ++i) {
    //         if (!(in >> g.color[i]))
    //             throw std::runtime_error("Missing color for vertex " + std::to_string(i));
    //     }
    //     in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    //     // now read adjacency-lines
    //     std::string line;
    //     while (std::getline(in, line)) {
    //         if (line.empty()) continue;
    //         std::istringstream iss(line);

    //         int u;
    //         if (!(iss >> u)) continue;    // skip malformed

    //         int v, w;
    //         while (iss >> v >> w) {
    //             g.adj_list[u].emplace_back(v, w);
    //         }
    //         // any leftover token is ignored
    //     }

    //     return g;
    // }


    // ----- in-place Kahn topological reorder -----
    void topologically_reorder(int sink) {

        const std::size_t n = adj_list.size();
        std::vector<int> indeg(n, 0);
        for (auto& nbrs : adj_list)
            for (auto& [v, w] : nbrs)
                ++indeg[v];

        std::queue<int> q;
        for (std::size_t v = 0; v < n; ++v)
            if (indeg[v] == 0 && (int)v != sink)        // ① never push sink now
                q.push((int)v);

        bool sink_ready = (indeg[sink] == 0);           // store flag instead
        std::vector<int> order;  order.reserve(n);

        while (!q.empty() || sink_ready) {
            int u;
            if (!q.empty()) {                           // ② process queue first
                u = q.front();  q.pop();
            } else {                                    // queue empty → only sink left
                u = sink;
                sink_ready = false;
            }

            order.push_back(u);

            for (auto& [v, w] : adj_list[u]) {
                if (--indeg[v] == 0) {
                    if (v == sink) sink_ready = true;   // ③ mark sink ready, don't push
                    else           q.push(v);
                }
            }
        }

        if (order.size() != n)
            throw std::runtime_error("Graph contains a cycle; topological order impossible");

        // build mapping old->new
        std::vector<int> new_idx(n);
        for (std::size_t i = 0; i < n; ++i)
            new_idx[ order[i] ] = (int)i;

        // reorder colors, original vertex, haplotype, and (if present) level
        std::vector<std::vector<int>> new_color(n);
        std::vector<std::vector<int>> new_original_vertex(n);
        std::vector<int>              new_haplotype(n);
        for (std::size_t i = 0; i < n; ++i){
            new_color[i]           = color[ order[i] ];
            new_original_vertex[i] = original_vertex[ order[i] ];
            new_haplotype[i]       = haplotype[ order[i] ];
        }
        color.swap(new_color);
        original_vertex.swap(new_original_vertex);
        haplotype.swap(new_haplotype);

        if (level.size() == n) {
            std::vector<int> new_level(n);
            for (std::size_t i = 0; i < n; ++i)
                new_level[i] = level[ order[i] ];
            level.swap(new_level);
        }

        // reorder & remap adjacency list
        std::vector<std::vector<std::pair<int,int>>> new_adj(n);
        for (std::size_t old_u = 0; old_u < n; ++old_u) {
            int u = new_idx[old_u];
            for (auto& [old_v, w] : adj_list[old_u]) {
                int v = new_idx[old_v];
                new_adj[u].emplace_back(v, w);
            }
        }
        adj_list.swap(new_adj);
    }

    // ----- print -----
    void print(std::ostream& out = std::cout) const {
        for (std::size_t u = 0; u < adj_list.size(); ++u) {
            out << 'u' << u << " (color =";
            for (auto c : color[u]) {
                out << ' ' << c;
            }
            out << ",\t\toriginal vertex =";
            for (auto v : original_vertex[u]) {
                out << ' ' << v;
            }
            if (haplotype.size() == adj_list.size()) {
                out << ",\t\thap = " << haplotype[u];
            }
            if (level.size() == adj_list.size()) {
                out << ",\t\tlevel = " << level[u];
            }
            out << ")\t\t->";
            for (auto& [v, w] : adj_list[u])
                out << ' ' << v << '(' << w << ')';
            out << '\n';
        }
    }


/* ------------------------------------------------------------------
 *  Graph::compactify()
 * ------------------------------------------------------------------*/
int compactify (int old_sink)
{
    const std::size_t n = adj_list.size();

    /* ---------- 1. degree vectors -------------------------------- */
    std::vector<int> indeg (n,0), outdeg (n,0);
    std::vector<int> indeg0(n,0), outdeg0(n,0);   // 0-weight only

    for (std::size_t u = 0; u < n; ++u)
        for (auto [v,w] : adj_list[u]) {
            ++outdeg[u]; ++indeg[v];
            if (w == 0) { ++outdeg0[u]; ++indeg0[v]; }
        }

    /* ---------- 2. new graph containers -------------------------- */
    std::vector<std::vector<int>>                new_color;
    std::vector<std::vector<int>>                new_orig;
    std::vector<std::vector<std::pair<int,int>>> new_adj;
    std::vector<int>                             new_hap;   // <-- carry haplotype
    new_color.reserve(n);
    new_orig .reserve(n);
    new_adj  .reserve(n);
    new_hap  .reserve(n);

    /* ---------- 3. helpers --------------------------------------- */
    std::vector<int>  id_map (n, -1);   // old id → new id (−1 = not yet)
    std::vector<char> done   (n, 0);    // 0 = not processed, 1 = finished
    std::vector<char> swallowed(n,0);   // interior vertices already removed

    auto add_vertex = [&](std::size_t old_id) -> int {
        int new_id           = static_cast<int>(new_adj.size());
        id_map[old_id]       = new_id;
        new_adj .emplace_back();
        new_color.emplace_back( color[old_id] );
        new_orig .emplace_back( original_vertex[old_id] );
        new_hap  .emplace_back( haplotype[old_id] );   // <-- set haplotype for the kept vertex
        return new_id;
    };

    /* ---------- 4. sweep every vertex ---------------------------- */
    for (std::size_t u0 = 0; u0 < n; ++u0)
    {
        if (done[u0]) continue;                          // already processed

        /* decide whether u0 must be *kept* as an explicit vertex     */
        bool keep =
              !color[u0].empty()           ||
              indeg0[u0] != indeg[u0]      || outdeg0[u0] != outdeg[u0] ||
              indeg0[u0] != 1              || outdeg0[u0] != 1;

        if (!keep) continue;                             // interior chain node

        /* obtain (or reuse) new id for u0 */
        int new_u = (id_map[u0] != -1) ? id_map[u0] : add_vertex(u0);
        done[u0] = 1;

        for (auto [v, w] : adj_list[u0])
        {
            /* weight-1 (recombination) edges are copied verbatim */
            if (w != 0) {
                int nv = (id_map[v] != -1) ? id_map[v] : add_vertex(v);
                new_adj[new_u].emplace_back(nv, w);
                continue;
            }

            /* -------- follow the 0-weight chain ---------------- */
            std::size_t cur = v;
            while (!swallowed[cur] &&
                   color[cur].empty()              &&
                   indeg0[cur] == 1 && outdeg0[cur] == 1 &&
                   indeg [cur] == 1 && outdeg [cur] == 1)
            {
                swallowed[cur] = 1;                          // consume
                new_orig[new_u].insert(new_orig[new_u].end(),
                                       original_vertex[cur].begin(),
                                       original_vertex[cur].end());
                cur = adj_list[cur][0].first;                // unique 0-edge
            }

            int nv = (id_map[cur] != -1) ? id_map[cur] : add_vertex(cur);
            new_adj[new_u].emplace_back(nv, 0);              // still weight 0
        }
    }

    /* ---------- 5. swap compacted graph in ---------------------- */
    adj_list.swap(new_adj);
    color    .swap(new_color);
    original_vertex.swap(new_orig);
    haplotype.swap(new_hap);                                  // <-- swap haplotype, too

    /* sanity check */
    assert(adj_list.size() <= n && "compactify must not increase |V|");

    return id_map[old_sink];
}


int strict_bfs_levelize_and_reorder()
{


    const std::size_t n = adj_list.size();
        std::vector<int> indeg(n, 0);
        for (auto& nbrs : adj_list)
            for (auto& [v, w] : nbrs)
                ++indeg[v];

    int source = -1;
    for (std::size_t v = 0; v < n; ++v)
        if (indeg[v] == 0){
            if(source == -1){
                source  = v;
            }else{
                std::cout << "Uh oh, multiple potential sources found while leveling" << std::endl;
                exit(-1);
            }
        }
    
    const int n0 = static_cast<int>(adj_list.size());
    if (n0 == 0) return 0;
    if (source < 0 || source >= n0) throw std::runtime_error("bad source index");

    // ---------- 1) BFS from source: shortest hop distances ----------
    std::vector<int> dist(n0, -1);
    std::queue<int> q;
    dist[source] = 0;
    q.push(source);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (auto [v, w] : adj_list[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }
    }

    // ---------- 2) Topological order (for DAG; throw if cycle) ----------
    std::queue<int> qk;
    for (int v = 0; v < n0; ++v) if (indeg[v] == 0) qk.push(v);

    std::vector<int> topo; topo.reserve(n0);
    while (!qk.empty()) {
        int u = qk.front(); qk.pop();
        topo.push_back(u);
        for (auto [v, w] : adj_list[u]) {
            if (--indeg[v] == 0) qk.push(v);
        }
    }
    if ((int)topo.size() != n0)
        throw std::runtime_error("Graph contains a cycle; strict leveling requires a DAG");

    // ---------- 3) Seed levels from BFS; then lift to make edges forward ----------
    // Unreached vertices start at 0; then lifting below will place them.
    std::vector<int> lvl(n0, 0);
    for (int v = 0; v < n0; ++v) if (dist[v] >= 0) lvl[v] = dist[v];

    // Longest-path style relaxation in topo order:
    for (int u : topo) {
        for (auto [v, w] : adj_list[u]) {
            if (lvl[v] <= lvl[u]) lvl[v] = lvl[u] + 1;
        }
    }

    // ---------- 4) Insert dummy vertices for edges that skip levels ----------
    // Build new arrays; start by copying originals (original nodes keep their data/levels).
    std::vector<std::vector<std::pair<int,int>>> next_adj(n0);
    std::vector<std::vector<int>>                next_color = color;
    std::vector<std::vector<int>>                next_orig  = original_vertex;
    std::vector<int>                             next_lvl   = lvl;
    std::vector<int>                             next_hap   = haplotype;   // <-- carry haplotype

    auto add_dummy = [&](int new_level, int hap) -> int {  // <-- pass hap to set for dummy
        int id = static_cast<int>(next_adj.size());
        next_adj.emplace_back();
        next_color.emplace_back();     // empty
        next_orig.emplace_back();      // empty
        next_lvl.push_back(new_level);
        next_hap.push_back(hap);       // <-- set haplotype for dummy (e.g., inherit from u)
        return id;
    };

    for (int u = 0; u < n0; ++u) {
        for (auto [v, w] : adj_list[u]) {
            int gap = next_lvl[v] - next_lvl[u] - 1;
            if (gap <= 0) {
                // gap==0: edge goes to immediate next or same/earlier level;
                // we already lifted levels so same/earlier won't happen; keep as-is.
                next_adj[u].emplace_back(v, w);
            } else {
                // create chain of 'gap' dummies at consecutive levels
                int prev = u;
                for (int step = 1; step <= gap; ++step) {
                    int dummy = add_dummy(next_lvl[u] + step, haplotype[u]); // <-- inherit hap from u
                    // carry original weight on the first segment; remaining are 0
                    next_adj[prev].emplace_back(dummy, (step == 1 ? w : 0));
                    prev = dummy;
                }
                // final segment to the original target with weight 0
                next_adj[prev].emplace_back(v, 0);
            }
        }
    }

    // Swap in strictly leveled graph (but not yet ordered by level)
    adj_list.swap(next_adj);
    color   .swap(next_color);
    original_vertex.swap(next_orig);
    level   .swap(next_lvl);
    haplotype.swap(next_hap);                                  // <-- swap haplotype with leveled graph

    // ---------- 5) Order vertices by (level, stable by old id) ----------
    const int n1 = static_cast<int>(adj_list.size());
    std::vector<int> order(n1);
    for (int i = 0; i < n1; ++i) order[i] = i;

    std::stable_sort(order.begin(), order.end(),
        [&](int a, int b){
            if (level[a] != level[b]) return level[a] < level[b];
            return a < b;  // stable tie-breaker
        });

    // Compute max width while we’re at it
    int max_level = 0;
    for (int v = 0; v < n1; ++v) if (level[v] > max_level) max_level = level[v];
    std::vector<int> width(max_level + 1, 0);
    for (int v = 0; v < n1; ++v) ++width[level[v]];
    int max_width = 0;
    for (int w : width) if (w > max_width) max_width = w;

    // Build mapping old->new
    std::vector<int> new_id(n1, -1);
    for (int i = 0; i < n1; ++i) new_id[ order[i] ] = i;

    // Reorder color/original/level/haplotype
    std::vector<std::vector<int>> new_color(n1);
    std::vector<std::vector<int>> new_orig (n1);
    std::vector<int>              new_level(n1);
    std::vector<int>              new_hap(n1);   // <-- reorder haplotype too
    for (int i = 0; i < n1; ++i) {
        int old = order[i];
        new_color[i] = std::move(color[old]);
        new_orig [i] = std::move(original_vertex[old]);
        new_level[i] = level[old];
        new_hap  [i] = haplotype[old];           // <-- carry hap in the same reorder
    }
    color.swap(new_color);
    original_vertex.swap(new_orig);
    level.swap(new_level);
    haplotype.swap(new_hap);                      // <-- swap in reordered haplotype

    // Remap adjacency to new ids
    std::vector<std::vector<std::pair<int,int>>> new_adj(n1);
    for (int old_u = 0; old_u < n1; ++old_u) {
        int u = new_id[old_u];
        for (auto [old_v, w] : adj_list[old_u]) {
            int v = new_id[old_v];
            new_adj[u].emplace_back(v, w);
        }
    }
    adj_list.swap(new_adj);

    // ---------- 7) Build per-level buckets in topological order ----------
    vertices_in_level.clear();
    vertices_in_level.resize(max_level + 1);
    for (int u = 0; u < n1; ++u) {
        vertices_in_level[level[u]].push_back(u);
    }
    return max_width;
}

//
};
