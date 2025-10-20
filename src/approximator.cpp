#include "approximator.h"
#include "ExpandedGraph.hpp"
#include <cstdint>
#include <functional>
#include <numeric>
#include <random>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <chrono>
#include <sstream>

//#include "edlib.h"
#include "misc.h"
#include "Classifier.hpp"
#include "Fitter.hpp"


// Constructor
Approximator::Approximator(gfa_t *g) : Solver(g) {}

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


// treat all colors as distinct
std::vector<int> Approximator::dp_approximation_solver(ExpandedGraph g, int R) {

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
inline std::vector<T> intersection_sorted(const std::vector<T>& a,
                                          const std::vector<T>& b) {
    std::vector<T> out;
    out.reserve(std::min(a.size(), b.size())); // upper bound
    std::size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] < b[j]) {
            ++i;
        } else if (b[j] < a[i]) {
            ++j;
        } else {
            out.push_back(a[i]); // equal → in intersection
            ++i; ++j;
        }
    }
    out.shrink_to_fit();
    return out;
}

template <class T>
inline std::vector<T> symdiff_sorted(const std::vector<T>& a,
                                     const std::vector<T>& b) {
    std::vector<T> out;
    out.reserve(a.size() + b.size()); // upper bound
    std::size_t i = 0, j = 0;
    while (i < a.size() || j < b.size()) {
        if (j == b.size() || (i < a.size() && a[i] < b[j])) {
            out.push_back(a[i]); // only in a
            ++i;
        } else if (i == a.size() || b[j] < a[i]) {
            out.push_back(b[j]); // only in b
            ++j;
        } else {
            // equal → not in symmetric diff
            ++i; ++j;
        }
    }
    out.shrink_to_fit();
    return out;
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



std::vector<std::tuple<int, int, std::string, std::string>> Approximator::diploid_dp_approximation_solver(ExpandedGraph g, int R, std::vector<bool> color_homo_bv, std::vector<std::vector<AnchorRec>> anchorsByHap) {
    
    
    // Build vertex -> position-in-its-level map once
    const int L = (int)g.vertices_in_level.size();
    std::vector<int> pos_in_level(g.adj_list.size(), -1);
    for (int l = 0; l < L; ++l) {
        const auto& Lv = g.vertices_in_level[l];
        for (int i = 0; i < (int)Lv.size(); ++i) pos_in_level[ Lv[i] ] = i;
    }

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
                                    dst.s_het = src.s_het + s_hets[idx];
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

    std::vector<std::unordered_map<int,int>> p1_color_freq_by_r;
    std::vector<std::unordered_map<int,int>> p2_color_freq_by_r;

    std::vector<std::vector<int>> p1_color_by_r;
    std::vector<std::vector<int>> p2_color_by_r;
    for( int r = 0; r <= R; r++ ) {


        std::cout << "Computing haplotype paths for r = " << r << std::endl;
        auto& sink_dp = at3(dp_cur, k_sink, 0, 0, r); // corresponds to single sink vertex
        int r1 = sink_dp.weighted_p1_edges.size()-1;
        int r2 = sink_dp.weighted_p2_edges.size()-1;
        std::cerr << "Recombinations in P1: " << r1 << std::endl;
        std::cerr << "Recombinations in P2: " << r2 << std::endl;

        int start_exp = g.vertices_in_level.at(0).at(0);
        int h; 

        std::unordered_map<int,int> p1_color_freq;
        std::vector<int> p1_colors;
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

            for(auto a : anchorsByHap[h]){
                if(a.startOrg > start_org && a.endOrg < end_org){
                    for(auto c : a.colours){
                        if(p1_color_freq.find(c) == p1_color_freq.end()){
                            p1_color_freq[c] = 1;
                            p1_colors.push_back(c);
                        }else{
                            p1_color_freq[c] += 1;
                        }
                    }
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

        std::unordered_map<int,int> p2_color_freq;
        std::vector<int> p2_colors;
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

            for(auto a : anchorsByHap[h]){
                if(a.startOrg > start_org && a.endOrg < end_org){
                    for(auto c : a.colours){
                        if(p2_color_freq.find(c) == p2_color_freq.end()){
                            p2_color_freq[c] = 1;
                            p2_colors.push_back(c);
                        }else{
                            p2_color_freq[c] += 1;
                        }
                    }
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
        p1_color_freq_by_r.push_back(p1_color_freq);
        p2_color_freq_by_r.push_back(p2_color_freq);
        p1_color_by_r.push_back(p1_colors);
        p2_color_by_r.push_back(p2_colors);
    }

    // Delta score approach
    // computes score and approx cert.
    std::vector<int> obj_by_r;

    for( int r = 0; r <= R; r++ ) {

        std::cout << "Computing score for r = " << r << std::endl;

        
        auto& sink_dp = at3(dp_cur, k_sink, 0, 0, r); // corresponds to single sink vertex

        int s_het = sink_dp.s_het;

        std::vector<int> path_1_homo_colors;
        std::vector<int> path_1_hetero_colors;

        for(auto c : p1_color_by_r[r]){
            if(color_homo_bv[c]){
                path_1_homo_colors.push_back(c);
            }else{
                path_1_hetero_colors.push_back(c);
            }
        }

        std::vector<int> path_2_homo_colors;
        std::vector<int> path_2_hetero_colors;

        for(auto c : p2_color_by_r[r]){
            if(color_homo_bv[c]){
                path_2_homo_colors.push_back(c);
            }else{
                path_2_hetero_colors.push_back(c);
            }
        }


        dedup_inplace(path_1_homo_colors);
        dedup_inplace(path_1_hetero_colors);
        dedup_inplace(path_2_homo_colors);
        dedup_inplace(path_2_hetero_colors);
    
        auto intersection = intersection_sorted(path_1_homo_colors, path_2_homo_colors);
        auto symdif = symdiff_sorted(path_1_hetero_colors, path_2_hetero_colors);
        int intersection_count = intersection.size();
        int symdiff_count = symdif.size();

        int m_G_hom = 0;
        int m_G_het = 0;
        auto p1_color_freq = p1_color_freq_by_r[r];
        auto p2_color_freq = p2_color_freq_by_r[r];

        for(auto c : intersection){
            auto it1 = p1_color_freq.find(c);
            auto it2 = p2_color_freq.find(c);
            int k1 = (it1 == p1_color_freq.end()) ? 0 : it1->second;
            int k2 = (it2 == p2_color_freq.end()) ? 0 : it2->second;
            m_G_hom += (k1 >= k2 ? k1 : k2);
        }

        for(auto c : symdif){
            auto it1 = p1_color_freq.find(c);
            auto it2 = p2_color_freq.find(c);
            int k1 = (it1 == p1_color_freq.end()) ? 0 : it1->second;
            int k2 = (it2 == p2_color_freq.end()) ? 0 : it2->second;
            m_G_het += k1 + k2;
        }


        // approximation cert.
        float m_G_hom_avg = m_G_hom / (float) intersection_count;
        float m_G_het_avg = m_G_het / (float) symdiff_count;
        float m_bar = std::max(m_G_hom_avg, m_G_het_avg);

        int loss_het = s_het-m_G_het;
        float additive_term = loss_het / (float)m_G_het_avg;

        int obj = intersection_count + symdiff_count;
        std::cout << "r: " << r << " obj: " << obj << std::endl;
        obj_by_r.push_back(obj);

        float opt_obj_upper_bound =  m_bar*(obj + additive_term);


        std::cout << "Approximation certificate: multiplicative factor: " << opt_obj_upper_bound/(float)obj  << std::endl;    

    }

    int best_r = 0;
    for (size_t r = 0; r + 1 < obj_by_r.size(); ++r) {
        int delta = obj_by_r[r + 1] - obj_by_r[r];
        std::cout << "r: " << r << " -> " << r + 1 << ", Δobj: " << delta << std::endl;
        if(delta <= 0){
            best_r = r;
            break;
        }
    }

    std::cerr << "Best Recombination count: " << best_r << std::endl;
    auto best_sol = solutions.at(best_r);
    solutions.clear();
    solutions.push_back(best_sol);

    return solutions;

}



void Approximator::solve(std::vector<std::pair<std::string, std::string>> &full_ip_reads, bool diploid)
{
    // View anchors for debug

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
                    vertex_to_haplotype_map.push_back(-1);
                    nodeID = nextID++;
                }

                /* store the record */
                anchorsByHap[h].push_back({
                                        startOrig,
                                        endOrig,
                                        startExp,
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

    //remove shunted simple chains (currently just removes shunted sigle edges)
    // for haplotype path finding this is not needed
    for (int h = 0; h < paths.size(); ++h)
    {
        //std::cout << "haplotype: " << h << std::endl;
        auto& vec = anchorsByHap[h];
        if (vec.empty()) continue;

        for (auto& anc : vec)
        {
            // remove chain if possible
            if(anc.startExp != anc.endExp){
                bool simple = true;
                for(int u = anc.startExp+1; u < anc.endExp; u++){
                    if(expanded_graph_adj_list[u].size() > 0){ //only removes single edge shunts
                        simple = false;
                    }
                }
                if(simple){
                    auto& nbrs = expanded_graph_adj_list[anc.startExp];
                    int v = anc.startExp+1; 
                    nbrs.erase( std::remove_if(nbrs.begin(), nbrs.end(),
                            [v](const std::pair<int,int>& p) {
                                return p.first == v;   // keep only p.first ≠ v
                                }), nbrs.end() );
                    for(int u = anc.startExp+1; u < anc.endExp; u++){
                        expanded_graph_adj_list[u].clear();
                        color[u].clear();
                    }

                }
            }
        }
    }


    // run our approximation algorithm
    ExpandedGraph g;
    g.adj_list = expanded_graph_adj_list;
    g.color = color;
    g.original_vertex = expanded_vertex_to_original_map;  
    g.haplotype = vertex_to_haplotype_map;  
    

    /* compactify and keep the mapping */
    //std::cout << "Compactifying..." << std::endl;
    //sink = g.compactify(sink);

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


        std::vector<std::tuple<int, int, std::string, std::string>> solutions = diploid_dp_approximation_solver(g, recombination_limit, color_homo_bv, anchorsByHap);
        
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
            std::cout << "No solution reported, output file not written." << std::endl;
        }
    }
}

