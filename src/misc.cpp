// misc.cpp
#include "misc.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

// ------------------- globals (define once here) -------------------
int    mg_verbose   = 1;
int    mg_dbg_flag  = 0;
double mg_realtime0 = 0.0;

// ------------------- local DFS helper (internal linkage) ---------
namespace {
void dfs(int u,
         std::vector<int>& visited,
         std::vector<int>& parent,
         std::vector<int>& low,
         std::vector<int>& disc,
         int& time,
         std::vector<std::vector<int>>& adj_list,
         std::vector<std::pair<int, int>>& bridges) {
  visited[u] = 1;
  low[u] = disc[u] = ++time;

  for (int v : adj_list[u]) {
    if (!visited[v]) {
      parent[v] = u;
      dfs(v, visited, parent, low, disc, time, adj_list, bridges);
      low[u] = std::min(low[u], low[v]);
      if (low[v] > disc[u]) {
        bridges.push_back({u, v});
      }
    } else if (v != parent[u]) {
      low[u] = std::min(low[u], disc[v]);
    }
  }
}
}  // namespace

// ------------------- public functions (definitions) ---------------

std::vector<std::pair<int, int>> find_bridges(std::vector<std::vector<int>>& adj_list) {
  std::vector<std::pair<int, int>> bridges;
  int n = static_cast<int>(adj_list.size());
  std::vector<int> visited(n, 0);
  std::vector<int> parent(n, -1);
  std::vector<int> low(n, 0);
  std::vector<int> disc(n, 0);
  int time = 0;

  std::cerr << "n = " << n << std::endl;

  for (int i = 0; i < n; ++i) {
    if (!visited[i]) {
      dfs(i, visited, parent, low, disc, time, adj_list, bridges);
    }
  }

  // reverse the order of the bridges
  std::reverse(bridges.begin(), bridges.end());

  return bridges;
}

// bool compare_T(const mg128_t& a, const mg128_t& b) {
//   return a.x < b.x;
// }

void get_hap_name(char* gfa_name, char* reads_name, std::string& hap_name) {
  // split gfa_name to get the hap_name by last '/' and extension
  std::string gfa_file = gfa_name;
  size_t found = gfa_file.find_last_of("/\\");
  if (found != std::string::npos) {
    hap_name = gfa_file.substr(found + 1);
  } else {
    hap_name = gfa_file;
  }
  found = hap_name.find_last_of(".");
  if (found != std::string::npos) {
    hap_name = hap_name.substr(0, found);
  }

  // split reads_name similarly and append
  std::string reads_file = reads_name;
  found = reads_file.find_last_of("/\\");
  if (found != std::string::npos) {
    hap_name += "_";
    hap_name += reads_file.substr(found + 1);
  } else {
    hap_name += "_";
    hap_name += reads_file;
  }
  found = hap_name.find_last_of(".");
  if (found != std::string::npos) {
    hap_name = hap_name.substr(0, found);
  }
}

std::string reverse_strand_(std::string seq) {
  std::string rev_seq;
  rev_seq.reserve(seq.size());
  for (int i = static_cast<int>(seq.size()) - 1; i >= 0; --i) {
    char c = seq[i];
    if (c == 'A' || c == 'a') rev_seq += 'T';
    else if (c == 'T' || c == 't') rev_seq += 'A';
    else if (c == 'C' || c == 'c') rev_seq += 'G';
    else if (c == 'G' || c == 'g') rev_seq += 'C';
    else rev_seq += c;
  }
  return rev_seq;
}

// ------------------- ksort / radix instantiations (once) ----------
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mg128_t, sort_key_128x, 8)
KSORT_INIT_GENERIC(uint32_t)
