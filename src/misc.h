#pragma once
#include <vector>
#include <string>
#include "PHIpriv.h"
#include "ksort.h"

// declare only
std::vector<std::pair<int,int>>
find_bridges(std::vector<std::vector<int>>& adj_list);

//bool  compare_T(const mg128_t& a, const mg128_t& b);
void get_hap_name(char* gfa_name, char* reads_name, std::string& hap_name);
std::string reverse_strand_(std::string seq);

// globals – declare with extern
extern int    mg_verbose;
extern int    mg_dbg_flag;
extern double mg_realtime0;
