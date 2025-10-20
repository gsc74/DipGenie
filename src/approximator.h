#pragma once
#include "ExpandedGraph.hpp"
#include "solver.h"


class Approximator : public Solver {
    public:

        // Data Structures

        struct AnchorRec {
            int startOrg;
            int endOrg;
            int startExp;                  // left-most vertex   (expanded graph)
            int endExp;                    // right-most vertex  (expanded graph)
            std::vector<int> colours;      // colour set of *this* anchor
            int nodeID;                    // super-node ID that represents it
        };

        // Constructor
        explicit Approximator(gfa_t *g);	// This is constructor

        // Functions
        void solve(std::vector<std::pair<std::string, std::string>> &ip_reads, bool diploid);
        std::vector<int> dp_approximation_solver(ExpandedGraph g, int R);
        std::vector<std::tuple<int, int, std::string, std::string>> diploid_dp_approximation_solver(ExpandedGraph g, int R, std::vector<bool> color_homo_bv, std::vector<std::vector<AnchorRec>> anchorsByHap);
};
