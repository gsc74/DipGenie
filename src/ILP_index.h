#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include "PHIpriv.h"
#include <assert.h>
#include <math.h>
#include "khashl.h"
#include "kalloc.h"
#include <zlib.h>
#include "kseq.h"
#include "kvec-km.h"
#include "sys.h"

#include "Classifier.hpp"
#include "Fitter.hpp"

// CPP
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <queue>
#include <stack>
#include <functional>
#include <tuple>
#include <map>
#include <sstream>
#include <unordered_set>
#include <cctype>
#include <string>

#include "gurobi_c++.h"
#include "MurmurHash3.h"
#include "solver.h"


class ILP_index : public Solver {
    public:

        bool is_low_cov;
        int32_t top_k;

        // Constructor
        explicit ILP_index(gfa_t *g);	// This is constructor

        // Functions
        void solve(std::vector<std::pair<std::string, std::string>> &ip_reads);

};
