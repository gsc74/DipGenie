#pragma once
#include "GindexGraph.hpp"
#include "PairedReadBinning.hpp"
#include <functional>

namespace hic {
struct GraphInput {
    std::vector<ComponentReference> references;
    // Shards are loaded only after binning; references are read first.
    std::function<NativeGraph(size_t)> load;
};
GraphInput graph_input(const std::map<std::string,std::string>&args,int threads=1);
void write_graph_index(const GraphInput&,const std::filesystem::path&);
}
