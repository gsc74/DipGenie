#pragma once
#include "ReadLibraries.hpp"
namespace hic {
// Single inference objective/search for SR, Hi-C, and their combination.
int infer_component(std::map<std::string,std::string> args,
                    std::vector<ReadLibrary> libraries, NativeGraph graph);
}
