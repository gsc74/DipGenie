#pragma once
#include <cstdlib>
#include <string>
// mkdtemp pattern under $TMPDIR (default /tmp).
inline std::string test_temp_pattern(const char*prefix) {
    const char*dir=std::getenv("TMPDIR");
    return std::string(dir&&*dir?dir:"/tmp")+"/"+prefix+"-XXXXXX";
}
