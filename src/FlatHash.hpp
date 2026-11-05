#pragma once
#include "../third_party/unordered_dense/ankerl/unordered_dense.h"

namespace hic {
// Only lookup/storage changes. These hashes never select minimizers or stand
// in for key equality. Do not retain references across insertion/rehashing.
template<class Key,class Value,class Hash=ankerl::unordered_dense::hash<Key>,class Equal=std::equal_to<Key>>
using FlatMap=ankerl::unordered_dense::map<Key,Value,Hash,Equal>;
}
