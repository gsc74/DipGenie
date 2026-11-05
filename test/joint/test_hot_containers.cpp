#include "FlatHash.hpp"
#include "IntegerSort.hpp"
#include <cassert>
#include <iostream>
#include <random>
#include <unordered_map>
using namespace hic;
struct Collision {size_t operator()(uint64_t)const{return 0;}};
template<class Hash> void table_test() {
    FlatMap<uint64_t,uint64_t,Hash> table;
    std::unordered_map<uint64_t,uint64_t> expected;
    std::mt19937_64 random(481);
    for(int i=0;i<12000;++i) {
        const uint64_t key=i%11==0?UINT64_MAX:i%7==0?0:random()%1024;
        if(i%5==0){assert(table.erase(key)==expected.erase(key));}
        else {auto [a,new_a]=table.try_emplace(key,0);auto [b,new_b]=expected.try_emplace(key,0);
            assert(new_a==new_b);a->second+=i;b->second+=i;}
        if(i%101==0)table.reserve(table.size()+71);
        assert(table.size()==expected.size());
    }
    for(auto[key,value]:expected)assert(table.at(key)==value);
    const auto&frozen=table;
    for(uint64_t key=0;key<2048;++key)assert((frozen.find(key)==frozen.end())==(expected.find(key)==expected.end()));
    auto copied=table;assert(copied==table);table.clear();assert(table.empty());
    table.swap(copied);for(auto[key,value]:expected)assert(table.at(key)==value);
}
int main() {
    table_test<ankerl::unordered_dense::hash<uint64_t>>();table_test<Collision>();
    std::mt19937_64 random(258);
    for(size_t n:std::vector<size_t>{0,1,2,63,64,65,4095,4096,4097,65537})for(int shape=0;shape<7;++shape) {
        std::vector<uint64_t> a(n);
        for(size_t i=0;i<n;++i)a[i]=shape==0?random():shape==1?i:shape==2?n-i:shape==3?7:
            shape==4?random()%17:shape==5?(UINT64_MAX-(i%9)):(random()&255)<<48;
        auto expected=a;std::sort(expected.begin(),expected.end());sort_words(a);assert(a==expected);
    }
    sort_words(nullptr,0);std::cout<<"hot container collision/growth/clear and integer-sort parity passed\n";
}
