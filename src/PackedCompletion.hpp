#pragma once
// Exact checkpointed storage: no change to the bound or candidate ordering.
#include "CompletionTransition.hpp"
#include "RowExecutor.hpp"
#include <cstring>
#include <numeric>

namespace hic {
class PackedScoreRows {
    struct Row {size_t offset;Score base,step;unsigned char bytes;};
    size_t columns_=0;
    std::vector<Row> rows_;
    std::vector<unsigned char> data_;
public:
    explicit PackedScoreRows(size_t columns=0):columns_(columns){}
    size_t columns()const{return columns_;}
    size_t rows()const{return rows_.size();}
    size_t bytes()const{return rows_.capacity()*sizeof(Row)+data_.capacity();}
    size_t payload_bytes()const{return rows_.size()*sizeof(Row)+data_.size();}
    unsigned width(size_t row)const{return rows_.at(row).bytes;}
    void append_rows(PackedScoreRows&&part) {
        if(columns_!=part.columns_)throw std::invalid_argument("Packed batch columns");
        if(part.data_.size()>SIZE_MAX-data_.size())throw std::overflow_error("Packed batch allocation");
        const size_t offset=data_.size();
        data_.insert(data_.end(),part.data_.begin(),part.data_.end());
        for(auto row:part.rows_){row.offset+=offset;rows_.push_back(row);}
        std::vector<Row>().swap(part.rows_);std::vector<unsigned char>().swap(part.data_);
    }
    void append(const std::vector<Score>&values) {
        if(!columns_||values.size()!=columns_)throw std::invalid_argument("Packed row dimensions");
        const auto limits=std::minmax_element(values.begin(),values.end());
        if(*limits.first<0||*limits.second>score_limit)throw std::invalid_argument("Packed nonnegative reward range");
        const Score base=*limits.first,range=*limits.second-base;
        Score step=0;
        for(Score value:values){Score a=step,b=value-base;while(b){Score r=a%b;a=b;b=r;}step=a;}
        if(!step)step=1;
        const unsigned __int128 largest=range/step;
        const unsigned char width=!range?0:largest<=UINT8_MAX?1:largest<=UINT16_MAX?2:largest<=UINT32_MAX?4:largest<=UINT64_MAX?8:16;
        if(width&&columns_>(SIZE_MAX-data_.size())/width)throw std::overflow_error("Packed row allocation");
        const size_t offset=data_.size();
        data_.resize(offset+columns_*width);
        for(size_t z=0;width&&z<columns_;++z) {
            const unsigned __int128 code=(values[z]-base)/step;
            auto* p=data_.data()+offset+z*width;
            switch(width) {
                case 1: {uint8_t v=uint8_t(code);std::memcpy(p,&v,1);break;}
                case 2: {uint16_t v=uint16_t(code);std::memcpy(p,&v,2);break;}
                case 4: {uint32_t v=uint32_t(code);std::memcpy(p,&v,4);break;}
                case 8: {uint64_t v=uint64_t(code);std::memcpy(p,&v,8);break;}
                case 16: std::memcpy(p,&code,16);break;
            }
        }
        rows_.push_back({offset,base,step,width});
    }
    void decode(size_t row,std::vector<Score>&output)const {
        const auto&r=rows_.at(row);output.resize(columns_);
        if(!r.bytes){std::fill(output.begin(),output.end(),r.base);return;}
        const auto* p=data_.data()+r.offset;
        for(size_t z=0;z<columns_;++z) {
            unsigned __int128 code=0;
            switch(r.bytes) {
                case 1: {uint8_t v;std::memcpy(&v,p+z,1);code=v;break;}
                case 2: {uint16_t v;std::memcpy(&v,p+2*z,2);code=v;break;}
                case 4: {uint32_t v;std::memcpy(&v,p+4*z,4);code=v;break;}
                case 8: {uint64_t v;std::memcpy(&v,p+8*z,8);code=v;break;}
                case 16: std::memcpy(&code,p+16*z,16);break;
            }
            output[z]=r.base+Score(code)*r.step;
        }
    }
};

class PackedCompletion {
    using Entry=CompletionTransition::Entry;
    bool free_swaps_=false;
    int H_=0,blocks_=0,stride_=0,cache_begin_=-1,cache_end_=-1;
    size_t domain_=0;
    PackedScoreRows local_;
    std::vector<Score> checkpoints_,cache_,other_;
public:
    Score penalty=0,free_upper=0;
    int relaxed_switches=0;
    size_t reconstructed_blocks=0;
    PackedCompletion()=default;
    PackedCompletion(int H,int blocks,int budget,PackedScoreRows local,std::vector<Score> other,int stride=0,int threads=1,bool free_swaps=false):
        free_swaps_(free_swaps),H_(H),blocks_(blocks),domain_(size_t(checked_pair_domain(size_t(H)))),
        local_(std::move(local)),other_(std::move(other)) {
        if(blocks<1||blocks>INT32_MAX/2||budget<0||stride<0||threads<1||
           local_.columns()!=domain_||local_.rows()!=size_t(blocks)||other_.size()!=size_t(blocks+1))
            throw std::invalid_argument("Packed completion dimensions");
        stride_=stride?stride:std::max(1,int(std::sqrt(double(blocks))));
        const size_t checkpoints=(size_t(blocks)+size_t(stride_)-1)/size_t(stride_);
        if(checkpoints>SIZE_MAX/domain_)throw std::overflow_error("Completion checkpoint dimensions");
        Score range=0,total=0;
        std::vector<Score> row;
        for(int b=0;b<blocks;++b) {
            local_.decode(size_t(b),row);
            const auto limits=std::minmax_element(row.begin(),row.end());
            if(total>score_limit-*limits.second)throw std::overflow_error("Completion reward range");
            total+=*limits.second;range+=*limits.second-*limits.first;
        }
        auto scan=[&](Score lambda,bool save) {
            std::vector<Score> row;
            std::vector<Entry> next(domain_,{0,0});
            if(save)checkpoints_.assign(checkpoints*domain_,0);
            for(int b=blocks-1;b>=0;--b) {
                if(save&&(b==blocks-1||(b+1)%stride_==0)) {
                    const size_t offset=size_t(b/stride_)*domain_;
                    for(size_t z=0;z<domain_;++z)checkpoints_[offset+z]=next[z].score;
                }
                local_.decode(size_t(b),row);
                for(size_t z=0;z<domain_;++z)next[z].score+=row[z];
                if(b)CompletionTransition::transition(next,H,lambda,free_swaps_);
            }
            Entry best=next[0];
            for(auto x:next)if(CompletionTransition::better(x,best))best=x;
            return best;
        };
        const auto zero=scan(0,false);
        __int128 dual=zero.score;
        if(zero.switches>budget&&range>0) {
            Score low=0,high=range;
            auto consider=[&](Score lambda) {
                const auto value=scan(lambda,false);
                const __int128 bound=__int128(value.score)+__int128(lambda)*budget;
                if(bound<dual){dual=bound;penalty=lambda;}
                return value.switches;
            };
            (void)consider(high);
            // Evaluate a few levels of the SAME binary search speculatively.
            // Only probes on its chosen branch update the bound. Thus parallel
            // scheduling cannot change the penalty, ties, or candidate order.
            int depth=1;
            while(depth<4&&((size_t(1)<<(depth+1))-1)<=size_t(threads))++depth;
            RowExecutor pool(std::min(threads,15));
            for(int iteration=0;iteration<16&&high-low>1;) {
                struct Probe {Score low,high,mid;Entry result{0,0};size_t lower=SIZE_MAX,upper=SIZE_MAX;};
                std::vector<Probe> probes;
                std::function<size_t(Score,Score,int)> make=[&](Score l,Score h,int levels) {
                    if(!levels||h-l<=1)return SIZE_MAX;
                    const size_t i=probes.size();const Score mid=l+(h-l)/2;
                    probes.push_back({l,h,mid,{0,0},SIZE_MAX,SIZE_MAX});
                    const size_t lower=make(l,mid,levels-1),upper=make(mid,h,levels-1);
                    probes[i].lower=lower;probes[i].upper=upper;return i;
                };
                make(low,high,std::min(depth,16-iteration));
                pool.for_each(probes.size(),[&](size_t i){probes[i].result=scan(probes[i].mid,false);});
                for(size_t i=0;i!=SIZE_MAX;) {
                    const auto&p=probes[i];const auto value=p.result;
                    const __int128 bound=__int128(value.score)+__int128(p.mid)*budget;
                    if(bound<dual){dual=bound;penalty=p.mid;}
                    ++iteration;
                    if(value.switches>budget){low=p.mid;i=p.upper;}
                    else {high=p.mid;i=p.lower;}
                }
            }
        }
        const auto selected=scan(penalty,true);relaxed_switches=selected.switches;
        if(dual>score_limit)throw std::overflow_error("Completion bound");
        free_upper=Score(dual);
    }
    // Called once, serially, before a block's queries. All subsequent queries
    // are read-only and may run in parallel. Never reconstruct inside a worker.
    void prepare(int block) {
        if(block<0||block>=blocks_)throw std::invalid_argument("Completion block");
        if(block>=cache_begin_&&block<cache_end_)return;
        cache_begin_=block/stride_*stride_;
        cache_end_=int(std::min(size_t(blocks_),size_t(cache_begin_)+size_t(stride_)));
        const size_t rows=size_t(cache_end_-cache_begin_);
        if(rows>SIZE_MAX/domain_)throw std::overflow_error("Completion tile dimensions");
        cache_.resize(rows*domain_);
        std::vector<Entry> next(domain_);
        const size_t offset=size_t(block/stride_)*domain_;
        for(size_t z=0;z<domain_;++z)next[z]={checkpoints_[offset+z],0};
        std::vector<Score> local;
        for(int b=cache_end_-1;b>=cache_begin_;--b) {
            for(size_t z=0;z<domain_;++z)cache_[size_t(b-cache_begin_)*domain_+z]=next[z].score;
            if(b==cache_begin_)break;
            local_.decode(size_t(b),local);
            for(size_t z=0;z<domain_;++z)next[z].score+=local[z];
            CompletionTransition::transition(next,H_,penalty,free_swaps_);
        }
        reconstructed_blocks+=rows;
    }
    Score remaining(int block,int pair,int budget_left,Score independent)const {
        if(block<0||block>=blocks_||pair<0||size_t(pair)>=domain_||budget_left<0||independent<0)
            throw std::invalid_argument("Completion query");
        if(block+1==blocks_)return 0;
        if(block<cache_begin_||block>=cache_end_)throw std::logic_error("Completion tile not prepared");
        const __int128 bound=__int128(cache_[size_t(block-cache_begin_)*domain_+size_t(pair)])+
                             __int128(penalty)*budget_left+other_[block+1];
        return Score(std::min(bound,__int128(independent)));
    }
    // All histories with the same remaining resource share this row. Check
    // and multiply once, then decode the contiguous destination bounds.
    void remaining_row(int block,int budget_left,Score independent,Score*output)const {
        if(block<0||block>=blocks_||budget_left<0||independent<0||!output)
            throw std::invalid_argument("Completion row query");
        if(block+1==blocks_){std::fill_n(output,domain_,Score(0));return;}
        if(block<cache_begin_||block>=cache_end_)throw std::logic_error("Completion tile not prepared");
        const Score extra=penalty*budget_left+other_[block+1];
        const auto*row=cache_.data()+size_t(block-cache_begin_)*domain_;
        for(size_t z=0;z<domain_;++z)output[z]=std::min(row[z]+extra,independent);
    }
    size_t storage_bytes()const {
        return local_.bytes()+(checkpoints_.capacity()+cache_.capacity()+other_.capacity())*sizeof(Score);
    }
    size_t packed_local_payload_bytes()const{return local_.payload_bytes();}
};
}
