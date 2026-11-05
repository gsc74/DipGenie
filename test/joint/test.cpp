#include "BuildModel.hpp"
#include <cassert>
#include <iostream>
#include <random>

static hic::Score independent(const hic::PresenceModel&m,const std::vector<int>&path) {
    int r=0;for(size_t b=1;b<path.size();++b)r+=(path[b]/m.H!=path[b-1]/m.H)+(path[b]%m.H!=path[b-1]%m.H);
    if(r>m.budget)return hic::impossible;
    std::vector<int> states;hic::Score value=0;
    for(const auto&t:m.targets) {
        bool a=false,b=false;
        for(const auto&[block,donors]:t.occurrences){a|=donors.test(path[block]/m.H);b|=donors.test(path[block]%m.H);}
        const auto h=t.hom?t.reward_units:t.other_reward_units;
        const auto k=t.hom?t.other_reward_units:t.reward_units;
        value+=2*m.weights.genotype*(h*(a&&b)+k*(a!=b));
        states.push_back(int(a)+2*int(b));
    }
    for(const auto&e:m.contacts){
        int a=states[e.x],b=states[e.y];
        int x=(a==1)-(a==2),y=(b==1)-(b==2);
        value+=m.weights.contact*e.units*(1+e.relation*x*y);
    }
    return value;
}

int main() {
    for(int relation:{-1,1})for(unsigned char a=0;a<4;++a)for(unsigned char b=0;b<4;++b) {
        auto swap=[](unsigned char x){return static_cast<unsigned char>(((x&1)<<1)|((x&2)>>1));};
        for(bool panel:{false,true}) {
            auto value=hic::phase_factor(a,b,relation,panel);
            assert(value>=0&&value<=2);
            assert(value==hic::phase_factor(swap(a),swap(b),relation,panel));
            assert(value+hic::phase_factor(swap(a),b,relation,panel)==2);
        }
    }
    // Confidence crossing 1/2 has no switch discontinuity in either dosage.
    for(int q=1;q<=1024;++q) {
        auto make=[](int h){return hic::PresenceTarget{h>=512,{{0,1}},h>=512?h:1024-h,h>=512?1024-h:h};};
        auto a=make(q-1),b=make(q);
        for(unsigned char s=0;s<4;++s)
            assert(hic::PresenceModel::genotype_reward(a,s)-hic::PresenceModel::genotype_reward(b,s)<=1&&hic::PresenceModel::genotype_reward(a,s)-hic::PresenceModel::genotype_reward(b,s)>=-1);
    }
    std::mt19937 rng(20261024);size_t paths=0;
    for(int trial=0;trial<96;++trial) {
        hic::PresenceModel m;m.H=3;m.blocks=3;m.budget=rng()%5;m.weights={3,2};
        for(int x=0;x<6;++x){
            hic::PresenceTarget t;t.hom=rng()%2;t.reward_units=1+rng()%10;t.other_reward_units=rng()%10;
            int first=rng()%3;
            for(int b=first;b<3;++b)if(b==first||rng()%2)t.occurrences.push_back({b,1+rng()%7});
            m.targets.push_back(t);
        }
        for(int e=0;e<12;++e){size_t x=rng()%6,y=rng()%6;if(x!=y)m.contacts.push_back({x,y,1+int(rng()%11),rng()%2?1:-1,bool(rng()%2)});}
        hic::Score optimum=hic::impossible;
        for(int a=0;a<9;++a)for(int b=0;b<9;++b)for(int c=0;c<9;++c){
            auto score=independent(m,{a,b,c});assert(score==m.score({a,b,c}));optimum=std::max(optimum,score);++paths;
        }
        for(size_t width:{size_t(0),size_t(3),size_t(16)}){
            auto one=hic::solve_presence_bidirectional(m,width,1);
            auto many=hic::solve_presence_bidirectional(m,width,4);
            assert(one.best.score==independent(m,one.best.path));
            assert(one.best.score<=optimum&&one.best.upper>=optimum&&m.universal_bound()>=optimum);
            assert(one.best.path==many.best.path&&one.best.upper==many.best.upper);
            if(!width)assert(one.best.score==optimum);
        }
        // v1 objective recovered exactly when phase coefficient is zero and
        // each target has one hard unit reward and zero alternate reward.
        m.weights={1,0};for(auto&t:m.targets){t.reward_units=1;t.other_reward_units=0;}
        auto result=hic::solve_presence_bidirectional(m,0,1);
        assert(result.best.score==independent(m,result.best.path)&&result.best.score%2==0);
    }
    std::cout<<"PASS soft confidence threshold continuity, phase symmetry/neutrality, "<<paths
             <<" exact score comparisons, joint DP bounds, parallel parity and v1 special case\n";
}
