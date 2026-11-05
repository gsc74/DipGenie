#include "NativeFactors.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace hic;
namespace nf=hic::native_factors;
static MoleculeMarkers observation(std::string id,uint32_t lib,MoleculeKind source,MoleculeGeometry geometry,int64_t distance) {
    MoleculeMarkers f;f.id=std::move(id);f.library=lib;f.kind=source;f.geometry=geometry;
    f.reference_distance=distance;f.localization_qc=1;f.ends[0]={0,1,1};f.ends[1]={2,3,4};return f;
}
static std::string text(const nf::Snapshot&s){std::ostringstream out;nf::write_snapshot(out,s);return out.str();}
int main() {
    PreparedModel p;auto&m=p.model;m.H=2;m.blocks=3;m.budget=4;m.free_swaps=true;
    SequenceIndex index;index.targets.resize(5);
    p.target_ids={0,1,1,2,3,4};p.linkage_remap={0,1,3,SIZE_MAX,5};
    // Word 1 has two copy levels; word 3 has no dosage evidence. Neither the
    // second level nor the zero-count word may become a molecule endpoint.
    for(size_t i=0;i<p.target_ids.size();++i) {
        const Score het[]={1024,512,1024,256,1024,128};
        m.targets.push_back({false,{{int(i%3),DonorMask(i==3?2:1)}},het[i],0,0});
        p.het_reward_units+=het[i];
    }
    m.contacts={{0,3,77,-1,true},{0,3,13,1,true}};
    RawEvidence raw;
    raw.markers.push_back(observation("a",0,MoleculeKind::ShortRead,MoleculeGeometry::Known,299));
    raw.markers.push_back(observation("b",1,MoleculeKind::HiC,MoleculeGeometry::Known,300));
    raw.markers.push_back(observation("c",1,MoleculeKind::HiC,MoleculeGeometry::Unknown,unknown_contact_distance));
    raw.markers.push_back(observation("d",2,MoleculeKind::HiC,MoleculeGeometry::Trans,unknown_contact_distance));
    raw.markers.back().localization_qc=2;
    for(uint64_t edge:{300,1000,3000,10000,30000,100000,300000,1000000,3000000})
        assert(nf::distance_bin(edge)==nf::distance_bin(edge-1)+1);
    assert(nf::classify(raw.markers[0]).source!=nf::classify(raw.markers[1]).source);
    auto unknown=raw.markers[0];unknown.geometry=MoleculeGeometry::Unknown;unknown.reference_distance=unknown_contact_distance;
    auto zero=raw.markers[0];zero.reference_distance=0;
    assert(!(nf::classify(unknown)==nf::classify(zero)));
    for(bool stage2:{false,true}) {
        const auto canonical=nf::prepare(p,index,raw,stage2,false,1);
        assert(canonical.molecules==4&&canonical.raw_units==4096&&canonical.read_units<=4096&&canonical.classes.size()==5);
        for(int threads:{2,4,8})assert(text(nf::prepare(p,index,raw,stage2,false,threads))==text(canonical));
        // Independent small arithmetic oracle: unique unordered pairs, FNV
        // remainder order, allocation, then ONE outer heterozygosity discount.
        std::map<std::tuple<nf::Class,size_t,size_t>,Score> expected;
        for(const auto&f:raw.markers) {
            std::set<std::pair<size_t,size_t>> pairs;
            for(size_t a:f.ends[0])for(size_t b:f.ends[1]) {
                auto x=p.linkage_remap[a],y=p.linkage_remap[b];if(x==SIZE_MAX||y==SIZE_MAX||x==y)continue;
                pairs.emplace(std::min(x,y),std::max(x,y));
            }
            std::vector<std::pair<size_t,size_t>> list(pairs.begin(),pairs.end());
            std::vector<Score> units(list.size());uint64_t h=1469598103934665603ULL;
            for(unsigned char c:f.id){h^=c;h*=1099511628211ULL;}
            Score den=0;for(auto[x,y]:list)den+=stage2?m.targets[x].reward_units*m.targets[y].reward_units:1;
            Score used=0;
            for(size_t i=0;i<list.size();++i){auto[x,y]=list[i];units[i]=1024*(stage2?m.targets[x].reward_units*m.targets[y].reward_units:1)/den;used+=units[i];}
            for(size_t j=0;used<1024;++j){++units[(h%list.size()+j)%list.size()];++used;}
            for(size_t i=0;i<list.size();++i){auto[x,y]=list[i];expected[{nf::classify(f),x,y}]+=units[i];}
        }
        size_t checked=0;
        for(const auto&g:canonical.classes)for(const auto&e:g.factors)if(!e.panel) {
            assert(e.x!=2&&e.y!=2&&e.x!=4&&e.y!=4);
            const Score want=expected.at({g.identity,e.x,e.y})*m.targets[e.x].reward_units*m.targets[e.y].reward_units/(1024*1024);
            assert(e.units==want);++checked;
        }
        assert(checked==expected.size());
        auto model=m;const auto total=nf::apply(canonical,model);
        assert(total==canonical.read_units+90);
        size_t signs=0;for(const auto&e:model.contacts)if(e.panel&&e.x==0&&e.y==3)++signs;
        assert(signs==2); // opposite panel signs must remain separate
    }
    // Synthetic calibration with known orientations; same endpoints but
    // distinct source/library/status/QC classes must not share confidence.
    nf::Snapshot snap;nf::Class sr{nf::Source::SR,0,MoleculeGeometry::Known,0,1};
    nf::Class hc{nf::Source::HiC,1,MoleculeGeometry::Known,0,1};
    nf::Class unk{nf::Source::HiC,1,MoleculeGeometry::Unknown,0,1};
    nf::Class panel;
    snap.classes={{sr,{{0,1,1024,1,false}}},{hc,{{0,3,1024,1,false}}},
                  {unk,{}},{panel,{{0,1,3072,1,true},{0,3,1024,1,true}}}};
    const std::vector<int> path(3,1);
    const auto est=nf::estimate(snap,m,path);
    assert(est.classes.at(sr).agree==1024&&est.classes.at(sr).confidence==1024);
    assert(est.classes.at(hc).oppose==1024&&est.classes.at(hc).confidence==0);
    assert(est.classes.at(unk).epsilon==.5L&&est.classes.at(unk).confidence==0);
    assert(std::fabs(est.classes.at(sr).epsilon-.25L)<1e-15L);
    const auto pc=uint32_t(std::llround(1024*std::log(7.L/3.L)/std::log(3.L)));
    assert(est.classes.at(panel).confidence==pc);
    const auto old=m.score(path);nf::apply(snap,m,&est);nf::normalize(p);
    assert(old!=m.score(path));
    // Exhaust all path pairs: independently sum literal dosage/factors and
    // check the certificate under the new normalization, not the old score.
    Score exact_best=impossible;
    for(int a=0;a<4;++a)for(int b=0;b<4;++b)for(int c=0;c<4;++c) {
        std::vector<int> q{a,b,c};Score independent=0;
        std::vector<int> orientation;
        for(const auto&t:m.targets) {
            unsigned presence=0;for(auto[block,mask]:t.occurrences) {
                presence|=(mask.test(q[block]/2)?1:0)|(mask.test(q[block]%2)?2:0);
            }
            const int copies=bool(presence&1)+bool(presence&2);
            const Score v=copies==0?t.absent_reward_units:copies==1?t.reward_units:t.other_reward_units;
            independent+=2*m.weights.genotype*v;
            orientation.push_back(int(bool(presence&1))-int(bool(presence&2)));
        }
        for(const auto&e:m.contacts)independent+=m.weights.contact*e.units*(1+e.relation*orientation[e.x]*orientation[e.y]);
        assert(m.score(q)==independent&&independent<=m.universal_bound());
        exact_best=std::max(exact_best,independent);
    }
    // New admitted observations change the fixed factor array, NOT the DP
    // state. Check its exact recurrence and pruned certificate independently
    // against the exhaustive joint-pair optimum after calibration/normalization.
    std::vector<int> canonical_exact,canonical_pruned;
    for(int threads:{1,2,4,8}) {
        const auto exact=solve_presence(m,0,threads);
        assert(exact.score==exact_best&&m.score(exact.path)==exact_best);
        const auto pruned=solve_presence(m,2,threads,{},true);
        assert(!pruned.path.empty()&&m.score(pruned.path)==pruned.score);
        assert(pruned.score<=exact_best&&exact_best<=pruned.upper&&pruned.upper<=m.universal_bound());
        assert(Result({pruned.path,pruned.score,m.universal_bound()}).ratio()>=1);
        if(threads==1){canonical_exact=exact.path;canonical_pruned=pruned.path;}
        else{assert(exact.path==canonical_exact&&pruned.path==canonical_pruned);}
    }
    nf::Snapshot empty;assert(nf::estimate(empty,m,path).maximum_llr==0);
    nf::apply(empty,m);nf::normalize(p);assert(m.weights.contact==0);
    p.het_reward_units=0;nf::apply(snap,m);nf::normalize(p);assert(m.weights.contact==0);
    auto bad=raw.markers[0];bad.reference_distance=-1;bool rejected=false;
    try{nf::classify(bad);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    bad=raw.markers[0];bad.geometry=MoleculeGeometry::Legacy;rejected=false;
    try{nf::classify(bad);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
    // Overflow is rejected before pooled scores can silently wrap.
    auto oversized=snap;oversized.classes[0].factors[0].units=score_limit/nf::scale+1;
    rejected=false;try{nf::apply(oversized,m,&est);}catch(const std::overflow_error&){rejected=true;}assert(rejected);
    // A half-confidence observation of one allocation unit is still
    // evidence, not zero. Retaining the shared denominator is exact for
    // all candidate paths and their factor-wise bounds.
    {
        nf::Snapshot tiny;tiny.classes={{sr,{{0,1,1,1,false}}},{hc,{{0,3,3,1,false}}}};
        nf::Estimate fixed;fixed.classes[sr].confidence=512;fixed.classes[hc].confidence=256;
        auto a=p,b=p;a.het_reward_units=b.het_reward_units=100;
        nf::apply(tiny,a.model,&fixed);nf::normalize(a);
        assert(a.model.contacts.size()==2&&a.model.contacts[0].units==512&&a.model.contacts[1].units==768);
        // Same rational weights in reduced integer units 2 and 3.
        b.model.contacts={{0,1,2,1,false},{0,3,3,1,false}};nf::normalize(b);
        for(int x=0;x<4;++x)for(int y=0;y<4;++y)for(int z=0;z<4;++z) {
            const std::vector<int> path{x,y,z};
            const auto sa=a.model.score(path),sb=b.model.score(path);
            assert(sa*b.model.weights.genotype==sb*a.model.weights.genotype);
        }
        assert(a.model.universal_bound()*b.model.weights.genotype==b.model.universal_bound()*a.model.weights.genotype);
    }
    // No heterozygous allocation support means no molecule factor, not 0/0.
    auto zero_het=p;for(auto&t:zero_het.model.targets){t.reward_units=0;t.other_reward_units=0;}
    const auto omitted=nf::prepare(zero_het,index,raw,true,false,4);
    assert(omitted.molecules==0&&omitted.raw_units==0&&omitted.read_units==0);
    std::cout<<"PASS native factors: source classes, strict bins, fixed eligibility, allocation oracle, mass, frozen calibration, score/bound replay and thread parity\n";
}
