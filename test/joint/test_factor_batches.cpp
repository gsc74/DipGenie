#include "NativeFactors.hpp"
#include <cassert>
#include <iostream>
#include <random>
namespace nf=hic::native_factors;
using namespace hic;
int main() {
    PreparedModel p;p.model.H=2;p.model.blocks=3;p.model.budget=4;
    SequenceIndex index;index.targets.resize(8);
    for(size_t x=0;x<8;++x) {
        p.target_ids.push_back(x);p.linkage_remap.push_back(x);
        p.model.targets.push_back({false,{{int(x%3),DonorMask(1)}},Score(128+64*x),Score(x%3),0});
    }
    RawEvidence raw;std::mt19937_64 random(3915);
    // More than a full 64-worker batch plus a partial final batch. Repeated
    // factor keys cross both task and batch boundaries; IDs stay distinct.
    for(size_t i=0;i<16403;++i) {
        MoleculeMarkers m;m.id=std::to_string(i);m.library=i%3;
        m.kind=i%2?MoleculeKind::HiC:MoleculeKind::ShortRead;
        m.geometry=MoleculeGeometry::Known;m.reference_distance=100+i%400;
        m.localization_qc=1+i%2;
        for(int end=0;end<2;++end)for(int j=0;j<4;++j)m.ends[end].push_back(random()%8);
        raw.markers.push_back(std::move(m));
    }
    for(bool stage2:{false,true}) {
        auto canonical=nf::prepare(p,index,raw,stage2,false,1);
        std::ostringstream expected;nf::write_snapshot(expected,canonical);
        for(int threads:{2,8,64}) {
            const auto actual=nf::prepare(p,index,raw,stage2,false,threads);
            std::ostringstream out;nf::write_snapshot(out,actual);
            assert(out.str()==expected.str());
        }
    }
    std::cout<<"native factor repeated-key, multi-batch and partial-batch parity passed\n";
}
