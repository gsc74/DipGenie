#pragma once
#include "PhaseFactor.hpp"
#include <tuple>

namespace hic {
// Two unresolved contacts sharing a future target cannot choose that target's
// orientation independently. Once their earlier endpoints have signs a,b,
// their joint maximum is w1+w2+|r1*w1*a+r2*w2*b|, not 2*(w1+w2).
// Disjoint pairs of original edges give an admissible completion tightening.
// Subtract the deficit when both earlier endpoints are known, and restore it
// when the common future endpoint is scored. Thus every complete path retains
// EXACTLY its original objective value; no evidence or tuning weight is added.
struct JointContactEvent {
    size_t x,y;
    Score a,b;
    int ra,rb,sign; // sign=-1 opens the pending deficit; +1 closes it
    Score value(unsigned char px,unsigned char py)const {
        Score field=a*ra*phase_orientation(px)+b*rb*phase_orientation(py);
        if(field<0)field=-field;
        return sign*(a+b-field);
    }
};
template<class Contacts>
std::vector<std::vector<JointContactEvent>> joint_contact_events(
        const Contacts&contacts,const std::vector<int>&first,
        const std::vector<int>&last,int blocks,Score weight,bool enabled) {
    std::vector<std::vector<JointContactEvent>> events(blocks);
    if(!enabled||weight==0)return events;
    struct Link{size_t future,past;Score units;int relation;};
    std::vector<Link> links;links.reserve(contacts.size());
    for(const auto&e:contacts) {
        if(first[e.x]!=last[e.x]||first[e.y]!=last[e.y]||first[e.x]==first[e.y])continue;
        const bool reverse=first[e.x]>first[e.y];
        links.push_back({reverse?e.x:e.y,reverse?e.y:e.x,e.units,e.relation});
    }
    std::sort(links.begin(),links.end(),[&](const auto&a,const auto&b){
        return std::make_tuple(a.future,first[a.past],a.past,a.relation,a.units)<
               std::make_tuple(b.future,first[b.past],b.past,b.relation,b.units);
    });
    for(size_t i=0;i<links.size();) {
        size_t end=i+1;while(end<links.size()&&links[end].future==links[i].future)++end;
        for(size_t j=i;j+1<end;j+=2) {
            const auto&a=links[j],&b=links[j+1];
            const int open=std::max(first[a.past],first[b.past]),close=first[a.future];
            JointContactEvent e{a.past,b.past,score_product(a.units,weight),score_product(b.units,weight),a.relation,b.relation,-1};
            events[open].push_back(e);e.sign=1;events[close].push_back(e);
        }
        i=end;
    }
    return events;
}
}
