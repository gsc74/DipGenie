#pragma once
// Exact reward representation, normalization and a posteriori certificate.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

// GCC's exact signed 128-bit arithmetic avoids chromosome-scale overflow
// from clearing the (fixed) contact denominator. This changes neither the
// input quantization, objective, candidate ordering nor certificate.
inline std::ostream& operator<<(std::ostream&out,__int128 value) {
    if(value<0){out<<'-';value=-(value+1);unsigned __int128 n=static_cast<unsigned __int128>(value)+1;
        char buf[40];size_t i=0;do{buf[i++]=char('0'+n%10);n/=10;}while(n);while(i)out<<buf[--i];return out;}
    unsigned __int128 n=value;char buf[40];size_t i=0;do{buf[i++]=char('0'+n%10);n/=10;}while(n);while(i)out<<buf[--i];return out;
}

namespace hic {
using Score=__int128; // Fixed-point INPUT weights; all DP/bound sums exact.
constexpr Score score_limit=std::numeric_limits<Score>::max()/4;
inline std::string score_text(Score s){std::ostringstream out;out<<s;return out.str();}
inline Score score_product(Score a,Score b){if(a<0||b<0||(b&&a>score_limit/b))throw std::overflow_error("score product");return a*b;}
inline Score score_sum(Score a,Score b){if(a<0||b<0||a>score_limit-b)throw std::overflow_error("score sum");return a+b;}
constexpr Score impossible=std::numeric_limits<Score>::min()/4;
struct Result {
    std::vector<int> path;
    Score score=0, upper=0;
    size_t generated=0, discarded=0, peak=0, max_frontier=0;
    Score completion_penalty=0;
    int completion_relaxed_switches=0;
    double ratio() const {
        if(path.empty() || (score==0 && upper>0)) return std::numeric_limits<double>::quiet_NaN();
        return score==0?1:double(upper)/double(score);
    }
    std::string ratio_text() const {
        if(path.empty() || (score==0 && upper>0)) return "NULL";
        if(score<0 || upper<score) throw std::logic_error("invalid reward certificate");
        if(score==0) return "1.000000";
        // Round UP so even the displayed decimal is a valid bound. Retain
        // the exact integer (upper,score) certificate in the run record.
        constexpr Score scale=1000000;
        using boost::multiprecision::cpp_int;
        const cpp_int numerator(score_text(upper)),denominator(score_text(score));
        const cpp_int exact=(numerator*scale+denominator-1)/denominator;
        if(exact>cpp_int(score_text(score_limit)))throw std::overflow_error("displayed certificate");
        Score scaled=exact.convert_to<Score>();
        std::ostringstream out;
        out<<score_text(scaled/scale)<<'.'<<std::setfill('0')<<std::setw(6)<<score_text(scaled%scale);
        return out.str();
    }
};
// Clear the contact-average denominator exactly; lambda=a/b is fixed INPUT.
// Contact units must have been allocated once per physical molecule, with
// their total equal to W. No candidate-dependent renormalization is allowed.
struct NormalizedWeights {
    Score genotype, contact;
};
inline NormalizedWeights normalized_weights(Score W, Score het, Score a=1, Score b=1) {
    if(W<0 || het<0 || a<0 || b<=0) throw std::invalid_argument("normalization inputs");
    if(W==0 || a==0 || het==0) return {1,0};
    Score g=score_product(b,W), c=score_product(a,het);
    auto gcd128=[](__int128 x,__int128 y) {while(y) {auto t=x%y;x=y;y=t;}return x;};
    auto d=gcd128(g,c); g/=d;c/=d;
    if(g>std::numeric_limits<Score>::max()/4 || c>std::numeric_limits<Score>::max()/4)
        throw std::overflow_error("normalization coefficient storage");
    return {Score(g),Score(c)};
}
inline Score weighted_score(Score g,Score c,NormalizedWeights w) {
    if(g<0 || c<0 || w.genotype<0 || w.contact<0) throw std::invalid_argument("score inputs");
    Score v=score_sum(score_product(g,w.genotype),score_product(c,w.contact));
    if(v>std::numeric_limits<Score>::max()/4) throw std::overflow_error("normalized score storage");
    return Score(v);
}

}
