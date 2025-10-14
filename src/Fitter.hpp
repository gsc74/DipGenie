#pragma once
#include "Classifier.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <functional>
#include <iostream>

// ============ Toggle this to enable the BayesOpt backend ============
// Define USE_BAYESOPT if you have the BayesOpt library installed
// (https://rmcantin.github.io/bayesopt/)
//
//   g++ ... -DUSE_BAYESOPT -lbayesopt -lNLOPT
//
// ====================================================================
#ifdef USE_BAYESOPT
  #include <bayesopt/bayesopt.hpp>
  #include <bayesopt/parameters.hpp>
#endif

struct HistBin { int multiplicity; double freq; };

struct KGFitOptions {
    // histogram usage
    int   max_copy        = 20;
    int   max_x_use       = 200;
    int   smooth_win      = 7;
    bool  fit_error       = true;
    bool  fit_varw        = true;

    // bounds
    double u_lo = 1,   u_hi = 20.0; // hom mean
    double sd_lo = 0.5,  sd_hi = 2.0; // hom sd
    double varw_lo = 0.71, varw_hi = 4.0;  // het variance
    double pd_lo = 0.1,  pd_hi = 1.0; // proportion het
    double pe_lo = 0.0,  pe_hi = 0.1; // proportion error (ambigous)
    double s_lo  = 1.01,  s_hi  = 4.0; // shape (exponent)
    double zp_lo = 1.01, zp_hi = 4.0; 

    // fallback coarse search (if BayesOpt not compiled)
    int grid_u   = 7, grid_sd = 7, grid_varw = 5;
    int grid_pd  = 7, grid_pe = 5, grid_s = 5, grid_zp = 7;
    int refine_iters = 2;
};

struct KGFitResult {
    KGParams P;
    double   nll;
    int      valley_x;
    int      peak_x;
};

// ---------- tiny helpers ----------
static inline std::vector<double> moving_avg(const std::vector<double>& y, int w) {
    if (w < 1) return y;
    int n=y.size(), h=w/2;
    std::vector<double> z(n,0.0);
    for (int i=0;i<n;i++){
        int lo=std::max(0,i-h), hi=std::min(n-1,i+h);
        double s=0; int c=0;
        for (int j=lo;j<=hi;j++){ s+=y[j]; c++; }
        z[i]= s / std::max(c,1);
    }
    return z;
}
static inline int argmax(const std::vector<double>& v, int lo, int hi){
    lo = std::max(lo,0); hi = std::min(hi,(int)v.size()-1);
    int k=lo; double best=-1;
    for (int i=lo;i<=hi;i++){ if (v[i]>best){best=v[i];k=i;} }
    return k;
}
static inline double derr_old_val(int c, double s) {
    if (c <= 0) return 0.0;
    double a = std::pow((double)c, -s);
    double b = std::pow((double)(c+1), -s);
    double v = a - b;
    return (v > 0.0 ? v : 1e-300);
}
static inline std::vector<double> zeta_weights(double zp, int C){
    std::vector<double> w(C+1,0.0); double S=0.0;
    for (int k=1;k<=C;k++){ w[k] = 1.0/std::pow((double)k, zp); S += w[k]; }
    for (int k=1;k<=C;k++) w[k] /= S;
    return w;
}
static inline double normal_pdf(double x, double mu, double sd){
    double s=std::max(sd,1e-12), z=(x-mu)/s;
    static const double INV = 0.3989422804014327;
    return INV/s * std::exp(-0.5*z*z);
}

struct CompCaches {
    std::vector<double> zeta_hom, zeta_het;
    int C=0;
};
static inline void fill_caches(CompCaches& cc, const KGParams& P){
    cc.C = P.max_copy;
    cc.zeta_hom = zeta_weights(P.zp_copy,     P.max_copy);
    cc.zeta_het = zeta_weights(P.zp_copy_het, P.max_copy);
}
static inline double f_hom_x(int x, const KGParams& P, const CompCaches& cc){
    double sum=0.0;
    for (int copy=1; copy<=cc.C; ++copy) {
        double mu = copy * P.u_v;
        double sd = std::sqrt((double)copy) * P.sd_v;
        sum += cc.zeta_hom[copy] * normal_pdf(x, mu, sd);
    }
    return std::max(sum, 1e-300);
}
static inline double f_het_x(int x, const KGParams& P, const CompCaches& cc){
    double u_base = 0.5*P.u_v;
    double sd_base= 0.5*std::sqrt(std::max(P.var_w,1e-12));
    double sum=0.0;
    for (int copy=1; copy<=cc.C; ++copy) {
        double mu = copy * u_base;
        double sd = std::sqrt((double)copy) * sd_base;
        sum += cc.zeta_het[copy] * normal_pdf(x, mu, sd);
    }
    return std::max(sum, 1e-300);
}
static inline double f_err_x(int x, const KGParams& P){
    return derr_old_val(x, P.err_shape);
}

// NLL with optional Beta prior on p_d
static inline double nll_hist(const KGParams& P,
                              const std::vector<HistBin>& H,
                              int max_x_use
                             ) {
    CompCaches cc; fill_caches(cc, P);
    double nll=0.0;
    int N = std::min((int)H.size()-1, max_x_use);
    for (int x=1; x<=N; ++x){
        double y = H[x].freq; if (y<=0) continue;
        double fe  = f_err_x(x, P);
        double fhet= f_het_x(x, P, cc);
        double fhom= f_hom_x(x, P, cc);
        double mix = P.p_e * fe
                   + (1.0 - P.p_e) * ( P.p_d * fhet + (1.0 - P.p_d) * fhom );
        nll += - y * std::log(mix + 1e-300);
    }
    return nll;
}

// KmerGenie-like valley & peak for initialization/reporting
static inline void estimate_valley_peak(const std::vector<HistBin>& H, int smooth_w,
                                        int& valley_x, int& peak_x) {
    int N = (int)H.size();
    std::vector<double> y(N,0.0);
    for (int i=0;i<N;i++) y[i]=H[i].freq;
    auto ys = moving_avg(y, smooth_w);
    valley_x = 2; double vmin = ys[2];
    for (int i=2;i<std::min(N-2, 50); ++i){
        if (ys[i] < vmin){ vmin=ys[i]; valley_x=i; }
        if (i>5 && ys[i] > ys[i-1] && ys[i-1] > ys[i-2]) break;
    }
    peak_x = argmax(ys, valley_x+1, std::min(N-1, valley_x + 6*(valley_x+1)));
}

// -------- mapping between vector params and KGParams, with bounds --------
struct ParamSpec {
    // order: [u_v, sd_v, var_w, zp_copy, zp_copy_het, p_d, p_e, err_shape]
    // Some can be frozen by setting lo==hi.
    std::vector<double> lo, hi;
    std::vector<std::string> name;

    explicit ParamSpec(const KGFitOptions& opt, bool fit_error, bool fit_varw) {
        name = {"u_v","sd_v","var_w","zp_copy","zp_copy_het","p_d","p_e","err_shape"};
        lo   = {opt.u_lo, opt.sd_lo, opt.varw_lo, opt.zp_lo, opt.zp_lo, opt.pd_lo, opt.pe_lo, opt.s_lo};
        hi   = {opt.u_hi, opt.sd_hi, opt.varw_hi, opt.zp_hi, opt.zp_hi, opt.pd_hi, opt.pe_hi, opt.s_hi};
        // Note: if !fit_error or !fit_varw, caller will set lo==hi for those entries.
    }

    size_t dim() const { return lo.size(); }

    void clamp(std::vector<double>& x) const {
        for (size_t i=0;i<x.size();++i)
            x[i] = std::clamp(x[i], lo[i], hi[i]);
    }
    bool frozen(size_t i) const { return std::fabs(hi[i]-lo[i]) < 1e-12; }
};

static inline void vecToParams(const std::vector<double>& v, const ParamSpec& ps,
                               KGParams& P, int max_copy)
{
    size_t i=0;
    P.max_copy     = max_copy;
    P.u_v          = v[i++]; 
    P.sd_v         = v[i++];
    P.var_w        = v[i++];
    P.zp_copy      = v[i++];
    P.zp_copy_het  = v[i++];
    P.p_d          = v[i++];
    P.p_e          = v[i++];
    P.err_shape    = v[i++];
    // classification knobs untouched (reject_cost, amb_margin, etc.)
}

static inline std::vector<double> paramsToVec(const KGParams& P){
    return {P.u_v, P.sd_v, P.var_w, P.zp_copy, P.zp_copy_het, P.p_d, P.p_e, P.err_shape};
}

// =========================== The fitter ===========================
struct KGFitterBO {

    static KGFitResult fit(const std::vector<HistBin>& rawH, KGFitOptions opt = {}) {
        // Dense histogram 0..max_x_use
        int Nmax = 0; for (auto& b : rawH) Nmax = std::max(Nmax, b.multiplicity);
        int N = std::min(Nmax, opt.max_x_use);
        std::vector<HistBin> H(N+1); for (int x=0;x<=N;x++) H[x] = {x,0.0};
        for (auto& b : rawH) if (b.multiplicity<=N) H[b.multiplicity].freq += b.freq;

        // Valley/peak for seeds & reporting
        int valley=2, peak=10;
        estimate_valley_peak(H, opt.smooth_win, valley, peak);

        // Seeds
        auto fwhm = [&](int cx)->double{
            double pk = H[cx].freq; double half = pk/2.0; int L=cx, R=cx;
            for (int i=cx; i>=std::max(1,cx-10); --i){ if (H[i].freq<=half){ L=i; break; } }
            for (int i=cx; i<=std::min(N,cx+10); ++i){ if (H[i].freq<=half){ R=i; break; } }
            return std::max(2, R-L)/2.35;
        };
        double u_seed   = std::clamp<double>(std::max(4, peak), opt.u_lo, opt.u_hi);
        double sd_seed  = std::clamp(fwhm(peak), opt.sd_lo, opt.sd_hi);
        double varw_seed= std::clamp(2.0*sd_seed*sd_seed, opt.varw_lo, opt.varw_hi);

        // error seeds from valley
        double total=0.0, left=0.0;
        for (int x=1;x<=N;x++){ total += H[x].freq; if (x<=valley) left += H[x].freq; }
        double pe_seed = (total>0? left/total : 0.05);
        pe_seed = std::clamp(pe_seed, opt.pe_lo, opt.pe_hi);
        double s_seed  = 2.0;


        // Build initial P and bounds spec
        KGParams P0;
        P0.max_copy   = opt.max_copy;
        P0.u_v        = u_seed;
        P0.sd_v       = sd_seed;
        P0.var_w      = varw_seed;
        P0.zp_copy    = 1.5;
        P0.zp_copy_het= 1.7;
        P0.p_d        = 0.05;
        P0.p_e        = pe_seed;
        P0.err_shape  = s_seed;

        // // ------------------ ADAPTIVE BOUNDS------------------
        // // Guard against tiny means (ultra-low coverage) making formulas negative
        // auto safe_pos = [](double x, double floor=1e-6){ return std::max(x, floor); };

        // // (A) sd_v bounds
        // // probabilistic tail cap: <=10% mass below 0.5 for hom copy=1
        // double sd_tail_cap = (u_seed - 0.5) / 1.2816;           // may be <0 if u_seed<=0.5
        // sd_tail_cap = std::max(0.3, sd_tail_cap);               // keep sane floor

        // // don't allow sd much wider than observed peak (~2.5x FWHM/2.35)
        // double sd_seed_cap = 2.5 * sd_seed;

        // // final adaptive bounds for sd_v
        // opt.sd_lo = std::max(0.3, 0.5 * sd_seed);
        // opt.sd_hi = std::max(opt.sd_lo + 1e-6, std::min(sd_tail_cap, sd_seed_cap));

        // // (B) var_w bounds (het base variance)
        // // probabilistic tail cap for het copy=1: mean = 0.5*u_seed, sd_base = 0.5*sqrt(var_w)
        // // => sqrt(var_w) <= (u_seed - 1.0)/1.2816
        // double varw_tail_cap = std::pow((u_seed - 1.0) / 1.2816, 2.0);
        // varw_tail_cap = std::max(0.1, varw_tail_cap);           // floor

        // // keep within a few x of the seed estimate
        // double varw_seed_est = 2.0 * sd_seed * sd_seed; // diploid heuristic
        // double varw_seed_cap = 4.0 * varw_seed_est;

        // // final adaptive bounds for var_w
        // opt.varw_lo = std::max(0.2, 0.5 * varw_seed_est);
        // opt.varw_hi = std::max(opt.varw_lo + 1e-6, std::min(varw_tail_cap, varw_seed_cap));

        // // (optional) clamp u_v bounds to something data-aware too
        // opt.u_lo = std::max(1.2, 0.5 * u_seed);
        // opt.u_hi = std::max(opt.u_lo + 1.0, 3.0 * u_seed);

        // // Debug print (once) so you can see what box was used
        // std::cerr << "[KGFitterBO] adaptive bounds:"
        //         << " u_v=[" << opt.u_lo << "," << opt.u_hi << "]"
        //         << " sd_v=[" << opt.sd_lo << "," << opt.sd_hi << "]"
        //         << " var_w=[" << opt.varw_lo << "," << opt.varw_hi << "]\n";

        ParamSpec ps(opt, opt.fit_error, opt.fit_varw);

        // Freeze entries by setting lo==hi (honors fit_error / fit_varw)
        if (!opt.fit_varw) { ps.lo[2]=ps.hi[2]=varw_seed; }
        if (!opt.fit_error){ ps.lo[6]=ps.hi[6]=pe_seed; ps.lo[7]=ps.hi[7]=s_seed; }

        // Objective: R^d -> NLL (with bounds clamp)
        auto objective = [&](const std::vector<double>& x_unclamped)->double{
            std::vector<double> x = x_unclamped; ps.clamp(x);
            KGParams P = P0; vecToParams(x, ps, P, opt.max_copy);
            return nll_hist(P, H, opt.max_x_use);
        };

#ifdef USE_BAYESOPT
        // --- BayesOpt backend ---
        std::cout << "using BayesOpt to fit model..." << std::endl;
        size_t dim = ps.dim();

        // bounds for [lb, ub]
        std::vector<double> lb, ub;
        lb.reserve(dim); ub.reserve(dim);
        for (size_t i=0;i<dim;i++){ lb.push_back(ps.lo[i]); ub.push_back(ps.hi[i]); }

        // BayesOpt functor
        struct BOFunc : public bayesopt::ContinuousModel {
            std::function<double(const std::vector<double>&)> f;
            std::vector<double> lb, ub;

            BOFunc(size_t dim,
                const std::vector<double>& lo,
                const std::vector<double>& hi,
                std::function<double(const std::vector<double>&)> fobj,
                bayesopt::Parameters param)
            : bayesopt::ContinuousModel(dim, param), f(std::move(fobj)), lb(lo), ub(hi) {}

            // BayesOpt works in [0,1]^d; remap to [lb,ub]
            double evaluateSample(const vectord &x) override {
                std::vector<double> xv(x.size());
                for (size_t i=0;i<x.size();++i){
                    double t = std::max(0.0, std::min(1.0, x(i)));
                    xv[i] = lb[i] + t*(ub[i]-lb[i]);
                }
                return f(xv); // this is your NLL
            }
            bool checkReachability(const vectord &) override { return true; }
        };

        // parameters (note: initializer is in **global** namespace)
        bayesopt::Parameters bop = initialize_parameters_to_default();
        bop.random_seed  = 42;

        BOFunc model(dim, lb, ub, objective, bop);

        // Run optimization; BayesOpt writes the best point into `best01`
        vectord best01(dim);
        model.optimize(best01);

        // Map optimum from [0,1]^d -> [lb, ub]
        std::vector<double> xstar(dim);
        for (size_t i=0;i<dim;i++){
            double t = std::max(0.0, std::min(1.0, best01(i)));
            xstar[i] = lb[i] + t * (ub[i]-lb[i]);
        }

        // Build params at the optimum and compute its NLL
        KGParams Pbest = P0;
        vecToParams(xstar, ps, Pbest, opt.max_copy);
        double bestNLL = nll_hist(Pbest, H, opt.max_x_use);

        return KGFitResult{ Pbest, bestNLL, valley, peak };


#else
        // ========= Fallback: coarse grid =========
        // Simple linspace helper
        auto linspace = [](double lo,double hi,int k){
            std::vector<double> v; v.reserve(std::max(k,1));
            if (k<=1){ v.push_back((lo+hi)/2.0); return v; }
            for (int i=0;i<k;i++){
                double t = (double)i/(double)(k-1);
                v.push_back(lo + t*(hi-lo));
            }
            return v;
        };

        // Build grids (respect frozen params)
        auto gridOrFreeze = [&](double lo, double hi, int n)->std::vector<double>{
            if (std::fabs(hi-lo) < 1e-12) return {lo};
            return linspace(lo,hi,n);
        };
        auto U   = gridOrFreeze(ps.lo[0], ps.hi[0], opt.grid_u);
        auto SD  = gridOrFreeze(ps.lo[1], ps.hi[1], opt.grid_sd);
        auto VW  = gridOrFreeze(ps.lo[2], ps.hi[2], opt.grid_varw);
        auto ZP  = gridOrFreeze(ps.lo[3], ps.hi[3], opt.grid_zp);
        auto ZPH = gridOrFreeze(ps.lo[4], ps.hi[4], opt.grid_zp);
        auto PD  = gridOrFreeze(ps.lo[5], ps.hi[5], opt.grid_pd);
        auto PE  = gridOrFreeze(ps.lo[6], ps.hi[6], opt.grid_pe);
        auto SS  = gridOrFreeze(ps.lo[7], ps.hi[7], opt.grid_s);

        KGParams best = P0;
        double bestNLL = std::numeric_limits<double>::infinity();

        for (double u : U)
        for (double sd : SD)
        for (double vw : VW)
        for (double zp : ZP)
        for (double zph : ZPH)
        for (double pd : PD)
        for (double pe : PE)
        for (double s : SS)
        {
            KGParams P = P0;
            std::vector<double> v = {u,sd,vw,zp,zph,pd,pe,s};
            vecToParams(v, ps, P, opt.max_copy);
            double nll = nll_hist(P, H, opt.max_x_use);
            if (nll < bestNLL){ bestNLL = nll; best = P; }
        }
        return KGFitResult{ best, bestNLL, valley, peak };
#endif
    }
};
