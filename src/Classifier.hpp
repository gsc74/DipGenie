#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <cstdint>

// Standalone C++ port of KmerGenie’s diploid “blue/green” component logic
// from model.r + model-diploid.r (with shape.v = 0, i.e., plain normal kernels).
// - Zeta/Zipf prior over copy number: dzeta(k; zp) ∝ 1 / k^zp (normalized over 1..max_copy)
// - Kernel at abundance x for a given copy: Normal(mean = copy * u_v, sd = sqrt(copy) * sd_v)
// - Homozygous good:   (zp.copy,     u_v,              sd_v)
// - Heterozygous good: (zp.copy.het, 0.5 * u_v,        0.5 * sqrt(var_w))
// - Mixture weight p_d: probs_good = p_d * probs_het_good + (1-p_d) * probs_hom_good

struct KGParams {
    // model-diploid.r parameters we actually need
    double zp_copy      = 1.3;   // zp.copy      (hom heavy-tail)
    double zp_copy_het  = 1.3;   // zp.copy.het  (het heavy-tail)
    double u_v          = 4.0;   // abundance at copy=1 (hom base)
    double sd_v         = 1.2;   // sd at copy=1 (hom base)
    double var_w        = 2.0;   // variance at copy=0.5 (het base): sd_het_base = 0.5*sqrt(var_w)
    double p_d          = 0.5;   // heterozygous proportion among good
    int    max_copy     = 5;     // truncate the zeta mixture at this copy number 
    // classification option
    double amb_margin   = 0.05;  // min posterior gap to avoid AMBIGUOUS
    double p_e = 0.01;       // small but non-zero error prior
    double err_shape = 2.0;  // s>1; 1.8–3.0 typical
    bool   treat_error_as_ambiguous = true; // label ERROR as AMB (or keep ERROR explicitly)
    double reject_cost = 0.4; // need at least 1 - reject cost probability to not mark as AMB
};

struct KGPosterior {
    double p_het, p_hom;  // normalized posteriors for a single abundance x
    enum Label { HET, HOM, AMBIGUOUS, ERR } label;
};

class KmerGenieDiploidLike {
public:
    explicit KmerGenieDiploidLike(const KGParams& P) : P_(P) {
        if (P_.max_copy < 1) throw std::invalid_argument("max_copy must be >= 1");
        if (P_.u_v <= 0 || P_.sd_v <= 0) throw std::invalid_argument("u_v, sd_v must be > 0");
        // precompute normalized zeta/zipf weights for both components
        zeta_hom_ = zetaWeights_(P_.zp_copy,     P_.max_copy);
        zeta_het_ = zetaWeights_(P_.zp_copy_het, P_.max_copy);
        sd_het_base_ = 0.5 * std::sqrt(std::max(P_.var_w, 1e-12)); // model-diploid.r
    }

    // Blue curve (het): probs.het.good over x = 1..max_x (unnormalized to histogram; same scale as model curves)
    std::vector<double> probsHetGood(int max_x) const {
        return pgood_(max_x, /*u_base=*/0.5*P_.u_v, /*sd_base=*/sd_het_base_, zeta_het_);
    }
    // Green curve (hom): probs.hom.good over x = 1..max_x
    std::vector<double> probsHomGood(int max_x) const {
        return pgood_(max_x, /*u_base=*/P_.u_v,     /*sd_base=*/P_.sd_v,       zeta_hom_);
    }

    // Posterior classification for a single multiplicity x (>=1), using p_d as prior weight
    KGPosterior classify(int x) const {
        double fe  = derr_old_val(x, P_.err_shape);                  // error pmf
        double fhet= valHet_(x);                                     // blue curve height
        double fhom= valHom_(x);                                     // green curve height

        double a = P_.p_e * fe;
        double b = (1.0 - P_.p_e) * P_.p_d * fhet;
        double c = (1.0 - P_.p_e) * (1.0 - P_.p_d) * fhom;

        double Z = std::max(a + b + c, 1e-300);
        double perr = a / Z, phet = b / Z, phom = c / Z;

        // simified, always hom or het
        KGPosterior::Label L;
        if (x == 1 || phet >= phom) {
            L = KGPosterior::HET;
        } else {
            L = KGPosterior::HOM;
        }
        //}
        return {/*p_het=*/phet, /*p_hom=*/phom, L};
    }

    // Partition (id,count) pairs into het/hom/amb
    template <class K>
    struct Partition { std::vector<K> het, hom, amb; };

    template <class K>
    Partition<K> partition(const std::vector<std::pair<K,uint32_t>>& kmers) const {
        Partition<K> P;
        P.het.reserve(kmers.size()/2);
        P.hom.reserve(kmers.size()/2);
        for (auto& kv : kmers) {
            auto post = classify((int)kv.second);
            switch (post.label) {
                case KGPosterior::HET: P.het.push_back(kv.first); break;
                case KGPosterior::HOM: P.hom.push_back(kv.first); break;
                default:               P.amb.push_back(kv.first); break;
            }
        }
        return P;
    }

    const KGParams& params() const { return P_; }

private:
    KGParams P_;
    std::vector<double> zeta_hom_, zeta_het_;
    double sd_het_base_{1.0};

    static inline double normal_pdf_(double x, double mu, double sd) {
        double s = std::max(sd, 1e-12);
        double z = (x - mu) / s;
        static const double INV_SQRT_2PI = 0.3989422804014327;
        return INV_SQRT_2PI / s * std::exp(-0.5 * z * z);
    }

    static inline double derr_old_val(int c, double s) {
        // p_err(c) = 1/c^s - 1/(c+1)^s
        if (c <= 0) return 0.0;
        double a = std::pow((double)c, -s);
        double b = std::pow((double)(c+1), -s);
        double v = a - b;
        return (v > 0.0 ? v : 1e-300);
    }

    // dzeta over 1..max_copy, normalized (Zipf with exponent p)
    static std::vector<double> zetaWeights_(double zp, int max_copy) {
        if (zp <= 1.0) throw std::invalid_argument("zp.copy must be > 1 (to be normalizable)");
        std::vector<double> w(max_copy+1, 0.0);
        double S=0.0;
        for (int k=1;k<=max_copy;++k) { w[k] = 1.0 / std::pow((double)k, zp); S += w[k]; }
        for (int k=1;k<=max_copy;++k) w[k] /= S;
        return w; // w[0] unused; sum_{k=1..max_copy} w[k] = 1
    }

    // pgood(x; zp, u_base, sd_base, shape=0) with mixture over copy=1..max_copy
    std::vector<double> pgood_(int max_x, double u_base, double sd_base, const std::vector<double>& zeta_w) const {
        int N = std::max(0, max_x);
        std::vector<double> out(N+1, 0.0);
        for (int x=1; x<=N; ++x) {
            double sum = 0.0;
            for (int copy=1; copy<=P_.max_copy; ++copy) {
                double mu = copy * u_base;
                double sd = std::sqrt((double)copy) * sd_base;
                sum += zeta_w[copy] * normal_pdf_(x, mu, sd);
            }
            out[x] = sum;
        }
        return out;
    }

    // Single x value, hom vs het
    double valHom_(int x) const {
        double sum=0.0;
        for (int copy=1; copy<=P_.max_copy; ++copy) {
            double mu = copy * P_.u_v;
            double sd = std::sqrt((double)copy) * P_.sd_v;
            sum += zeta_hom_[copy] * normal_pdf_(x, mu, sd);
        }
        return std::max(sum, 1e-300);
    }
    double valHet_(int x) const {
        double u_base = 0.5 * P_.u_v;
        double sd_base= sd_het_base_;
        double sum=0.0;
        for (int copy=1; copy<=P_.max_copy; ++copy) {
            double mu = copy * u_base;
            double sd = std::sqrt((double)copy) * sd_base;
            sum += zeta_het_[copy] * normal_pdf_(x, mu, sd);
        }
        return std::max(sum, 1e-300);
    }
};
