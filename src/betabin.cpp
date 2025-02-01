#include <cassert>
#include <cmath>
#include <utility>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"
#include "betabin.hpp"
#include "constants.hpp"
using namespace std;


namespace Betabinomial {

    const std::pair<int, int> missing_pair{MISSING_AD_VALUE, MISSING_AD_VALUE};

    struct pair_hash
    {
        size_t operator()(const pair<int, int> &p) const
        {
            size_t h = (size_t(p.first) << 32) ^ size_t(p.second);
            return h;
        }
    };


    MatrixXd calc_prob_matrix(
        const ArrayXi& refs,
        const ArrayXi& alts,
        const ArrayXd& pis,
        const double v,
        bool as_loglikelihood)
    {   
        // Create the matrix
        int n_sites = refs.size();
        int n_pis = pis.size();
        MatrixXd prob_matrix(MatrixXd::Constant(n_sites, n_pis, DEFAULT_FLOAT));

        // Iterate over all (REF, ALT) pairs
        unordered_map<pair<int, int>, int, pair_hash> ra_pair_map;
        for (int i = 0; i < n_sites; ++i) {

            // Check if you have already computed for this (REF, ALT) pair
            pair<int, int> ra_pair{refs(i), alts(i)};
            
            // Handle case of missing data
            if (ra_pair == missing_pair) {
                if (as_loglikelihood) {
                    prob_matrix.row(i).setZero();
                } else {
                    prob_matrix.row(i).setOnes();
                }
                continue;
            }

            // Check if already computed
            auto found = ra_pair_map.find(ra_pair);
            if (found != ra_pair_map.end()) {
                prob_matrix.row(i) = prob_matrix.row(found->second); // TODO: is this best?
                continue;
            };

            // If not, compute
            double tmp0 = lgamma(alts(i) + refs(i) + 1) - lgamma(alts(i) + 1) - lgamma(refs(i) + 1);
            double tmp1 = lgamma(v) - lgamma(alts(i) + refs(i) + v);

            for (int j = 0; j < n_pis; ++j) {
                double pi = pis(j);
                double tmp2 = lgamma(alts(i) + pi * v) - lgamma(pi * v);
                double tmp3 = lgamma(refs(i) + (1.0 - pi) * v) - lgamma((1.0 - pi) * v);

                if (as_loglikelihood) {
                    prob_matrix(i, j) = tmp0 + tmp1 + tmp2 + tmp3;
                } else {
                    prob_matrix(i, j) = exp(tmp0 + tmp1 + tmp2 + tmp3);
                }
                
                ra_pair_map[ra_pair] = i;
            }
        }
        return prob_matrix; // TODO: would be much better to not copy.
    }


    LookupMatrix::LookupMatrix(
        const VCFData& data, 
        const int n_pi_bins, 
        const double e_0, 
        const double e_1, 
        const double v,
        bool as_loglikelihood)
    : n_pi_bins(n_pi_bins),
    e_0(e_0),
    e_1(e_1),
    v(v),
    pi_bin_midpoints(ArrayXd::LinSpaced(n_pi_bins, e_0, 1 - e_1)), // [e_0, 1 - e_1] evenly spaced.
    prob_matrix(calc_prob_matrix(data.refs, data.alts, pi_bin_midpoints, v, as_loglikelihood))
    {
        check_valid();
    };

    void LookupMatrix::check_valid() const 
    {
        if ((prob_matrix.array() == DEFAULT_FLOAT).any()) {
            throw std::invalid_argument("Failed to completely initialise Betabinomial lookup matrix.");
        }
    }

    MatrixXd LookupMatrix::subset(ArrayXd pi_vals) const
    {
        ArrayXi indices = (
           (n_pi_bins - 1) / (1 - e_1 - e_0) * (pi_vals - e_0)
        ).round().cast<int>();
        return prob_matrix(Eigen::all, indices);
    }

}

