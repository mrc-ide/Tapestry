#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"
#include "betabin.hpp"
using namespace std;


BetabinomialArray::BetabinomialArray(const Data& data, const Parameters& params)
    : params(params),
    loglookup(calc_loglookup_matrix(data, params))
{};


MatrixXd BetabinomialArray::calc_loglookup_matrix(const Data& data, const Parameters& params)
{
    
    // Create storage
    MatrixXd loglookup = MatrixXd::Constant(data.n_loci, params.n_pi_bins, 9999.0);
    
    // Create array of error adjusted WSAF, at which we will precompute probs from Betabin
    ArrayXd pi_bin_midpoints = ArrayXd::LinSpaced(params.n_pi_bins, params.e_0, 1 - params.e_1); // TODO: double check length

    // Iterate over all (REF, ALT) pairs
    unordered_map<pair<int, int>, int, pair_hash> ra_pair_map;
    for (int i = 0; i < data.n_loci; ++i) {

        // Check if you have already computed for this (REF, ALT) pair
        pair<int, int> ra_pair(data.refs[i], data.alts[i]);
        auto found = ra_pair_map.find(ra_pair);
        if (found != ra_pair_map.end()) {
            loglookup.row(i) = loglookup.row(found->second); // TODO: is this best?
            continue;
        };

        // If not, compute probabilities
        // TODO: there could exist faster ways to compute
        double tmp0 = lgamma(data.alts[i] + data.refs[i] + 1) - lgamma(data.alts[i] + 1) - lgamma(data.refs[i] + 1);
        double tmp1 = lgamma(params.v) - lgamma(data.alts[i] + data.refs[i] + params.v);

        for (int j = 0; j < params.n_pi_bins; ++j) {
            double pi = pi_bin_midpoints(j);
            double tmp2 = lgamma(data.alts[i] + pi * params.v) - lgamma(pi * params.v);
            double tmp3 = lgamma(data.refs[i] + (1.0 - pi) * params.v) - lgamma((1.0 - pi) * params.v);

            loglookup(i, j) = tmp0 + tmp1 + tmp2 + tmp3;
            ra_pair_map[ra_pair] = i;
        }
    }

    // TODO: double check that all the values have been chaged for 9999.0
    // I think two step bool is_zero = (loglookup.array() > 0).any(); something like this exists

    return loglookup;
}


double BetabinomialArray::operator()(int locus, double pi_val) const
{
    // TODO: there are other ways to round in C++
    // this function called *millions of times*; optimise is good
    // Using the params reference might also slow down

    int ix = round(
        (params.n_pi_bins - 1) / (1 - params.e_1 - params.e_0) * (pi_val - params.e_0)
    );
    
    return loglookup(locus, ix);
}

