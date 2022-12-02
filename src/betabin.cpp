#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"
#include "betabin.hpp"
using namespace std;


// --------------------------------------------------------------------------------
// Abstract base class
// We are allowinng for multiple implementations of the betabinomial look-up
// array.
// --------------------------------------------------------------------------------


BetabinomialArray::BetabinomialArray() {}


BetabinomialArray::BetabinomialArray(
        int n_loci,
        vector<int> refs,
        vector<int> alts
)
    : n_loci(n_loci), 
    refs(refs),
    alts(alts),
    n_pi_bins(100),
    e_0(0.01),
    e_1(0.05),
    v(500),
    betabin_array(boost::extents[n_loci][n_pi_bins])
{
}


BetabinomialArray::BetabinomialArray(
        int n_loci,
        vector<int> refs,
        vector<int> alts,
        int n_pi_bins,
        double e_0,
        double e_1,
        double v
)
    : n_loci(n_loci), 
    refs(refs),
    alts(alts),
    n_pi_bins(n_pi_bins),
    e_0(e_0),
    e_1(e_1),
    v(v),
    betabin_array(boost::extents[n_loci][n_pi_bins])
{
}


// TODO: 
// - Probably should also be pure virtual
// - Bounds checks needed, e.g.
//  - If pi < e_0 gives a negative index
//  - If pi = 1 - e_1 gives index = size, instead of size -1
double BetabinomialArray::operator()(int locus, double pi_val) const
{
    int ix = round(
        n_pi_bins * (pi_val - e_0)/(1 - e_1 - e_0)
    );
    cout << "Index: " << ix << endl;
    return betabin_array[locus][ix];
}


// --------------------------------------------------------------------------------
// Concrete implementations
//
// 
// --------------------------------------------------------------------------------


BetabinomialArrayByHash::BetabinomialArrayByHash() {}


BetabinomialArrayByHash::BetabinomialArrayByHash(
        int n_loci,
        vector<int> refs,
        vector<int> alts
)
    : BetabinomialArray(n_loci, refs, alts)
{
}


BetabinomialArrayByHash::BetabinomialArrayByHash(
        int n_loci,
        vector<int> refs,
        vector<int> alts,
        int n_pi_bins,
        double e_0,
        double e_1,
        double v
)
    : BetabinomialArray(n_loci, refs, alts, n_pi_bins, e_0, e_1, v)
{
}


void BetabinomialArrayByHash::compute_array(bool as_loglikelihood)
{
    // Compute pi values for each bin
    // TODO: Don't quite reach top bin
    // Might just want to let this be [0, 1] -- wasting some compute though
    vector<double> pi_vals(n_pi_bins);
    for (int i = 0; i < n_pi_bins; ++i) {
        pi_vals[i] = i * (1 - e_1 - e_0) / n_pi_bins + e_0;
    }

    // Avoid recalculation for identical (REF, ALT) pairs
    // Key: (REF, ALT) pair
    // Value: Row index in this->betabin_array
    unordered_map<pair<int, int>, int, pair_hash> ra_pair_map;

    for (int i = 0; i < n_loci; ++i) {        
        // If (REF, ALT) pair already observed, use those values
        pair<int, int> ra_pair(refs[i], alts[i]);
        auto found = ra_pair_map.find(ra_pair);
        if (found != ra_pair_map.end()) {
            betabin_array[i] = betabin_array[found->second]; // TODO: double check this
            continue;
        }

        // Otherwise, compute
        double tmp0 = lgamma(alts[i] + refs[i] + 1) - lgamma(alts[i] + 1) - lgamma(refs[i] + 1);
        double tmp1 = lgamma(v) - lgamma(alts[i] + refs[i] + v);
        for (int j = 0; j < n_pi_bins; ++j) {
            
            const double& pi = pi_vals[j];
            double tmp2 = lgamma(alts[i] + pi*v) - lgamma(pi*v);
            double tmp3 = lgamma(refs[i] + (1.0 - pi)*v) - lgamma((1.0 - pi)*v);

            if (as_loglikelihood) {
                betabin_array[i][j] = tmp0 + tmp1 + tmp2 + tmp3;
            } else {
                betabin_array[i][j] = exp(tmp0 + tmp1 + tmp2 + tmp3);
            }

            // Store this  (REF, ALT) pair
            ra_pair_map[ra_pair] = i;
        }
    }
}

