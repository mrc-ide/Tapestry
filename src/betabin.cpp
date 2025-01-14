#include <cassert>
#include <cmath>
#include <utility>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"
#include "betabin.hpp"
#include "constants.hpp"
using namespace std;


// BetabinomialArray::BetabinomialArray(const Parameters& params, const VCFData& data)
//     : params(params),
//     data(data),
//     as_loglikelihood(true),
//     lookup_matrix(MatrixXd::Constant(data.n_sites, params.n_pi_bins, DEFAULT_FLOAT))
// {
//     calc_lookup_matrix(as_loglikelihood);
//     check_valid();
// };


// BetabinomialArray::BetabinomialArray(const Parameters& params, const VCFData& data, bool as_loglikelihood)
//     : params(params),
//     data(data),
//     as_loglikelihood(as_loglikelihood),
//     lookup_matrix(MatrixXd::Constant(data.n_sites, params.n_pi_bins, DEFAULT_FLOAT))
// {
//     calc_lookup_matrix(as_loglikelihood);
//     check_valid();
// };


BetabinomialArray::BetabinomialArray(
        const VCFData& data, 
        const int n_pi_bins, 
        const double e_0, 
        const double e_1, 
        const double v)
    : data(data),
    n_pi_bins(n_pi_bins),
    e_0(e_0),
    e_1(e_1),
    v(v),
    as_loglikelihood(true),
    lookup_matrix(MatrixXd::Constant(data.n_sites, n_pi_bins, DEFAULT_FLOAT))
{
    calc_lookup_matrix(as_loglikelihood);
    check_valid();
};

BetabinomialArray::BetabinomialArray(
        const VCFData& data, 
        const int n_pi_bins, 
        const double e_0, 
        const double e_1, 
        const double v,
        bool as_loglikelihood)
    : data(data),
    n_pi_bins(n_pi_bins),
    e_0(e_0),
    e_1(e_1),
    v(v),
    as_loglikelihood(as_loglikelihood),
    lookup_matrix(MatrixXd::Constant(data.n_sites, n_pi_bins, DEFAULT_FLOAT))
{
    calc_lookup_matrix(as_loglikelihood);
    check_valid();
};

void BetabinomialArray::calc_lookup_matrix(bool as_loglikelihood)
{

    // Create array of error adjusted WSAF, at which we will precompute probs from Betabin
    ArrayXd pi_bin_midpoints = ArrayXd::LinSpaced(
        n_pi_bins, 
        e_0, 1 - e_1
        ); // TODO: double check length

    // Iterate over all (REF, ALT) pairs
    unordered_map<pair<int, int>, int, pair_hash> ra_pair_map;
    for (int i = 0; i < data.n_sites; ++i) {

        // Check if you have already computed for this (REF, ALT) pair
        pair<int, int> ra_pair{data.refs(i), data.alts(i)};
        
        // Handle case of missing data
        if (ra_pair == missing_pair) {
            if (as_loglikelihood) {
                lookup_matrix.row(i).setZero();
            } else {
                lookup_matrix.row(i).setOnes();
            }
            continue;
        }

        // Check if already computed
        auto found = ra_pair_map.find(ra_pair);
        if (found != ra_pair_map.end()) {
            lookup_matrix.row(i) = lookup_matrix.row(found->second); // TODO: is this best?
            continue;
        };

        // If not, compute probabilities
        // TODO: is this math construction best?
        // TODO: Here is where we are computing the likelihood;
        double tmp0 = lgamma(data.alts(i) + data.refs(i) + 1) - lgamma(data.alts(i) + 1) - lgamma(data.refs(i) + 1);
        double tmp1 = lgamma(v) - lgamma(data.alts(i) + data.refs(i) + v);

        for (int j = 0; j < n_pi_bins; ++j) {
            double pi = pi_bin_midpoints[j];
            double tmp2 = lgamma(data.alts(i) + pi * v) - lgamma(pi * v);
            double tmp3 = lgamma(data.refs(i) + (1.0 - pi) * v) - lgamma((1.0 - pi) * v);

            if (as_loglikelihood) {
                lookup_matrix(i, j) = tmp0 + tmp1 + tmp2 + tmp3;
            } else {
                lookup_matrix(i, j) = exp(tmp0 + tmp1 + tmp2 + tmp3);
            }
            
            ra_pair_map[ra_pair] = i;
        }
    }
}

void BetabinomialArray::check_valid() const
{
    if ((lookup_matrix.array() == DEFAULT_FLOAT).any()) {
        throw std::invalid_argument("Failed to completely initialise Beta-binomial lookup.");
    }
}


// TODO: this is *not* the fastest way to get what I need out of this array
// - Much better would be to have a function:
//  VectorXd pi_val_ixs = convert_pi_vals_to_ixs(pi_vals);
// - Then use the indices to slice the array
// - Or just directly returned a slice version of the array given some pi_vals
// - One I have that array, I *can* martix multiply
double BetabinomialArray::operator()(int locus, double pi_val) const
{
    // TODO: there are other ways to round in C++
    // this function called *millions of times*; optimise is good
    // Using the params reference might also slow down

    int ix = round(
        (n_pi_bins - 1) / (1 - e_1 - e_0) * (pi_val - e_0)
    );
    
    return lookup_matrix(locus, ix);
}

MatrixXd BetabinomialArray::subset(ArrayXd pi_vals) const
{
    ArrayXi indices = (
        (n_pi_bins - 1) / (1 - e_1 - e_0) * (pi_vals - e_0)
    ).round().cast<int>();
    
    return lookup_matrix(Eigen::all, indices);
}

