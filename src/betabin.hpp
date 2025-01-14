#pragma once
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "constants.hpp"
#include "data.hpp"
#include "parameters.hpp"
#include "typedefs.hpp"
using namespace std;


class BetabinomialArray
{
private:

    /*
    * Hash a pair of integers;
    * used to hash the (REF, ALT) read counts
    */
    struct pair_hash
    {
        size_t operator()(const pair<int, int> &p) const
        {
            size_t h = (size_t(p.first) << 32) ^ size_t(p.second);
            return h;
        }
    };

    //const Parameters& params;
    const VCFData& data;
    const int n_pi_bins;    // Number of bins to approximate
    const double e_0;       // REF->ALT error rate
    const double e_1;       // ALT->REF error rate
    const double v;         // Overdispersion parameter
    

    MatrixXd lookup_matrix;
    const pair<int, int> missing_pair{MISSING_AD_VALUE, MISSING_AD_VALUE};
    void calc_lookup_matrix(bool as_loglikelihood);
    void check_valid() const; // Check that all values have been initialised.

public:

    const bool as_loglikelihood; // defaults to true, first constructor

    // BetabinomialArray(const Parameters& params, const VCFData& data);
    // BetabinomialArray(const Parameters& params, const VCFData& data, bool as_loglikelihood);

    BetabinomialArray(
        const VCFData& data, 
        const int n_pi_bins, 
        const double e_0, 
        const double e_1, 
        const double v);
    BetabinomialArray(
        const VCFData& data, 
        const int n_pi_bins, 
        const double e_0, 
        const double e_1, 
        const double v,
        bool as_loglikelihood
        );

    /*
    * The () operator is overloaded to act like Eigen array indices, e.g.
    * array(row_index, col_index); with an important exception:
    * Instead of taking two integers, the column index is a double that
    * gets rounded to the approproriate array index
    * 
    */
    double operator()(int locus, double pi_val) const;

    /*
    * Given a vector of WSAF values (pi_vals), return a subsetted array
    * TODO: what to generalise to ArrayBase
    */
    MatrixXd subset(ArrayXd pi_vals) const;
};

