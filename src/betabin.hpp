#pragma once
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "data.hpp"
#include "parameters.hpp"
#include "typedefs.hpp"
using namespace std;


// Create a class that encapsulates a betabinomial array
// - during development, we may produce different types of arrays
// - we want them all to behave similiarly
// - could even profile lookup vs. no lookup
// - what would API be?

// betabin.compute_emission(pi_val);
// - But that is different from overloading the () operator

// - We have Factories that produce these different arrays



// Create a class that is a betabinomial array
// - I want it to behave close to a normal Eigen array
// - But it takes float values


/*
* Create a `Betabinomial` array
*
*
*/
// class BetabinomialArrayFactory
// {};



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

    const Parameters& params;
    const VCFData& data;

    MatrixXd lookup_matrix;


    void calc_lookup_matrix(bool as_loglikelihood);

public:

    const bool as_loglikelihood; // defaults to true, first constructor

    BetabinomialArray(const Parameters& params, const VCFData& data);
    BetabinomialArray(const Parameters& params, const VCFData& data, bool as_loglikelihood);

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

