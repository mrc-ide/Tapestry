#pragma once
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "typedefs.hpp"
using namespace std;


// --------------------------------------------------------------------------------
// Abstract base class
// We are allowinng for multiple implementations of the betabinomial look-up
// array.
// --------------------------------------------------------------------------------


class BetabinomialArray
{
public:
    // DATA
    int n_loci;
    vector<int> refs;
    vector<int> alts;

    // PARAMETERS
    int n_pi_bins;
    double e_0;
    double e_1;
    double v;

    // CONSTRUCTORS
    BetabinomialArray();
    BetabinomialArray(
        int n_loci, 
        vector<int> refs, 
        vector<int> alts);
    BetabinomialArray(
        int n_loci,
        vector<int> refs, 
        vector<int> alts, 
        int n_pi_bins, 
        double e_0, 
        double e_1, 
        double v);

    // CONCRETE METHODS
    double operator()(int locus, double pi_val) const;

    // ABSTRACT METHODS
    virtual void compute_array(bool as_loglikelihood = true) = 0;

protected:
    // LOOKUP ARRAY
    matrix_2d_double betabin_array;
};



// --------------------------------------------------------------------------------
// Concrete implementations
//
//
// --------------------------------------------------------------------------------


class BetabinomialArrayByHash : public BetabinomialArray
{
private:
    // TODO: 
    // - Ensure no possible collisions here
    // - Is it okay to define this here?
    struct pair_hash
    {
        size_t operator()(const pair<int, int> &p) const
        {
            size_t h = (size_t(p.first) << 32) ^ size_t(p.second);
            return h;
        }
    };

public:
    // CONSTRUCTORS
    BetabinomialArrayByHash();
    BetabinomialArrayByHash(
        int n_loci, 
        vector<int> refs, 
        vector<int> alts);
    BetabinomialArrayByHash(
        int n_loci,
        vector<int> refs, 
        vector<int> alts, 
        int n_pi_bins, 
        double e_0, 
        double e_1, 
        double v);

    // METHODS
    void compute_array(bool as_loglikelihood = true) override;
};

