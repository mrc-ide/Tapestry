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


class BetabinomialArray
{
private:
    // Store (REF, ALT) pairs as hash
    struct pair_hash
    {
        size_t operator()(const pair<int, int> &p) const
        {
            size_t h = (size_t(p.first) << 32) ^ size_t(p.second);
            return h;
        }
    };

    // Parameters
    const Parameters& params;

    // Betabinomial probability lookup matrix
    const MatrixXd loglookup;

    // Generate the matrix
    MatrixXd static calc_loglookup_matrix(const Data& data, const Parameters& params);

public:

    // Constructor
    BetabinomialArray(const Data& data, const Parameters& params);

    /*
    * The () operator is overloaded to act like Eigen array indices, e.g.
    * array(row_index, col_index); with an important exception:
    * Instead of taking two integers, the column index is a double that
    * gets rounded to the approproriate array index
    * 
    */
    double operator()(int locus, double pi_val) const;
};
