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
#include "libs/eigen-3.4.0/Dense"
#include "typedefs.hpp"
using namespace std;




namespace Betabinomial {

    /* Identify missing data */
    extern const std::pair<int,int> missing_pair;

    /* Used to hash the (REF, ALT) counts */
    struct pair_hash;

    MatrixXd calc_prob_matrix(
        const ArrayXi& refs,
        const ArrayXi& alts,
        const ArrayXd& pis,
        const double v,
        bool as_loglikelihood = false
    );

    class LookupMatrix {
    private:
        const int n_pi_bins;
        const double e_0;
        const double e_1;
        const double v;

        ArrayXd pi_bin_midpoints;
        MatrixXd prob_matrix;

        void check_valid() const; // check all values initialised
    public:
        LookupMatrix(
            const VCFData& data,
            const int n_pi_bins,
            const double e_0,
            const double e_1,
            const double v,
            bool as_loglikelihood = true);
        
        // Subset to the columns with midpoints closest to
        // indicated pi_vals
        MatrixXd subset(ArrayXd pi_vals) const;
    };
}

