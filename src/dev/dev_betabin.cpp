#include "parameters.hpp"
#include "typedefs.hpp"
#include "libs/eigen-3.4.0/Dense"
#include <iostream>
using namespace std;
using Eigen::ArrayXd;


// Goals:
// -------------------------
// - Create this large array with a minimum number of:
// (i) copies
// (ii) computations
// - Use this array with cleanest possible API
// (i) acts like any other normal array (it is an array type)
// (ii) is very fast
// (iii) Minimise approximation error

// In order to minimise approximation error, we probably *do* want to directly include
// e_0 and 1 - e_1; since these betabinomials are present regardless of the proportion updaet



// Two options for this array:
// ---------------------------
// (1)
// - Vector of bin *midpoints*
// - Round to the nearest *midpoint*
// (2)
// - Vector of bin *boundaries*
// - Find which bin value falls within
//
// In both cases, the answers would be equivalent
// Is one much faster than the other


int main(int argc, char* argv[])
{
    // INPUT PARAMETERS
    // array creation
    double e_0 = 0.01;
    double e_1 = 0.05;
    int n_pi_bins = 20;

    // array population
    // ref, alt
    double v = 20;  // helps parameterise betabinomial

    // Compute bin midpoints
    ArrayXd bin_midpoints = ArrayXd::LinSpaced(n_pi_bins, e_0, 1 - e_1);
    cout << bin_midpoints << endl;
    cout << bin_midpoints.size() << endl;

    

}

