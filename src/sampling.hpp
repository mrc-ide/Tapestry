#pragma once
#include <vector>
#include "typedefs.hpp"
#include "libs/eigen-3.4.0/Dense"
using namespace std;


/*
* For a biallelic SNP with PLAF `p`, calculate the probability of 
* sampling different IBD states and allele configurations, given an
* assumed COI.
*/
MatrixXd calc_sampling_probs(
    double p,
    const MatrixXi& allele_configs,
    const vector<vector<vector<int>>>& ibd_states
);


/*
*  Same as above, but compute matrices for an array of PLAF values
*  and store output in a vector
*  TODO: could also implement with plafs as a vector
*/
vector<MatrixXd> create_sampling_probs(
    const ArrayXd& plafs,
    const MatrixXi& allele_configs,
    const vector<vector<vector<int>>>& ibd_states
);
