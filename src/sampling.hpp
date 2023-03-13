#pragma once
#include <vector>
#include "typedefs.hpp"
using namespace std;


/*
* For a biallelic SNP with PLAF `p`, calculate the probability of 
* sampling different IBD states and allele configurations, given an
* assumed COI.
*/
MatrixXd calc_sampling_probs(
    double p,
    const vector<vector<vector<int>>>& ibd_states,
    const MatrixXi& allele_configs);

