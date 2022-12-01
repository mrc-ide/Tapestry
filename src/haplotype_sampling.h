#pragma once
#include <vector>
using namespace std;

matrix_2d_double calc_sampling_probs(
    double p,
    const vector<vector<vector<int>>>& ibd_configs,
    const matrix_2d_int& hap_configs);