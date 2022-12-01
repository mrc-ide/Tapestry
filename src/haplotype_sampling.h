#pragma once
#include <vector>
using namespace std;

matrix_2d_double calc_sampling_probs(
    double p,
    vector<vector<vector<int>>> ibd_configs,
    matrix_2d_int hap_configs);