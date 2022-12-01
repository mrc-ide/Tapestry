#pragma once
#include <vector>
#include "typedefs.hpp"
using namespace std;

matrix_2d_int create_powerset(int K);
vector<vector<vector<int>>> create_all_partitions(vector<int> collection);