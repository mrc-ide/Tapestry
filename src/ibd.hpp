#pragma once
#include <utility>
#include <vector>
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXi;
using namespace std;


class IBDContainer
{
private:
    vector<int> static create_strains(int K);
    MatrixXi static create_pair_to_column_index(int K);
    vector<pair<int, int>> create_column_index_to_pair(int K);
    vector<vector<int>> static create_state_index_to_column_indices(
        const vector<vector<vector<int>>>& states,
        const MatrixXi& pair_to_column_index
    );

public:
    // Constructor
    IBDContainer(int K);

    // Members

    // Integer vector of strains, e.g. if K=3, {0, 1, 2}
    const vector<int> strains;

    // All possible IBD states
    // size = Bell(k)
    // correspond to hidden states in HMM
    // first index is *state*
    // second inndex is *strain group*; 
    // e.g. if K=3 ...
    // {{{0},{1},{2}}, {{0,1}, {2}}, {{0,2},{1}}, {{0},{1,2}}, {{0,1,2}}
    const vector<vector<vector<int>>> states;

    // We also have a notion of (K choose 2) IBD pairs
    // these pairs have a (semi) arbitrary order (the 'column index')
    // so here... pair_to_column_index(0, 1) = 2
    // ... would mean {0, 1} pair has column index 2
    const MatrixXi pair_to_column_index;

    // Return the pair given a column index
    const vector<pair<int, int>> column_index_to_pair;

    // Go from IBD state index to all implied pairwise column indices
    // Note this will be ragged
    const vector<vector<int>> state_index_to_column_indices;
};

