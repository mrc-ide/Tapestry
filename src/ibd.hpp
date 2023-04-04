#pragma once
#include <string>
#include <utility>
#include <vector>
#include "bed.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXi;
using Eigen::VectorXi;
using namespace std;


// ================================================================================
// Container for IBD information
//
// ================================================================================


class IBDContainer
{
private:
    vector<int> static create_strains(int K);
    vector<string> static create_pair_names(const vector<pair<int, int>>& pairs);
    MatrixXi static create_pair_to_column_index(const vector<int>& strains);
    vector<vector<int>> static create_state_index_to_column_indices(
        const vector<vector<vector<int>>>& states,
        const MatrixXi& pair_to_column_index
    );

public:
    IBDContainer(int K);

    /* Strain indices
    *  e.g., if K=3, {0, 1, 2}
    */
    const vector<int> strains;

    /* All possible IBD states
    *  Here, by IBD state, we mean a particular *hidden state* in the
    *  HMM that defines our model. For higher K, these IBD states 
    *  define *pairs* or *groups* of strains that are in IBD.
    *  There are Bell(K) number of IBD states.
    *  e.g. if K=3...
    *
    *  {{{0},{1},{2}}, {{0,1}, {2}}, {{0,2},{1}}, {{0},{1,2}}, {{0,1,2}}
    *
    */
    const vector<vector<vector<int>>> states;

    /* For each column in a pairwise IBD matrix, which *pair* of strains
    *  is in IBD?
    */
    const vector<pair<int, int>> column_index_to_pair;

    /* Names for pairs of strains in IBD
    *  e.g., strains0-1
    */
    const vector<string> pair_names;


    /* Inverse of `column_index_to_pair`
    *  We have  a notion  of (K choose 2) IBD pairs. We give these pairs
    *  a (semi-arbitary) order, defined in `column_index_to_pair`.
    *  Here, pair_to_column_index(0, 1) = 2 would mean IBD state information
    *  about strains {0, 1} is stored in column 2 of the pairwise matrix.
    */
    const MatrixXi pair_to_column_index;

    /* Map from IBD state (i.e. from HMM) to column indices of
    *  all pairs of strains that are in IBD.
    *
    */
    const vector<vector<int>> state_index_to_column_indices;
};


// ================================================================================
// IBD Segment Detection
//
// ================================================================================


std::vector<BEDRecord> get_ibd_segments(
    std::vector<std::string> chroms,
    VectorXi pos,
    VectorXi ibd_states, // TODO: do I need this to be explicitly bool?
    std::string name
);