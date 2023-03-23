#include <utility>
#include <numeric>
#include <vector>
#include "combinatorics.hpp"
#include "ibd.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXi;
using namespace std;


IBDContainer::IBDContainer(int K)
    : strains(create_strains(K)),
    states(create_all_partitions(strains)),
    pair_to_column_index(create_pair_to_column_index(K)),
    column_index_to_pair(create_column_index_to_pair(K)),
    state_index_to_column_indices(create_state_index_to_column_indices(states, pair_to_column_index))
{};


vector<int> IBDContainer::create_strains(int K)
{
    vector<int> strains(K);
    iota(strains.begin(), strains.end(), 0);
    return strains;
}


MatrixXi IBDContainer::create_pair_to_column_index(int K)
{
    MatrixXi m = MatrixXi::Constant(K, K, -1);

    int z = 0;
    for (int i = 0; i < K; ++i) {
        for (int j = i + 1; j < K; ++j) {
            m(i, j) = z;
            m(j, i) = z;
            ++z;
        }
    }

    return m;
}


// TODO: code duplication with above
vector<pair<int, int>> IBDContainer::create_column_index_to_pair(int K)
{
    vector<pair<int,int>> mapping;

    int z = 0;
    for (int i = 0; i < K; ++i) {
        for (int j = i + 1; j < K; ++j) {
            mapping.push_back(pair<int,int>(i, j));
            ++z;
        }
    }

    return mapping;
}


vector<vector<int>> IBDContainer::create_state_index_to_column_indices(
    const vector<vector<vector<int>>>& states,
    const MatrixXi& pair_to_column_index
)
{
    // Prepare the empty mapping
    vector<vector<int>> mapping;

    // Iterate over IBD stats
    // TODO: this works but horrifying quadruple loop
    for (int j = 0; j < states.size(); ++j) {
        
        vector<int> state_indices;
        
        for (const vector<int>& ibd_group : states[j]) {
            
            // Check if this IBD group has  > 1 strain
            int n_strains = ibd_group.size();
            if (n_strains == 1) {
                continue;
            }
            
            // If so, iterate over all pairs, get index, and add to vector
            for (int a = 0; a < n_strains; ++a) {
                for (int b = a + 1; b < n_strains; ++b) {
                    state_indices.push_back(pair_to_column_index(ibd_group[a], ibd_group[b]));
                }
            }
        }

        mapping.push_back(state_indices);
    }

    return mapping;
}


