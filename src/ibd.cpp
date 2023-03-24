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
    column_index_to_pair(get_unique_pairs(strains)),
    pair_to_column_index(create_pair_to_column_index(strains)),
    state_index_to_column_indices(create_state_index_to_column_indices(states, pair_to_column_index))
{}


vector<int> IBDContainer::create_strains(int K)
{
    vector<int> strains(K);
    iota(strains.begin(), strains.end(), 0);
    return strains;
}


MatrixXi IBDContainer::create_pair_to_column_index(const vector<int>& strains)
{
    int K = strains.size();
    vector<pair<int, int>> pairs = get_unique_pairs(strains);
    MatrixXi m = MatrixXi::Constant(K, K, -1);
    
    for (int c = 0; c < pairs.size(); ++c) {
        m(pairs[c].first, pairs[c].second) = c;
        m(pairs[c].second, pairs[c].first) = c;
    }
    
    return m;
}


vector<vector<int>> IBDContainer::create_state_index_to_column_indices(
    const vector<vector<vector<int>>>& states,
    const MatrixXi& pair_to_column_index
)
{
    // Iterate over all IBD states
    vector<vector<int>> mapping;
    for (int j = 0; j < states.size(); ++j) {
        
        vector<int> state_indices;
        
        for (const vector<int>& ibd_group : states[j]) {
            
            // No pairs if one strain
            if (ibd_group.size() == 1) {
                continue;
            }
            
            // Otherwise add pairs
            for (const pair<int, int>& pair : get_unique_pairs(ibd_group)) {
                state_indices.push_back(
                    pair_to_column_index(pair.first, pair.second)
                );
            }
        }
        mapping.push_back(state_indices);
    }
    
    return mapping;
}


