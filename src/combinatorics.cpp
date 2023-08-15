#include <algorithm>
#include <vector>
#include <utility>
#include <cmath>
#include <stdexcept>
#include "combinatorics.hpp"
using namespace std;




MatrixXi create_powerset(int K)
{
    if (K < 1) {
        throw std::invalid_argument("Number of strains, `K`, must be one or greater.");
    }

    // Initialise
    int N = std::pow(2, K);  // Number of sets in powerset
    MatrixXi powerset = MatrixXi::Constant(N, K, 0);

    // Populate
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < K; ++j) {
            int val = i / pow(2, j);
            if (val % 2) {
                powerset(i, j) = 1;
            }
        }
    }

    return powerset;
}


vector<vector<vector<int>>> create_all_partitions(vector<int> collection)
{
    // A set of one element has only one partition
    if (collection.size() == 1) {
        vector<vector<vector<int>>> n1_partitions(
            1, vector<vector<int>>(
            1, vector<int>(
            1)
            )
        );
        n1_partitions[0][0][0] = collection[0];
        return n1_partitions;
    }

    // A set of n + 1 elements is created recurisvely
    //int elem = collection.front();
    //collection.erase(collection.begin());
    int elem = collection.back();
    collection.pop_back();

    vector<vector<vector<int>>> n_partitions = create_all_partitions(collection);
    vector<vector<vector<int>>> n1_partitions;
    
    for (vector<vector<int>> partition : n_partitions) {
        
        for (int i = 0; i < partition.size(); ++i) {
            vector<vector<int>> n1_partition = partition;
            n1_partition[i].push_back(elem);
            n1_partitions.push_back(n1_partition);
        }

        vector<vector<int>> n1_partition = partition;
        vector<int> only_elem_subset;
        only_elem_subset.push_back(elem);
        n1_partition.push_back(only_elem_subset);
        n1_partitions.push_back(n1_partition);
    }

    return n1_partitions;
}


vector<pair<int, int>> get_unique_pairs(vector<int> elements)
{
    // First, reduce to unique elements
    sort(elements.begin(), elements.end());
    auto it = unique(elements.begin(), elements.end());
    elements.resize(distance(elements.begin(), it));

    // Then get pairs
    int n_pairs = elements.size() * (elements.size() - 1) / 2;
    vector<pair<int,int>> pairs;
    pairs.reserve(n_pairs);
    for (int i = 0; i < elements.size(); ++i) {
        for (int j = i + 1; j < elements.size(); ++j) {
            pairs.emplace_back(pair<int,int>(elements[i],elements[j]));
        }
    }

    return pairs;
}

