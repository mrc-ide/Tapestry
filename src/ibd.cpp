#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "bed.hpp"
#include "combinatorics.hpp"
#include "ibd.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;


// ================================================================================
// Container for IBD information
//
// ================================================================================


IBDContainer::IBDContainer()
{}

IBDContainer::IBDContainer(int K)
    : strains(create_strains(K)),
    states(create_all_partitions(strains)),
    column_index_to_pair(get_unique_pairs(strains)),
    pair_names(create_pair_names(column_index_to_pair)),
    pair_to_column_index(create_pair_to_column_index(strains)),
    state_index_to_column_indices(create_state_index_to_column_indices(states, pair_to_column_index))
{}


vector<int> IBDContainer::create_strains(int K)
{
    vector<int> strains(K);
    iota(strains.begin(), strains.end(), 0);
    return strains;
}

vector<string>  IBDContainer::create_pair_names(const vector<pair<int, int>>& pairs)
{
    //string pair_name;
    vector<string> pair_names;
    for (const pair<int, int>& pair : pairs) {
        ostringstream pair_name;
        pair_name << "strains" << pair.first << "-" << pair.second;
        pair_names.push_back(pair_name.str());
    }
    return pair_names;
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


// ================================================================================
// Conversion from HMM IBD states to pairwise IBD matrix
//
// ================================================================================



MatrixXi convert_ibd_state_path_to_pairwise(
    const VectorXi& ibd_path, 
    const IBDContainer& ibd
)
{
    // Define size
    int n_sites = ibd_path.rows();
    int n_pairs = ibd.column_index_to_pair.size();

    // Initialise with zeros
    MatrixXi m = MatrixXi::Constant(n_sites, n_pairs, 0);

    // Iterate along vector and populate
    vector<int> indices;
    for (int i = 0; i < n_sites; ++i) {
        indices = ibd.state_index_to_column_indices[ibd_path[i]];
        m(i, indices) = VectorXi::Constant(indices.size(), 1);
    }

    return m;
}


// ================================================================================
// IBD Segment Detection
//
// ================================================================================



// TODO:
// - What is the most efficient way to populate the vector of records?
// - What about copying it out of the scope?
std::vector<BEDRecord> get_ibd_segments(
    std::vector<std::string> chroms,
    VectorXi pos,
    VectorXi ibd_states, // TODO: do I need this to be explicitly bool?
    std::string name
)
{
    // Init
    int t = 0;
    int start, end;
    std::string chrom = chroms[t];
    int ibd_state = ibd_states(t);
    if (ibd_state) {
        start = pos(t);
    }

    // Iterate
    //BEDRecord record;
    vector<BEDRecord> bed_records;  // we don't know the size a priori
    for (; t < pos.size(); ++t) {
        
        // Handle end of chromosome
        if (chroms[t] != chrom) {
            if (ibd_state) {
                end = pos(t-1);
                bed_records.push_back(BEDRecord(chrom, start, end, name));
            }
            if (ibd_states(t)) {
                start = pos(t);
            }
            chrom = chroms[t];
            ibd_state=ibd_states(t);
            continue;
        }

        // Handle chromosome internal
        if (ibd_state) {
            if (!ibd_states(t)) {  // segment ending
                end = (pos(t-1) + pos(t)) / 2;
                bed_records.push_back(BEDRecord(chrom, start, end, name));
            } 
        } else if (ibd_states(t)) { // segment starting
            start = (pos(t-1) + pos(t)) / 2;
        }

        // Update memory
        chrom = chroms[t];
        ibd_state = ibd_states(t);
    }

    // Terminate
    if (ibd_state) {
        end = pos(t-1); // As with end of chromosome.
        bed_records.push_back(BEDRecord(chrom, start, end, name));
    }

    return bed_records;
}

