#include "io.hpp"
#include "ibd.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXi;


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

