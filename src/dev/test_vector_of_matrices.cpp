#include <iostream>
#include <numeric>
#include <vector>
// #include "combinatorics.hpp"
// #include "data.hpp"
// #include "models.hpp"
// #include "parameters.hpp"
// #include "particles.hpp"
// #include "sampling.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXd;
using Eigen::MatrixXi;

int main()
{
    int n_sites = 10;
    int n_allele_configs = 8;
    int n_ibd_states = 5;
    std::vector<MatrixXd> sampling_probs(
        n_sites, 
        MatrixXd::Constant(n_allele_configs, n_ibd_states, -1.0)
    );
    std::cout << sampling_probs.size() << std::endl;
    for (int i=0;i<n_sites;++i)
    {   
        std::cout << "Array: " << i << std::endl;
        std::cout << sampling_probs[i] << std::endl;
        //std:cout << "" << std::endl;
    }
}