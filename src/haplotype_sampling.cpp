#include <vector>
#include <numeric>
#include "haplotype_sampling.h"
#include <iostream>
using namespace std;


/**
 * For a bi-allelic SNP with a PLAF `p`, calculate the probability
 * of observing a given haplotype configuration, given the IBD 
 * configuration
 * 
 * Corresponds to P(h | p, S) for all h and S.
 * 
 * @params p : the PLAF of the SNP in [0, 1].
 * @params ibd_configs : all possible IBD configurations amongst
 * the strains.
 * @params hap_configs : all possible haplotype configurations
 * amongst the strains. 
*/
vector<vector<double>> calc_sampling_probs(
    double p,
    vector<vector<vector<int>>> ibd_configs,
    vector<vector<int>> hap_configs)
{
    // Count nunmber of IBD and haplotypes configurations
    int n_ibd_configs = ibd_configs.size();
    int n_hap_configs = hap_configs.size();

    // Create sampling probability matrix for each combination
    vector<vector<double>> sampling_probs(
        n_hap_configs, vector<double>(
        n_ibd_configs, 1.0)
    );

    // Populate
    for (int i = 0; i < n_hap_configs; ++i) {
        for (int j = 0; j < n_ibd_configs; ++j) {

            // Collect IBD configuration
            vector<vector<int>>& ibd_config = ibd_configs[j];

            // Iterate over IBD groups within configuration
            for (vector<int>& ibd : ibd_config) {

                // Compute number of strains and ALT alleles in IBD group  
                int n_strains = ibd.size();
                int n_alt = 0;
                for (int& ix : ibd) {
                    n_alt += hap_configs[i][ix];
                }

                // Compute sampling probability of this IBD group
                if (n_strains == n_alt) {
                    sampling_probs[i][j] *= p;
                } else if (n_alt == 0) {
                    sampling_probs[i][j] *= (1 - p);
                } else {
                    sampling_probs[i][j] = 0;
                    break; // IBD and haplotype configs incompatible
                }

            }
        }
    }

    return sampling_probs;
}