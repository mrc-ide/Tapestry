#include <vector>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include "typedefs.hpp"
#include "sampling.hpp"
#include "constants.hpp"
using namespace std;



MatrixXd calc_sampling_probs(
    const double p,
    const MatrixXi& allele_configs,
    const vector<vector<vector<int>>>& ibd_states
)
{
    int n_ibd_states = ibd_states.size();
    int n_allele_configs = allele_configs.rows();
    int k = allele_configs.cols(); // COI is implicity in shape of `allele_configs`

    if (p < 0 || p > 1) {
        throw std::invalid_argument("PLAF `p` must be in [0, 1].");
    }
    if (n_ibd_states != BELL_NUMBERS[k]) {
        throw std::invalid_argument("Allele configurations and IBD states do not correspond to same COI.");
    }

    // Initiate
    MatrixXd sampling_probs = MatrixXd::Constant(n_allele_configs, n_ibd_states, 1.0);

    // Populate
    for (int i = 0; i < n_allele_configs; ++i) {
        for (int j = 0; j < n_ibd_states; ++j) {

            // Collect IBD state
            const vector<vector<int>>& ibd_state = ibd_states[j];

            // Iterate over IBD groups within state
            for (const vector<int>& ibd : ibd_state) {

                // Compute number of strains and ALT alleles in IBD group  
                int n_strains = ibd.size();
                int n_alt = 0;
                for (const int& ix : ibd) {
                    n_alt += allele_configs(i, ix);
                }

                // Compute sampling probability of this IBD group
                if (n_strains == n_alt) {
                    sampling_probs(i, j) *= p;
                } else if (n_alt == 0) {
                    sampling_probs(i, j) *= (1 - p);
                } else {
                    sampling_probs(i, j) = 0;
                    break; // IBD and haplotype configs incompatible
                }
            }
        }
    }

    return sampling_probs;

}

