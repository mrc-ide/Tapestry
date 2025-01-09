#pragma once
#include "betabin.hpp"
#include "combinatorics.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "parameters.hpp"
#include "particles.hpp"
#include "sampling.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;


/*
* Capture output of the viterbi algorithm
* TODO: where to put this?
*/
struct ViterbiResult
{
    VectorXi path;
    double logposterior;
    ViterbiResult() {}  // TODO: is this really necessary?
    ViterbiResult(int n_sites)
    : path(VectorXi::Constant(n_sites, 9999)),
    logposterior(0)
    {}
};

// --------------------------------------------------------------------------------
// Naive IBD model
// --------------------------------------------------------------------------------

class Model
{
private:

    // MEMBERS
    // Parameters and data
    const Parameters& params;       // Model hyperparameters
    const VCFData& data;            // Data from VCF

    // Supporting arrays, computed upon initialisation
    const MatrixXi allele_configs;                // All possible allele configurations
    const IBDContainer ibd;                       // IBD state information
    const vector<MatrixXd> sampling_probs;        // Prob. of IBD/allele combo given PLAF
    const BetabinomialArray betabin_lookup;       // Precomputed WSAF ~ Betabin(...)
    const vector<MatrixXd> transition_matrices;   // Distance-dependent HMM trans. probs.

    // FUNCTIONS
    MatrixXi static create_allele_configs(int K);
    vector<MatrixXd> static create_sampling_probs(
        const VCFData& data,
        const MatrixXi& allele_configs,
        const vector<vector<vector<int>>>& ibd_states
    );
    MatrixXd static calc_transition_matrix(int d_ij, const Parameters& params);
    vector<MatrixXd> static create_transition_matrices( 
        const Parameters& params,
        const VCFData& data
    );

public:

    // Constructors
    Model(const Parameters& params, const VCFData& data);

    /*
    * Calculate the log prior probability, 
    * give the parameters defined by `particle`
    */
    double calc_logprior(const Particle& particle) const;

    /*
    * Calculate loglikelihood of the model,
    * give the parameters defined by `particle`
    */
    double calc_loglikelihood(const Particle& particle) const;

    /*
    * Calculate log posterior probability,
    * give the parameters defined by `particle`
    */
    double calc_logposterior(const Particle& particle) const;

    // IBD path inference
    ViterbiResult get_viterbi_path(const Particle& particle) const;

    // Print method
    void print() const;

    // Destructor
    ~Model() {}
};

