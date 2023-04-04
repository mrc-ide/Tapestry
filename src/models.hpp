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


// ================================================================================
// Interface for different models
//
// ================================================================================


class Model
{
protected:
    // All models will have parameters and data
    const Parameters& params;       // Model hyperparameters
    const VCFData& data;            // Data from VCF

public:
    // Constructors
    Model(const Parameters& params, const VCFData& data);


    /*
    * Calculate the log prior probability, 
    * give the parameters defined by `particle`
    */
    virtual double calc_logprior(const Particle& particle) const = 0;

    /*
    * Calculate loglikelihood of the model,
    * give the parameters defined by `particle`
    */
    virtual double calc_loglikelihood(const Particle& particle) const = 0;

    /*
    * Calculate log posterior probability,
    * give the parameters defined by `particle`
    */
    double calc_logposterior(const Particle& particle) const;

    virtual ~Model();
};



// ================================================================================
// Concrete model implementations
//
// ================================================================================

// --------------------------------------------------------------------------------
// No IBD model
// --------------------------------------------------------------------------------


class NoIBDModel : public Model
{
private:
    // MEMBERS
    const vector<int> strains;                    // Indices for within-sample strains
    const MatrixXi allele_configs;                // TODO: Here, needed as doubles
    const MatrixXd sampling_probs;
    const BetabinomialArray betabin_lookup;

    // FUNCTIONS
    vector<int> static create_strains(int K);
    MatrixXi static create_allele_configs(int K);
    MatrixXd static create_sampling_probs(
        const VCFData& data,
        const MatrixXi& allele_configs
    );

public:
    NoIBDModel(const Parameters& params, const VCFData& data);
    virtual double calc_logprior(const Particle& particle) const override;
    virtual double calc_loglikelihood(const Particle& particle) const override;
    void print() const;
};


// --------------------------------------------------------------------------------
// Naive IBD model
// --------------------------------------------------------------------------------


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


class NaiveIBDModel : public Model
{
private:
    // MEMBERS
    // Supporting arrays, computed upon initialisation
    const MatrixXi allele_configs;                // All possible allele configurations
    const IBDContainer ibd;                       // IBD state information
    const vector<MatrixXd> sampling_probs;        // Prob. of IBD/allele combo given PLAF
    const BetabinomialArray betabin_lookup;       // Precomputed WSAF ~ Betabin(...)
    const vector<MatrixXd> transition_matrices;   // Distance-dependent HMM trans. probs.
    // Pre-allocate memory for Forward Algorithm
    // TODO: this did not improve performance
    // MatrixXd F;             // Forward matrix
    // VectorXd scales;        // Scaling factors for forward matrix normalisation                                

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

    // TODO: Could split viterbi into methods below...
    // MatrixXi calc_traceback_matrix(const Particle& particle) const;
    // VectorXi run_traceback(const MatrixXi& tb) const;

public:
    NaiveIBDModel(const Parameters& params, const VCFData& data);
    virtual double calc_logprior(const Particle& particle) const override;
    virtual double calc_loglikelihood(const Particle& particle) const override;
    void print() const;

    // IBD path inference
    ViterbiResult get_viterbi_path(const Particle& particle) const;
};

