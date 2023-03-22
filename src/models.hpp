#pragma once
#include "betabin.hpp"
#include "combinatorics.hpp"
#include "data.hpp"
#include "parameters.hpp"
#include "particles.hpp"
#include "sampling.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;



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


class NaiveIBDModel : public Model
{
private:
    // MEMBERS
    // Supporting arrays, computed upon initialisation
    const vector<int> strains;                    // Indices for within-sample strains
    const MatrixXi allele_configs;                // All possible allele configurations
    const vector<vector<vector<int>>> ibd_states; // All possible IBD states
    const vector<MatrixXd> sampling_probs;        // Prob. of IBD/allele combo given PLAF
    const BetabinomialArray betabin_lookup;
    const vector<MatrixXd> transition_matrices;

    // FUNCTIONS
    vector<int> static create_strains(int K);
    MatrixXi static create_allele_configs(int K);
    vector<vector<vector<int>>> static create_ibd_states(vector<int> strains);
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
    NaiveIBDModel(const Parameters& params, const VCFData& data);
    virtual double calc_logprior(const Particle& particle) const override;
    virtual double calc_loglikelihood(const Particle& particle) const override;
    void print() const;
};

