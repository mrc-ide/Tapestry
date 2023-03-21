#include <iostream>
#include <numeric>
#include <vector>
#include "betabin.hpp"
#include "combinatorics.hpp"
#include "data.hpp"
#include "models.hpp"
#include "parameters.hpp"
#include "particles.hpp"
#include "sampling.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::ArrayXi;


// ================================================================================
// Interface for different models
//
// ================================================================================


Model::Model(const Parameters& params, const VCFData& data)
    : params(params),
    data(data)
{};


double Model::calc_logposterior(const Particle& particle) const
{
    return calc_logprior(particle) + calc_loglikelihood(particle);
}


Model::~Model()
{};


// ================================================================================
// Concrete model implementations
//
// ================================================================================

// --------------------------------------------------------------------------------
// No IBD model
// --------------------------------------------------------------------------------


NoIBDModel::NoIBDModel(const Parameters& params, const VCFData& data)
    : Model(params, data),
    strains(create_strains(params.K)),
    allele_configs(create_allele_configs(params.K)),
    sampling_probs(create_sampling_probs(data, allele_configs)),
    betabin_lookup(params, data, false) // want linear Betabin probs
{};


std::vector<int> NoIBDModel::create_strains(int K)
{
    std::vector<int> strains(K);
    std::iota(strains.begin(), strains.end(), 0);
    return strains;
}


MatrixXi NoIBDModel::create_allele_configs(int K)
{
    return create_powerset(K);
}


MatrixXd NoIBDModel::create_sampling_probs(
        const VCFData& data,
        const MatrixXi& allele_configs
    )
{
    
    int n_strains = allele_configs.cols();  // elsewhere K
    int n_allele_configs = allele_configs.rows();
    
    // Prepare storage
    MatrixXd sampling_probs = MatrixXd::Constant(data.n_sites, n_allele_configs, -1.0);

    // Compute total ALT and REF carrying strains for each allele config
    ArrayXd alt_counts = allele_configs.rowwise().sum().cast<double>();
    ArrayXd ref_counts = ArrayXd::Constant(n_allele_configs, n_strains) - alt_counts;

    // Compute sampling probability p^ALT(1 - p)^REF
    // TODO: possibly vectorised implementation of this entire routine
    for (int i = 0; i < data.n_sites; ++i) {
        double p = data.plafs(i);
        sampling_probs.row(i) = Eigen::pow(p, alt_counts) * Eigen::pow(1-p, ref_counts);
    }

    // TODO:
    // - we want sites as columns for subsequent matrix multi;
    // - but looping over rows should be faster, hence transposing rather than building this way
    return sampling_probs.transpose();
}


double NoIBDModel::calc_logprior(const Particle& particle) const
{
    return 0;  // Uniform
}


double NoIBDModel::calc_loglikelihood(const Particle& particle) const
{
    ArrayXd wsaf = (allele_configs.cast<double>() * particle.ws.transpose()).array();
    ArrayXd wsaf_adj = (1 - wsaf) * params.e_0 + (1 - params.e_1) * wsaf;

    MatrixXd wsaf_betabin_probs = betabin_lookup.subset(wsaf_adj);

    VectorXd emission_probs(data.n_sites);
    for (int i = 0; i < data.n_sites; ++i) {
        emission_probs(i) = wsaf_betabin_probs.row(i) * sampling_probs.col(i);
    }
    
    return emission_probs.array().log().sum(); // TODO: base?
}


void NoIBDModel::print() const
{
    std::cout << "Strains" << std::endl;
    for (int i=0; i<params.K; ++i) std::cout << strains[i] << "\t";
    std::cout << endl;
    std::cout << "Allele Configurations:" << std::endl;
    std::cout << allele_configs << std::endl;
    std::cout << "Sampling probabilities: " << std::endl;
    std::cout << sampling_probs.col(0) << std::endl;
}


// --------------------------------------------------------------------------------
// Naive IBD model
// --------------------------------------------------------------------------------


NaiveIBDModel::NaiveIBDModel(const Parameters& params, const VCFData& data)
    : Model(params, data),
    strains(create_strains(params.K)),
    allele_configs(create_allele_configs(params.K)),
    ibd_states(create_ibd_states(strains)),
    sampling_probs(create_sampling_probs(data, allele_configs, ibd_states)),
    betabin_lookup(params, data)
{};


std::vector<int> NaiveIBDModel::create_strains(int K)
{
    std::vector<int> strains(K);
    std::iota(strains.begin(), strains.end(), 0);
    return strains;
}


MatrixXi NaiveIBDModel::create_allele_configs(int K)
{
    return create_powerset(K);
}


vector<vector<vector<int>>> NaiveIBDModel::create_ibd_states(vector<int> strains)
{
    return create_all_partitions(strains);
}

vector<MatrixXd> NaiveIBDModel::create_sampling_probs(
        const VCFData& data,
        const MatrixXi& allele_configs,
        const vector<vector<vector<int>>>& ibd_states
        )
{
    // Initialise
    vector<MatrixXd> sampling_probs(
        data.n_sites,
        MatrixXd::Constant(allele_configs.rows(), ibd_states.size(), -1.0)
    );

    for (int i = 0; i < data.n_sites; ++i) {
        // TODO: 
        // - This is a copy step, which is bad
        // - Better to use references to avoid
        sampling_probs[i] = calc_sampling_probs(
            data.plafs(i),
            allele_configs,
            ibd_states
        );
    }

    return sampling_probs;
}


void NaiveIBDModel::print() const
{
    std::cout << "Strains" << std::endl;
    for (int i=0; i<params.K; ++i) std::cout << strains[i] << "\t";
    std::cout << endl;
    std::cout << "Allele Configurations:" << std::endl;
    std::cout << allele_configs << std::endl;
    std::cout << "No. IBD States: " << ibd_states.size() << std::endl;
    std::cout << "Sampling probabilities: " << std::endl;
    std::cout << "First array: " << std::endl;
    std::cout << sampling_probs[0] << std::endl;
    std::cout << sampling_probs[data.n_sites - 1] << std::endl;
}


double NaiveIBDModel::calc_logprior(const Particle& particle) const
{
    return 0;  // Uniform
}


double NaiveIBDModel::calc_loglikelihood(const Particle& particle) const
{
    return 0; // TODO: Implement forward algorithm
}

