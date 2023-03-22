#include <iostream>
#include <numeric>
#include <vector>
#include "betabin.hpp"
#include "combinatorics.hpp"
#include "constants.hpp"
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
    betabin_lookup(params, data, false),  // want linear probabilities
    transition_matrices(create_transition_matrices(params, data))
    // F(MatrixXd::Constant(data.n_sites, BELL_NUMBERS[params.K], -1.0)),
    // scales(VectorXd::Constant(data.n_sites, -1.0))
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

MatrixXd NaiveIBDModel::calc_transition_matrix(int d_ij, const Parameters& params)
{
    
    // For now, fix
    int G = 20;

    // Co-efficient
    static double bp_per_M = params.rho * 100 * 1000; // convert from kbp per cM
    static double lambda = (G / bp_per_M) * params.K * (params.K - 1) / 2;

    // Compute possible matrix values
    double stay_prob = exp(-d_ij * lambda);
    double transition_prob = (1 - stay_prob)/(BELL_NUMBERS[params.K] - 1);

    MatrixXd tran_matrix = MatrixXd::Constant(
        BELL_NUMBERS[params.K], 
        BELL_NUMBERS[params.K],
        transition_prob
    );
    tran_matrix.diagonal() = VectorXd::Constant(BELL_NUMBERS[params.K], stay_prob);

    return tran_matrix;
}

vector<MatrixXd> NaiveIBDModel::create_transition_matrices(
    const Parameters& params, 
    const VCFData& data
    )
{
    // Initialise
    vector<MatrixXd> transition_matrices(
        data.n_sites - 1,
        MatrixXd::Constant(
            BELL_NUMBERS[params.K], 
            BELL_NUMBERS[params.K], 
            -1.0
        )
    );

    // Chromosome transitions
    MatrixXd chrom_transition_matrix = MatrixXd::Constant(
            BELL_NUMBERS[params.K], 
            BELL_NUMBERS[params.K], 
            1.0 / BELL_NUMBERS[params.K]  // All transitions equal
    );

    // TODO:
    // - the first matrix is the initialisation probabilities
    // - how to define here?
    // - it is really just a row vector
    // - need to get *indexing* correct for transitions
    // - think about forward algorithm
    // - at i, have probs i to i+1
    // - the last matrix ends up being blank

    // We need distances *between* positions
    // We need to handle chromosome ends
    int d_ij;
    for (int i = 0; i < data.n_sites - 1; ++i) {

        // Check if moving to new chromosome
        if (data.chroms(i) != data.chroms(i+1)) {
            transition_matrices[i] = chrom_transition_matrix;
            continue;
        }

        // Otherwise, compute distannce
        d_ij = data.pos(i+1) - data.pos(i);
        transition_matrices[i] = calc_transition_matrix(d_ij, params); // TODO: this is a copy, should use reference
    }

    return transition_matrices;
}


double NaiveIBDModel::calc_logprior(const Particle& particle) const
{
    return 0;  // Uniform
}


double NaiveIBDModel::calc_loglikelihood(const Particle& particle) const
{
    // Get adjusted WSAF values based on proportions, error parameters
    ArrayXd wsaf = (allele_configs.cast<double>() * particle.ws.transpose()).array();
    ArrayXd wsaf_adj = (1 - wsaf) * params.e_0 + (1 - params.e_1) * wsaf;

    // Extract relevant columns of Betabinomial array
    MatrixXd wsaf_betabin_probs = betabin_lookup.subset(wsaf_adj);

    // Create storage; have no made these class members to allocate once
    // This was just as fast as make an instance variable..
    MatrixXd F = MatrixXd::Constant(data.n_sites, BELL_NUMBERS[params.K], -1.0);
    VectorXd scales = VectorXd::Constant(data.n_sites, -1.0);

    // Initialise
    double loglike =  0;
    int i = 0;
    F.row(i) = RowVectorXd::Constant(BELL_NUMBERS[params.K], 1.0/BELL_NUMBERS[params.K]);   // Initiation
    F.row(i).array() *= (wsaf_betabin_probs.row(i) * sampling_probs[i]).array(); // Emission
    scales(i) = F.row(i).sum();
    F.row(i) /= scales(i);
    loglike += log(scales(i));
    ++i;

    // Iterate
    for (; i < data.n_sites; ++i) {
        F.row(i) = F.row(i-1) * transition_matrices[i-1];           // Transition
        F.row(i).array() *= (wsaf_betabin_probs.row(i) * sampling_probs[i]).array(); // Emission
        scales(i) = F.row(i).sum();
        F.row(i) /= scales(i);
        loglike += log(scales(i));
    }

    return loglike;
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
    std::cout << "Transition matrices: " << std::endl;
    std::cout << transition_matrices[0] << std::endl;
    std::cout << transition_matrices[data.n_sites - 2] << std::endl;
}
