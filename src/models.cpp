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
using Eigen::VectorXi;
using Eigen::RowVectorXd;
using Eigen::ArrayXi;
using Eigen::ArrayXd;


// --------------------------------------------------------------------------------
// Naive IBD model
// --------------------------------------------------------------------------------

// Constructor
Model::Model(const Parameters& params, const VCFData& data)
    : params(params),
    data(data),
    allele_configs(create_allele_configs(params.K)),
    ibd(params.K),
    sampling_probs(create_sampling_probs(data.plafs, allele_configs, ibd.states)),
    betabin_lookup(data, params.n_pi_bins, params.e_0, params.e_1, params.v, false),
    transition_matrices(create_transition_matrices(params, data))
{};

Model::Model(const Parameters& params, const VCFData& data, const Betabinomial::LookupMatrix& betabin_lookup)
    : params(params),
    data(data),
    allele_configs(create_allele_configs(params.K)),
    ibd(params.K),
    sampling_probs(create_sampling_probs(data.plafs, allele_configs, ibd.states)),
    betabin_lookup(betabin_lookup),
    transition_matrices(create_transition_matrices(params, data))
{};

MatrixXi Model::create_allele_configs(int K)
{
    return create_powerset(K);
}


MatrixXd Model::calc_transition_matrix(int d_ij, const Parameters& params)
{
    
    // Co-efficients
    // can't declare these `static`! params changes during program execution!
    double bp_per_M = params.rho * 100 * 1000; // convert from kbp per cM
    double lambda = (params.G / bp_per_M) * params.K * (params.K - 1) / 2;

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

vector<MatrixXd> Model::create_transition_matrices(
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

    // TODO: Chromosome transitions
    // - Probably should also depend on distance to next pos;
    // - Two processes ocurring; segregeation and recombination
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
        // TODO: Is this a bug? Verify correct for entering new chromosomes
        d_ij = data.pos(i+1) - data.pos(i);
        transition_matrices[i] = calc_transition_matrix(d_ij, params); // TODO: this is a copy, should use reference
    }

    return transition_matrices;
}


double Model::calc_logprior(const Particle& particle) const
{
    return  lgamma(params.K);  // Uniform over unit simplex
}


double Model::calc_loglikelihood(const Particle& particle) const
{
    // Get adjusted WSAF values based on proportions, error parameters
    // TODO: shouldn't allele configs just be a double?
    // TODO: is there EVER benefit of particle being a row vector?
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
    int t = 0;
    F.row(t) = RowVectorXd::Constant(BELL_NUMBERS[params.K], 1.0/BELL_NUMBERS[params.K]);   // Initiation
    F.row(t).array() *= (wsaf_betabin_probs.row(t) * sampling_probs[t]).array(); // Emission
    scales(t) = F.row(t).sum();
    F.row(t) /= scales(t);
    loglike += log(scales(t));
    ++t;

    // Iterate
    for (; t < data.n_sites; ++t) {
        F.row(t) = (F.row(t-1) * transition_matrices[t-1]);   // Transition
        F.row(t).array() *= (wsaf_betabin_probs.row(t) * sampling_probs[t]).array(); // Emission
        scales(t) = F.row(t).sum();
        F.row(t) /= scales(t);
        loglike += log(scales(t));
    }

    return loglike;
}


// TODO: would be quite benefical to cache the logposterior for already used proportions
double Model::calc_logposterior(const Particle& particle) const
{
    return calc_logprior(particle) + calc_loglikelihood(particle);
}


// Compute Viterbi path of IBD given proportions in a Particle
ViterbiResult Model::get_viterbi_path(const Particle& particle) const
{
    // Compute error adjusted WSAF
    ArrayXd wsaf = (allele_configs.cast<double>() * particle.ws.transpose()).array();
    ArrayXd wsaf_adj = (1 - wsaf) * params.e_0 + (1 - params.e_1) * wsaf;

    // Extract relevant columns of Betabinomial array
    MatrixXd wsaf_betabin_probs = betabin_lookup.subset(wsaf_adj);

    // Prepare storage
    MatrixXi TB = MatrixXi::Constant(data.n_sites, BELL_NUMBERS[params.K], -1.0);
    // We never to matrix algebra with these, keep as arrays
    ArrayXd V = ArrayXd::Constant(BELL_NUMBERS[params.K], 9999.0); // V(i)_{t}
    ArrayXd Vt = ArrayXd::Constant(BELL_NUMBERS[params.K], 9999.0); // V(i)_{t+1}

    // Initialise
    int t = 0;
    // NB: Very important to do double division 1.0 / BELL ; NOT 1 / BELL
    V = ArrayXd::Constant(BELL_NUMBERS[params.K], log(1.0/BELL_NUMBERS[params.K])); 
    V += log((wsaf_betabin_probs.row(t) * sampling_probs[t]).array());
    ++t;
    

    // Iterate
    ptrdiff_t i;  // argmax index
    for (; t < data.n_sites; ++t) {
        for (int j = 0; j < BELL_NUMBERS[params.K]; ++j) {
            double maxV = (V + transition_matrices[t-1].col(j).array().log()).maxCoeff(&i);
            Vt(j) = maxV + log(wsaf_betabin_probs.row(t) * sampling_probs[t].col(j));
            TB(t, j) = i;
        }
        V = Vt;
    }

    // Terminate
    ViterbiResult viterbi(data.n_sites);
    viterbi.logposterior = Vt.maxCoeff(&i);  // for now, not using this
    int last_state = i;

    // Traceback
    t = data.n_sites - 1;  // should already have this value, if in same function
    viterbi.path(t) = last_state;

    // Recur
    while (t  > 0) {
        viterbi.path(t-1) = TB(t, viterbi.path(t));
        --t;
    }

    return viterbi;
}


void Model::print() const
{
    std::cout << "Strains" << std::endl;
    for (int i=0; i<params.K; ++i) std::cout << ibd.strains[i] << "\t";
    std::cout << endl;
    std::cout << "Allele Configurations:" << std::endl;
    std::cout << allele_configs << std::endl;
    std::cout << "No. IBD States: " << ibd.states.size() << std::endl;
    std::cout << "Sampling probabilities: " << std::endl;
    std::cout << "First array: " << std::endl;
    std::cout << sampling_probs[0] << std::endl;
    std::cout << sampling_probs[data.n_sites - 1] << std::endl;
    std::cout << "Transition matrices: " << std::endl;
    std::cout << transition_matrices[0] << std::endl;
    std::cout << transition_matrices[data.n_sites - 2] << std::endl;
}
