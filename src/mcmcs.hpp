#pragma once
#include <random>
#include <string>
#include <vector>
#include "models.hpp"
#include "proposals.hpp"
#include "particles.hpp"
#include "particle_writers.hpp"
#include "random.hpp"


// --------------------------------------------------------------------------------
// Parallel Tempering MCMC
// --------------------------------------------------------------------------------

class MCMC
{
private:

    // MEMBERS

    struct TemperatureLevel
    {
        double beta;
        Particle* particle_ptr;
        double loglike;
        double logprior;
        double w_prop_sd;

        TemperatureLevel();
    };

    // Randomness
    RNG rng;
    std::uniform_real_distribution<double> U{0.0, std::nextafter(1.0, 2.0)};
    
    // Model specific
    const Parameters& params;
    const Model& model;
    ProposalEngine& proposal_engine;

    // FUNCTIONS
    std::vector<MCMC::TemperatureLevel>  static create_temp_levels(
        std::vector<Particle>& particles, 
        double beta_skew = 5.0
    );

    void run_iterations(int n, bool adaptive_on = false);
    void run_burn();
    void run_sampling();

public:
    // MEMBERS
    // Iterations
    int ix;                                     // Iteration index
    const int n_burn_iters;                     // Number of burn-in iterations
    const int n_sample_iters;                   // Number of sampling iterations
    const int n_total_iters;
    double acceptance_rate_cumul;               // Cumulative acceptance rate until `ix` (cold chain only)
    
    // Storage
    // TODO: is there a reason not to use Eigen for acceptance / logposterior?
    // - Initialisation a  bit trickier
    // - Change to circular brackets
    // - Otherwise don't see issue
    std::vector<double> acceptance_trace;        // Rolling E[acceptance rate]
    std::vector<double> loglike_trace;           // Trace of log-likelihood
    std::vector<double> logprior_trace;          // Trace of log-prior
    std::vector<Particle> particle_trace;       // Trace of particles (i.e. updates)

    const int n_temps;                    // Number of temperature levels
    const int swap_freq;                  // Number of iterations per swap attempt
    std::vector<Particle> particles;      // Particle for each temperature
    std::vector<TemperatureLevel> temps;  // Temperature information
    MatrixXd loglikelihoods;              // Loglikelihoods; for TI; TODO: should be array? Probably, see ModelEvidence 104-106

    //  Recording swaps
    int n_swap_attempts;                  // No. swap attempts; same for all pairs of rungs
    ArrayXd swap_rate_cumul;              // Cumulative swap rate for each pair, up to current `ix`
    MatrixXd swap_rates;                  // Rate of swapping for each pair of rungs
    

    // FUNCTIONS
    // Constructor
    MCMC(
        const Parameters& params, 
        const Model& model, 
        ProposalEngine& proposal_engine,
        int n_burn_iters,
        int n_sample_iters,
        int n_temps
    );

    void run();

    void write_output(
        const string& output_dir, 
        const ParticleWriter& particle_writer) const;

    Particle get_map_particle() const;

    // Destructor
    ~MCMC() {}
};

