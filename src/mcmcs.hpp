#pragma once
#include <random>
#include <string>
#include <vector>
#include "models.hpp"
#include "proposals.hpp"
#include "particles.hpp"
#include "particle_writers.hpp"
#include "random.hpp"


// ================================================================================
// Interface for different MCMC methods
//
// ================================================================================


class MCMC
{
protected:
    // MEMBERS
    // Randomness
    RNG rng;
    std::uniform_real_distribution<double> U{0.0, std::nextafter(1.0, 2.0)};
    
    // Model specific
    const Parameters& params;
    const Model& model;
    ProposalEngine& proposal_engine;

    // Parameters
    int ix;                                     // Iteration index
    const int n_burn_iters;                           // Number of burn-in iterations
    const int n_sample_iters;                         // Number of sampling iterations
    const int n_total_iters;
    double acceptance_rate;                     // Accept rate until `ix`
    
    // Storage
    // TODO: is there a reason not to use Eigen for acceptance / logposterior?
    // - Initialisation a  bit trickier
    // - Change to circular brackets
    // - Otherwise don't see issue
    std::vector<double> acceptance_trace;        // Rolling E[acceptance rate]
    std::vector<double> logposterior_trace;     // Trace of log-posterior values
    std::vector<Particle> particle_trace;       // Trace of particles (i.e. updates)

public:
    MCMC(
        const Parameters& params, 
        const Model& model, 
        ProposalEngine& proposal_engine
    );

    // Abstract
    virtual void run() = 0;

    void virtual write_output(
        const string& output_dir, 
        const ParticleWriter& particle_writer) const;

    Particle get_map_particle() const;

    virtual ~MCMC();
};


// ================================================================================
// Concrete MCMC methods
//
// ================================================================================

// --------------------------------------------------------------------------------
// Metropolis-Hastings
// --------------------------------------------------------------------------------


class MetropolisHastings : public MCMC
{
private:
    void run_iterations(int n);
    void run_burn();
    void run_sampling();

public:
    MetropolisHastings(
        const Parameters& params, 
        const Model& model, 
        ProposalEngine& proposal_engine
    );
    void run() override;
};



// --------------------------------------------------------------------------------
// Parallel Tempering
// --------------------------------------------------------------------------------


class ParallelTempering : public MCMC
{
private:
    /**
    * Define a temperature level
    */
    struct TemperatureLevel
    {
        double beta;
        Particle* particle_ptr;
        double loglikelihood;  // Needed for swaps; Thermodynamic Integration (TI)
        double beta_logposterior;   // Needed for within-level MH updates

        TemperatureLevel();    // TODO: Hmm.. what goes in constructor?
    };

    std::vector<ParallelTempering::TemperatureLevel>  static create_temp_levels(
        std::vector<Particle>& particles, 
        double lambda=0.5
    );

    void run_iterations(int n);
    void run_burn();
    void run_sampling();

public:
    // DATA MEMBERS
    const int n_temps;                    // Number of temperature levels
    const int swap_freq;                  // Number of iterations per swap attempt
    std::vector<Particle> particles;      // Particle for each temperature
    std::vector<TemperatureLevel> temps;  // Temperature information
    MatrixXd loglikelihoods;              // Loglikelihoods; for TI; TODO: should be array?

    //  Recording swaps
    double n_swap_attempts;               // No. swap attempts; same for all pairs of levels
    ArrayXd n_swaps;                      // No. swaps for each pair, up to current `ix`
    MatrixXd swap_rates;                  // Rate of swapping for each pair of temp. levels
    
    ParallelTempering(
        const Parameters& params, 
        const Model& model, 
        ProposalEngine& proposal_engine,
        int n_temps
    );

    void run() override;

    // Concrete
    void write_output(
        const std::string& output_dir, 
        const ParticleWriter& particle_writer) const override;
};

