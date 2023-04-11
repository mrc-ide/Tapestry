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
    int n_burn_iters;                           // Number of burn-in iterations
    int n_sample_iters;                         // Number of sampling iterations
    int n_total_iters;
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

    // Concrete
    void write_output(
        const string& output_dir, 
        const ParticleWriter& particle_writer) const;

    Particle get_map_particle() const;
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

