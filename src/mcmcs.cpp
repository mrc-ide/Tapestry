#include <algorithm>
#include <iostream>
#include <iterator>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <random>
#include <vector>
#include "mcmcs.hpp"
#include "models.hpp"
#include "proposals.hpp"
#include "particles.hpp"
#include "particle_writers.hpp"
#include "random.hpp"
namespace fs = std::filesystem;


// ================================================================================
// MCMC Interface
//
// ================================================================================


MCMC::MCMC(
    const Parameters& params,
    const Model& model,
    ProposalEngine& proposal_engine) 
    : params(params),
    model(model),
    proposal_engine(proposal_engine),
    ix(0),
    n_burn_iters(100),
    n_sample_iters(900),
    n_total_iters(n_burn_iters + n_sample_iters),
    acceptance_rate(-1.0),
    acceptance_trace(n_total_iters),
    logposterior_trace(n_total_iters),
    particle_trace(n_total_iters, params.K)  // TODO: Here is where space gets allocated. Double check.
{};


void MCMC::write_output(
    const string& output_dir, 
    const ParticleWriter& particle_writer) const
{
    
    // Create the directory if it doesn't exist
    fs::create_directories(output_dir);  // TODO: not sure this is cross-platform
    
    // Prepare file paths
    std::string mcmc_csv = output_dir + "/mcmc.diagnostics.csv";
    std::string particles_csv = output_dir + "/mcmc.parameters.csv";

    // Write MCMC diagnostics
    std::ofstream csv_file(mcmc_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }
    csv_file << "iter,phase,logposterior,acceptance_rate\n";
    for (int i = 0; i < ix; ++i) {
        csv_file << i << ",";
        csv_file << (i < n_burn_iters ? "burn" : "sample") << ","; 
        csv_file << logposterior_trace[i] << ",";
        csv_file << acceptance_trace[i] << "\n";
    }
    csv_file.close();

    // Write MCMC particles
    particle_writer.write_particle_trace(
        particles_csv,
        particle_trace
    );
}


Particle MCMC::get_map_particle() const
{
    // NB: we exclude the burn-in
    auto mapIter = std::max_element(
        logposterior_trace.begin() + n_burn_iters, 
        logposterior_trace.end()
    );
    int ix = std::distance(logposterior_trace.begin(), mapIter);
    
    return particle_trace[ix];
}


MCMC::~MCMC()
{}


// ================================================================================
// Concrete MCMC methods
//
// ================================================================================

// --------------------------------------------------------------------------------
// Metropolis-Hastings
// --------------------------------------------------------------------------------


MetropolisHastings::MetropolisHastings(
    const Parameters& params, 
    const Model& model, 
    ProposalEngine& proposal_engine)
    : MCMC(params, model, proposal_engine)
{};


void MetropolisHastings::run_burn()
{
    // Random initialisation
    ix = 0;
    particle_trace[ix] = proposal_engine.create_particle();
    logposterior_trace[ix] = model.calc_logposterior(particle_trace[ix]);
    acceptance_rate = 1.0;
    acceptance_trace[ix] = 1.0;
    ++ix;

    // Run burn-in iterations
    run_iterations(n_burn_iters - 1);
}


void MetropolisHastings::run_sampling()
{
    run_iterations(n_sample_iters);
}


void MetropolisHastings::run_iterations(int n)
{
    // Check that ix > 0
    if (ix == 0) {
        throw std::invalid_argument("Initialise the particles before iterating.");
    }

    // Current state
    Particle particle = particle_trace[ix - 1];
    double logposterior = logposterior_trace[ix - 1];
    
    // Iterate
    int N = ix + n;
    for (; ix < N; ++ix) {

        // Propose
        Particle proposed_particle = proposal_engine.propose_particle(particle);
        double proposed_logposterior = model.calc_logposterior(proposed_particle);

        // Compute acceptance probability
        // TODO: Important to verify in the same base, else biased
        double A = std::exp(proposed_logposterior - logposterior); 
        double u = U(rng.engine);

        // Accept
        if (u < A) {
            particle = proposed_particle; // TODO: No longer need proposed particle, should be a move operation
            logposterior = proposed_logposterior;
        }

        // Store
        particle_trace[ix] = particle;
        logposterior_trace[ix] = logposterior;

        // Track expected acceptance rate
        acceptance_rate += (A < 1.0 ? A : 1.0);
        acceptance_trace[ix] = acceptance_rate / ix;
    }
}


void MetropolisHastings::run()
{
    run_burn();
    run_sampling();
}

