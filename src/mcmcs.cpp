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


// --------------------------------------------------------------------------------
// Parallel Tempering
// --------------------------------------------------------------------------------


ParallelTempering::ParallelTempering(
    const Parameters& params, 
    const Model& model, 
    ProposalEngine& proposal_engine,
    int n_temps)
    : MCMC(params, model, proposal_engine),
    n_temps(n_temps),
    swap_freq(10),  // For now just set as constant
    particles(n_temps),  // Not exactly sure how this looks
    temps(create_temp_levels(particles)),
    loglikelihoods(MatrixXd::Constant(n_total_iters, n_temps, -9999)),
    n_swap_attempts(0.0),
    n_swaps(ArrayXd::Constant(n_temps - 1, 0.0)),
    swap_rates(MatrixXd::Constant(n_total_iters, n_temps - 1, 0.0))
{}



ParallelTempering::TemperatureLevel::TemperatureLevel()
    : beta(0.0),
    particle_ptr(NULL),
    loglikelihood(-9999),
    beta_logposterior(-9999)
{}



std::vector<ParallelTempering::TemperatureLevel> ParallelTempering::create_temp_levels(
    std::vector<Particle>& particles, 
    double lambda)
{

    // Initialise
    int n_temps = particles.size();
    std::vector<ParallelTempering::TemperatureLevel> temp_levels(n_temps);

    // Populate
    double beta = 1.0;
    for (int j = n_temps - 1; j >= 0; --j) {
        temp_levels[j].particle_ptr = &particles[j];
        temp_levels[j].beta = beta;
        beta *= lambda;
    }

    return temp_levels;
}


void ParallelTempering::run_burn()
{
    // Initialise
    for (int j = 0; j < n_temps; ++j) {
        particles[j] = proposal_engine.create_particle();
        temps[j].loglikelihood = model.calc_loglikelihood(particles[j]);
        temps[j].beta_logposterior = temps[j].beta * temps[j].loglikelihood + model.calc_logprior(particles[j]);
    }
    // Store cold chain
    particle_trace[ix] = *temps[n_temps - 1].particle_ptr;
    logposterior_trace[ix] = temps[n_temps - 1].beta_logposterior;
    acceptance_rate = 1.0;
    acceptance_trace[ix] = acceptance_rate;
    ++ix;

    run_iterations(n_burn_iters - 1);
}


void ParallelTempering::run_sampling()
{
    run_iterations(n_sample_iters);
}



void ParallelTempering::run_iterations(int n)
{
    int N = ix + n;
    for (; ix < N; ++ix) {

        // Within temperature-level MH
        for (int j = 0; j < n_temps; ++j) {

            // Propose a particle
            TemperatureLevel& temp_level = temps[j];
            Particle proposed_particle = proposal_engine.propose_particle(*temp_level.particle_ptr);  // TODO: does this work?

            // Compute proposed likelihood and posterior
            double proposed_loglikelihood = model.calc_loglikelihood(proposed_particle);
            double proposed_logposterior = temp_level.beta * proposed_loglikelihood + model.calc_logprior(proposed_particle);

            // Compute acceptance rate
            double A = proposed_logposterior - temp_level.beta_logposterior; // TODO: Not that this way, it is not E[acceptance]
            double u = std::log(U(rng.engine));

            // Accept
            if (u < A) {
                *temp_level.particle_ptr = proposed_particle;  // change the particle's value
                temp_level.loglikelihood = proposed_loglikelihood;
                temp_level.beta_logposterior = proposed_logposterior;
                ++acceptance_rate;
            }

            // Store the loglikelihood for TI
            loglikelihoods(ix, j) = temp_level.loglikelihood;
        }

        // Between-temperature swaps
        if (ix % swap_freq == 0) {
            for (int j = 1; j < n_temps; ++j) {
                double a1 = (temps[j].beta - temps[j-1].beta) * temps[j-1].loglikelihood;
                double a2 = (temps[j-1].beta - temps[j].beta) * temps[j].loglikelihood;
                double A = a1 + a2;
                double u = std::log(U(rng.engine));
                
                // Accept
                if (u < A) {
                    std::swap(temps[j].particle_ptr, temps[j-1].particle_ptr);
                    std::swap(temps[j].loglikelihood, temps[j-1].loglikelihood);

                    // Recompute
                    temps[j].beta_logposterior = temps[j].beta * temps[j].loglikelihood + model.calc_logprior(*temps[j].particle_ptr);
                    temps[j-1].beta_logposterior = temps[j-1].beta * temps[j-1].loglikelihood + model.calc_logprior(*temps[j-1].particle_ptr);

                    // Increment if swapped
                    ++n_swaps(j-1);
                }
            }
            ++n_swap_attempts;
        }

        // Compute running swap rates
        swap_rates.row(ix) = n_swaps / n_swap_attempts;

        // Store cold chain for this iteration
        particle_trace[ix] = *temps[n_temps - 1].particle_ptr;
        logposterior_trace[ix] = temps[n_temps - 1].beta_logposterior;
        acceptance_trace[ix] = acceptance_rate / (n_temps*ix); // TODO: just average over all chains, should be cold!
    }
}


void ParallelTempering::run()
{
    run_burn();
    run_sampling();
}


void ParallelTempering::write_output(
    const std::string& output_dir,
    const ParticleWriter& particle_writer) const
{
    
    // Call parent method
    MCMC::write_output(output_dir, particle_writer);

    // Prepare file paths
    std::string llk_csv = output_dir + "/mcmc.likelihood.csv";
    std::string beta_csv = output_dir + "/mcmc.betas.csv";
    std::string swap_csv = output_dir + "/mcmc.swap_rates.csv";

    // Write loglikelihood matrix
    const static Eigen::IOFormat CSVFormat(6, Eigen::DontAlignCols, ",",  "\n");
    std::ofstream llk_file(llk_csv);
    if (!llk_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }
    for (int j = 0; j < n_temps; ++j) {
        llk_file << "level" << j << (j < n_temps - 1 ? "," : "\n");
    }
    llk_file << loglikelihoods.format(CSVFormat);
    llk_file.close();

    // Write beta values
    std::ofstream beta_file(beta_csv);
    if (!beta_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }
    beta_file << "level,beta\n";
    for (int j = 0; j < n_temps; ++j) {
        beta_file << j << "," << temps[j].beta << "\n";
    }
    beta_file.close();

    // Quick write of swap rates; should add to mcmc.diagnostics.csv
    std::ofstream swap_file(swap_csv);
    if (!swap_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }
    for (int j = 0; j < n_temps - 1; ++j) {
        swap_file << "swap_" << j << "-" << j+1 << (j < n_temps - 2 ? "," : "\n");
    }
    swap_file << swap_rates.format(CSVFormat);
    swap_file.close();
}

