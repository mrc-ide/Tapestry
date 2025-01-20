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


// --------------------------------------------------------------------------------
// Parallel Tempering
// --------------------------------------------------------------------------------

MCMC::MCMC(
    const Parameters& params,
    const Model& model,
    ProposalEngine& proposal_engine,
    int n_temps) 
    : params(params),
    model(model),
    proposal_engine(proposal_engine),
    ix(0),
    n_burn_iters(100),
    n_sample_iters(900),
    n_total_iters(n_burn_iters + n_sample_iters),
    acceptance_rate_cumul(0.0),
    acceptance_trace(n_total_iters),
    loglike_trace(n_total_iters),
    logprior_trace(n_total_iters),
    particle_trace(n_total_iters, params.K),  // TODO: Here is where space gets allocated. Double check.
    n_temps(n_temps),
    swap_freq(10),  // For now just set as constant
    particles(n_temps),  // Not exactly sure how this looks
    temps(create_temp_levels(particles)),
    loglikelihoods(MatrixXd::Constant(n_total_iters, n_temps, -9999)),
    n_swap_attempts(0.0),
    n_swaps(ArrayXd::Constant(n_temps - 1, 0.0)),
    swap_rates(MatrixXd::Constant(n_total_iters, n_temps - 1, 0.0))
{};


void MCMC::write_output(
    const string& output_dir, 
    const ParticleWriter& particle_writer) const
{
    
    // Create the directory if it doesn't exist
    fs::create_directories(output_dir);  // TODO: not sure this is cross-platform
    
    // Prepare file paths
    std::string mcmc_csv = output_dir + "/mcmc.trace.csv";
    std::string particles_csv = output_dir + "/mcmc.parameters.csv";

    // Write MCMC diagnostics
    std::ofstream csv_file(mcmc_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }
    csv_file << "iter,phase,loglike,logprior,acceptance_rate\n";
    for (int i = 0; i < ix; ++i) {
        csv_file << i << ",";
        csv_file << (i < n_burn_iters ? "burn" : "sample") << ","; 
        csv_file << loglike_trace[i] << ",";
        csv_file << logprior_trace[i] << ",";
        csv_file << acceptance_trace[i] << "\n";
    }
    csv_file.close();

    // Write MCMC particles
    particle_writer.write_particle_trace(
        particles_csv,
        particle_trace
    );

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

    // Quick write of swap rates; should add to mcmc.trace.csv
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


Particle MCMC::get_map_particle() const
{
    // Find iteration that achieved the highest joint probability
    int i_max = 0;
    double logjoint_max = loglike_trace[0] + logprior_trace[0];
    for (int i = 0; i < loglike_trace.size(); ++i) {
        double logjoint = loglike_trace[i] + logprior_trace[i];
        if (logjoint > logjoint_max) {
            i_max = i;
            logjoint_max = logjoint;
        }
    }
    Particle map_particle = particle_trace[i_max];
    std::sort(map_particle.ws.begin(), map_particle.ws.end());
    return map_particle;
}


MCMC::TemperatureLevel::TemperatureLevel()
    : beta(0.0),
    particle_ptr(NULL),
    loglike(-9999),
    logprior(-9999)
{}


std::vector<MCMC::TemperatureLevel> MCMC::create_temp_levels(
    std::vector<Particle>& particles, 
    double beta_skew)
{

    // Initialise
    int n_temps = particles.size();
    std::vector<MCMC::TemperatureLevel> temp_levels(n_temps);

    // Populate with a sequence from 0 to 1, raised to the power beta_skew
    for (int j = 0; j < n_temps; ++j) {
        temp_levels[j].particle_ptr = &particles[j];
        temp_levels[j].beta = pow(j / double(n_temps - 1), beta_skew);
    }

    return temp_levels;
}


void MCMC::run_burn()
{
    // Initialise
    for (int j = 0; j < n_temps; ++j) {
        particles[j] = proposal_engine.create_particle();
        temps[j].loglike = model.calc_loglikelihood(particles[j]);
        temps[j].logprior = model.calc_logprior(particles[j]);
    }
    // Store cold chain
    particle_trace[ix] = *temps[n_temps - 1].particle_ptr;
    loglike_trace[ix] = temps[n_temps - 1].loglike;
    logprior_trace[ix] = temps[n_temps - 1].logprior;
    acceptance_rate_cumul = 1.0;
    acceptance_trace[ix] = 1.0;
    ++ix;

    run_iterations(n_burn_iters - 1);
}


void MCMC::run_sampling()
{
    run_iterations(n_sample_iters);
}



void MCMC::run_iterations(int n)
{
    int N = ix + n;
    for (; ix < N; ++ix) {

        // Within temperature-level MH
        for (int j = 0; j < n_temps; ++j) {

            // Propose a particle
            TemperatureLevel& temp_level = temps[j];
            Particle proposed_particle = proposal_engine.propose_particle(*temp_level.particle_ptr, params.w_proposal_sd);

            // Compute proposed likelihood and prior
            double proposed_loglike = model.calc_loglikelihood(proposed_particle);
            double proposed_logprior = model.calc_logprior(proposed_particle);

            // Compute acceptance rate under power-posterior
            double A = temp_level.beta * (proposed_loglike - temp_level.loglike) + (proposed_logprior - temp_level.logprior);
            double u = std::log(U(rng.engine));

            // Accept
            if (u < A) {
                *temp_level.particle_ptr = proposed_particle;  // change the particle's value
                temp_level.loglike = proposed_loglike;
                temp_level.logprior = proposed_logprior;
            }

            // Store the loglikelihood for TI
            loglikelihoods(ix, j) = temp_level.loglike;

            // Track expected acceptance rate in cold chain
            if (j == (n_temps - 1)) {
                acceptance_rate_cumul += (A < 0.0 ? std::exp(A) : 1.0);
            }
        }

        // Between-temperature swaps
        if (ix % swap_freq == 0) {
            for (int j = 1; j < n_temps; ++j) {
                double a1 = (temps[j].beta - temps[j-1].beta) * temps[j-1].loglike;
                double a2 = (temps[j-1].beta - temps[j].beta) * temps[j].loglike;
                double A = a1 + a2;
                double u = std::log(U(rng.engine));
                
                // Accept
                if (u < A) {
                    std::swap(temps[j].particle_ptr, temps[j-1].particle_ptr);
                    std::swap(temps[j].loglike, temps[j-1].loglike);
                    std::swap(temps[j].logprior, temps[j-1].logprior);

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
        loglike_trace[ix] = temps[n_temps - 1].loglike;
        logprior_trace[ix] = temps[n_temps - 1].logprior;
        acceptance_trace[ix] = acceptance_rate_cumul / double(ix + 1);
    }
}


void MCMC::run()
{
    run_burn();
    run_sampling();
}



