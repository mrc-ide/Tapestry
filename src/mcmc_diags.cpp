#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include "mcmc_diags.hpp"

// --------------------------------------------------------------------------------
// MCMC diagnostics
// --------------------------------------------------------------------------------

// Compute mean of a vector
double MCMC_diags::mean(const std::vector<double>& x) {
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    return sum / x.size();
}

// Compute variance of a vector
double MCMC_diags::variance(const std::vector<double>& x) {
    double mean_val = mean(x);
    double sum = 0.0;
    for (double xi : x) {
        sum += (xi - mean_val) * (xi - mean_val);
    }
    return sum / (x.size() - 1);  // Sample variance
}

// Compute autocorrelation at a given lag
double MCMC_diags::autocorrelation(const std::vector<double>& x, int lag) {

    double x_mean = mean(x);
    double num = 0.0, denom = 0.0;

    for (size_t i = 0; i < x.size(); i++) {
        denom += (x[i] - x_mean) * (x[i] - x_mean);
        if (i + lag < x.size()) {
            num += (x[i] - x_mean) * (x[i + lag] - x_mean);
        }
    }

    return num / denom;
}

// Compute ESS using Geyer's initial positive sequence rule
double MCMC_diags::compute_ESS_geyer(const std::vector<double>& x)
{
    size_t N = x.size();
    if (N <= 1) {
        return 1.0;  // Avoid division by zero
    }

    double rho_sum = 0.0;
    std::vector<double> rho_values;

    // Compute autocorrelations and apply Geyer’s truncation
    for (int lag = 1; lag < (int)N / 2; lag++) {  // Limit to N/2 to be conservative
        double rho = autocorrelation(x, lag);
        rho_values.push_back(rho);

        // Stop summing when the sum of two consecutive autocorrelations turns negative
        if (lag > 1 && (rho_values[lag - 2] + rho_values[lag - 1] < 0)) {
            break;
        }

        rho_sum += rho;
    }

    return N / (1.0 + 2 * rho_sum);
}

// Geweke's diagnostic function. Ess_scale accounts for autocorrelation in the chain when estimating the variance of the statistic
double MCMC_diags::geweke_diagnostic(const std::vector<double>& x, double ESS_scale, double fracA, double fracB) {
    int N = x.size();
    if (N < 2) {
        return 0.0;  // Avoid division by zero for small samples
    }

    // Define segment sizes
    int nA = static_cast<int>(fracA * N);
    int nB = static_cast<int>(fracB * N);

    // Ensure nA and nB are valid
    if (nA == 0 || nB == 0 || nA + nB > N) {
        return 0.0;
    }

    // Extract first and last parts of the vector
    std::vector<double> xA(x.begin(), x.begin() + nA);
    std::vector<double> xB(x.end() - nB, x.end());

    // Compute means and variances
    double meanA = mean(xA);
    double meanB = mean(xB);
    double varA = variance(xA);
    double varB = variance(xB);

    // Compute Geweke's Z-score
    double Z = (meanA - meanB) / std::sqrt((varA / (nA * ESS_scale)) + (varB / (nB * ESS_scale)));
    return Z;
}

// Extract a single parameter from the vector of particles
std::vector<double> MCMC_diags::extract_parameter(const std::vector<Particle>& particle_trace, std::function<double(const Particle&)> extractor) {

    std::vector<double> ret;
    for (const auto& particle : particle_trace) {
        ret.push_back(extractor(particle));
    }
    return ret;
}

// Calculate diagnostics for all parameters
void MCMC_diags::calc_diagnostics(const std::vector<Particle>& particle_trace, const int n_burn_iters)
{
    // Initialise member objects
    int k = particle_trace[0].ws.size();
    Geweke_w.resize(k);
    convergence_pass = std::vector<bool>(k, false);
    ess_w.resize(k);

    // Create trimmed version that is sampling phase only
    std::vector<Particle> particle_trace_sampling(particle_trace.begin() + n_burn_iters, particle_trace.end());

    // Compute ESS for all free parameters
    for (int i = 0; i < k; ++i) {

        // Extract parameter vector and subset to sampling phase
        std::vector<double> param_vec = extract_parameter(particle_trace, [i](const Particle& p) { return p.ws[i];} );
        std::vector<double> param_vec_sampling(param_vec.begin() + n_burn_iters, param_vec.end());

        // Get ESS
        ess_w[i] = compute_ESS_geyer(param_vec_sampling);

        // Get stationarity (convergence)
        double ESS_scale = ess_w[i] / param_vec_sampling.size();
        Geweke_w[i] = geweke_diagnostic(param_vec_sampling, ESS_scale);

        // True/False measure of convergence by comparing Geweke statistic against z-distribution. This is a flipped
        // hypothesis test, i.e. convergence is reached if we fail to reject the null. Hence, we opt for a significance
        // level of alpha=0.2 rather than the traditional alpha=0.05, as higher values are more conservative.
        if ((Geweke_w[i] > -1.282) && (Geweke_w[i] < 1.282)) {
            convergence_pass[i] = true;
        }
    }
}

// Write values to file
void MCMC_diags::write_diagnostics(std::ofstream& out) const
{
    out << "parameter,Geweke,converged,ESS\n";
    for (int i = 0; i < ess_w.size(); ++i) {
        out << "w_" << i << ", " << Geweke_w[i] << ", " << convergence_pass[i] << ", " << ess_w[i] << std::endl;
    }
    
}