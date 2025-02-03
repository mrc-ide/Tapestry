#pragma once
#include <vector>
#include "particles.hpp"


// --------------------------------------------------------------------------------
// MCMC diagnostics
// --------------------------------------------------------------------------------

class MCMC_diags
{
private:

    // FUNCTIONS
    double mean(const std::vector<double>& x);
    double variance(const std::vector<double>& x);
    double autocorrelation(const std::vector<double>& x, int lag);
    double compute_ESS_geyer(const std::vector<double>& x);
    double geweke_diagnostic(const std::vector<double>& x, double ESS_scale, double fracA = 0.1, double fracB = 0.5);

public:
    // MEMBERS
    RowVectorXd Geweke_w;
    std::vector<bool> convergence_pass;
    RowVectorXd ess_w;
    
    // FUNCTIONS
    // Constructor
    MCMC_diags() {};

    std::vector<double> extract_parameter(const std::vector<Particle>& particle_trace, std::function<double(const Particle&)> extractor);
    void calc_diagnostics(const std::vector<Particle>& particle_trace, const int n_burn_iters);
    void write_diagnostics(std::ofstream& out) const;

    // Destructor
    ~MCMC_diags() {}
};

