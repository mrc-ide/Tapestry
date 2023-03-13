#pragma once


struct Parameters
{
    // Model Hyperparameters
    const int K;                      // COI
    const double e_0;                 // REF -> ALT error probability
    const double e_1;                 // ALT -> REF error probability
    const double v;                   // Beta-binomial dispersion
    const double rho;                 // Recombination rate

    // MCMC Parameters
    const double prop_proposal_sd;    // Proportion-titre update standard deviation

    // Beta-binomial Approx.
    const int n_pi_bins;              // No. discrete error-adjusted WSAF bins

    Parameters(
        int K, 
        double e_0, 
        double e_1, 
        double v, 
        double rho, 
        double prop_proposal_sd,
        int n_pi_bins
    );

    void print() const;
};
