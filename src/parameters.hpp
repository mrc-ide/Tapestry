#pragma once

/*
* Encapsulate all fixed model and mcmc parameters
*/
struct Parameters
{
    
    // Model Hyperparameters
    const int K;
    const double e_0;
    const double e_1;
    const double v;
    const double rho;
    const double G;
    const int n_pi_bins;

    // MCMC Parameters
    const double w_proposal_sd;

    Parameters(
        int K, 
        double e_0, 
        double e_1, 
        double v, 
        double rho, 
        double G,
        double w_proposal_sd,
        int n_pi_bins
    );
    
    void print();
};
