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
    double target_acceptance;
    int swap_freq;

    Parameters(
        int K, 
        double e_0, 
        double e_1, 
        double v, 
        double rho, 
        double G,
        int n_pi_bins,
        double target_acceptance,
        int swap_freq
    );
    
    void print();
};
