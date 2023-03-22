#include <iostream>
#include "parameters.hpp"
using namespace std;


Parameters::Parameters(
    int K, 
    double e_0, 
    double e_1, 
    double v, 
    double rho, 
    double w_proposal_sd,
    const int n_pi_bins
)
    : K(K),
    e_0(e_0),
    e_1(e_1),
    v(v),
    rho(rho),
    w_proposal_sd(w_proposal_sd),
    n_pi_bins(n_pi_bins)
{};


void Parameters::print()
{
    cout << "Model hyperparameters:" << endl;
    cout << "  K: " << K << endl;
    cout << "  e_0: " << e_0 << endl;
    cout << "  e_1: " << e_1 << endl;
    cout << "  v: " << v << endl;
    cout << "  rho: "  << rho << endl;
    cout << "  n_pi_bins: " << n_pi_bins << endl;
    cout << "MCMC Parameters:" << endl;
    cout << "  w_proposal: " << w_proposal_sd << endl;
};

