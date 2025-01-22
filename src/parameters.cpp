#include <iostream>
#include "parameters.hpp"
using namespace std;


Parameters::Parameters(
    int K, 
    double e_0, 
    double e_1, 
    double v, 
    double rho,
    double G,
    const int n_pi_bins,
    double target_acceptance,
    int swap_freq
)
    : K(K),
    e_0(e_0),
    e_1(e_1),
    v(v),
    rho(rho),
    G(G),
    n_pi_bins(n_pi_bins),
    target_acceptance(target_acceptance),
    swap_freq(swap_freq)
{};


void Parameters::print()
{
    cout << "Model hyperparameters:" << endl;
    cout << "  K: " << K << endl;
    cout << "  e_0: " << e_0 << endl;
    cout << "  e_1: " << e_1 << endl;
    cout << "  v: " << v << endl;
    cout << "  rho: "  << rho << endl;
    cout << "  G: " << G << endl;
    cout << "  n_pi_bins: " << n_pi_bins << endl;
    cout << "MCMC Parameters:" << endl;
    cout << "  target_acceptance: " << target_acceptance << endl;
    cout << "  swap_freq: " << swap_freq << endl;
};

