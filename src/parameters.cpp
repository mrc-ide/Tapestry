#include <iostream>
#include "parameters.hpp"
using namespace std;


Parameters::Parameters(
    int K, 
    double e_0, 
    double e_1, 
    double v, 
    double rho, 
    double prop_proposal_sd,
    int n_pi_bins
) : K(K),
    e_0(e_0),
    e_1(e_1),
    v(v),
    rho(rho),
    prop_proposal_sd(prop_proposal_sd),
    n_pi_bins(n_pi_bins)
{}


void Parameters::print() const
{
    cout << "Parameters:" << endl;
    cout << "  K: " << K << endl;
    cout << "  e_0: " << e_0 << endl;
    cout << "  e_1: " << e_1 << endl;
    cout << "  v: " << v << endl;
    cout << "  rho: "  << rho << endl;
    cout << "  Prop SD: " << prop_proposal_sd << endl;
    cout << "  No. pi bins: " << n_pi_bins << endl;
}

