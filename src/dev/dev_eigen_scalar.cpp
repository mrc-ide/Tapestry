#include <iostream>
#include "betabin.hpp"
#include "parameters.hpp"
#include "particles.hpp"
#include "libs/eigen-3.4.0/Dense"
using namespace std;
using Eigen::ArrayXi;
using Eigen::ArrayXd;
using Eigen::VectorXd;

int main(int argc, char* argv[])
{

    // Instantiate parameters
    Parameters params(3, 0.01, 0.05, 100.0, 10.0, 0.5, 100);
    params.print();

    // Instantiate VCFData
    VCFData data(argv[1], argv[2]);
    data.print();

    // Create Betabinomial lookup
    cout << "Creating betabinomial array LOG..." << endl;
    BetabinomialArray betabin_lookup(params, data);
    cout << "Example values: " << endl;
    double q = 0.5;
    for (int i = 0; i < 10; ++i) {
        cout << betabin_lookup(i, q) << "\t";
    }

    // Create Betabinomial lookup
    cout << "Creating betabinomial array LIN..." << endl;
    BetabinomialArray betabin_lookup_lin(params, data, false);
    cout << "Example values: " << endl;
    for (int i = 0; i < 10; ++i) {
        cout << betabin_lookup_lin(i, q) << "\t";
    }

    //  Next... we want to slice...
    


}