#include <iostream>
#include "parameters.hpp"
using namespace std;


Parameters::Parameters(int K)
    : K(K),
    e_0(0.01),
    e_1(0.05),
    v(500),
    rho(10)
{
    print();
}


Parameters::Parameters(int K, double e_0, double e_1, double v, double rho)
    : K(K),
    e_0(e_0),
    e_1(e_1),
    v(v),
    rho(rho)
{
    print();
}


void Parameters::print()
{
    cout << "Parameters:" << endl;
    cout << "  K: " << K << endl;
    cout << "  e_0: " << e_0 << endl;
    cout << "  e_1: " << e_1 << endl;
    cout << "  v: " << v << endl;
    cout << "  rho: "  << rho << endl;
}

