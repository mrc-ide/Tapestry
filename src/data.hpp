#pragma once
#include <string>
#include <vector>
#include <Rcpp.h>
using namespace std;
    

class Data
{
public:
    const int n_loci;
    const vector<string> chroms;
    const vector<int> pos;
    const vector<int> refs;
    const vector<int> alts;
    const vector<double> plafs;
    const vector<double> wsafs;

    Data(
        int n_loci,
        Rcpp::StringVector chroms,
        Rcpp::NumericVector pos,
        Rcpp::NumericVector refs,
        Rcpp::NumericVector alts,
        Rcpp::NumericVector plafs,
        Rcpp::NumericVector wsafs
    );

    void print() const;
};

