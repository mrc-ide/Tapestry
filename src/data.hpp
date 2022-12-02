#pragma once
#include <string>
#include <vector>
#include <Rcpp.h>
using namespace std;
    

class Data
{
public:
    int n_loci;
    vector<string> chroms;
    vector<int> pos;
    vector<int> refs;
    vector<int> alts;
    vector<double> plafs;
    vector<double> wsafs;

    Data(
        int n_loci,
        Rcpp::StringVector chroms,
        Rcpp::NumericVector pos,
        Rcpp::NumericVector refs,
        Rcpp::NumericVector alts,
        Rcpp::NumericVector plafs,
        Rcpp::NumericVector wsafs
    );

    void print();
};

