
#pragma once
#include <Rcpp.h>
#include <vector>
#include "typedefs.hpp"
using namespace std;


class SystemVCF
{
private:
  void create_strains();
  void create_hap_configs();
  void create_ibd_configs();
  void create_hap_sampling_probs();

public:
  // HYPER PARAMETERS
  int K;                                  // COI
  double e_0, e_1;                        // Sequencing error rates
  double v;                               // WSAF dispersion

  // DATA
  int n_loci;
  vector<string> chroms;             // Chromosome names; TODO as map
  vector<int> pos;                   // SNP positions
  vector<int> refs;                  // Reference allele counts
  vector<int> alts;                  // Alternative allele counts
  vector<double> plafs;              // Population-level allele frequencies
  vector<double> wsafs;              // Within-sample allele frequencies

  // MODEL FEATURES
  vector<int> strains;
  matrix_2d_int hap_configs;
  vector<vector<vector<int>>> ibd_configs;
  matrix_3d_double hap_sampling_probs;

  // CONSTRUCTORS
  SystemVCF(
    int K, 
    int n_loci,
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
  );
  SystemVCF(
    int K, 
    int n_loci,
    double e_0,
    double e_1,
    double v,
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
  );

  // METHODS
  void print();
  void precompute_arrays();
};

