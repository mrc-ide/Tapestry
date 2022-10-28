
#pragma once

#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// class holding all common objects that every particle needs access to
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  std::vector<int> a;
  std::vector<int> r;
  int n_loci;
  
  // other parameters
  std::vector<double> p;
  double c;
  
  // gamma function lookup table
  std::vector<std::vector<double>> betabinom_lookup;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args_data, Rcpp::List args_params);
  void init_betabinom_lookup();
  double get_betabinom(int i, double pi_);
  
};
