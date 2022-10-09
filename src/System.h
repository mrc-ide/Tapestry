
#pragma once

#include <Rcpp.h>
#include <vector>

//------------------------------------------------
// class holding all common objects that every particle needs access to
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  Rcpp::NumericVector x;
  int n_loci;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::NumericVector x, Rcpp::List args_params);
  
};
