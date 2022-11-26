
#pragma once

#include <Rcpp.h>
#include <vector>

/**
 * Outstanding questions:
 * - When instantiating system, can I make sure the sizes of the vectors
 * are correct from the beginning (e.g. length n_loci), to ensure there
 * is no needless copying?
*/

//------------------------------------------------
// class holding all common objects that every particle needs access to
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  Rcpp::NumericVector x;
  int n_loci;
  
  // allele frequencies
  std::vector<double> p;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::NumericVector x, Rcpp::List args_params);
  
};


class SystemVCF
{
public:
  int n_loci;
  std::vector<std::string> chroms; // might want a map
  std::vector<int> pos;
  std::vector<int> refs;
  std::vector<int> alts;
  std::vector<double> plafs;
  std::vector<double> wsafs;

  SystemVCF() {};

  void load(
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs
  );
};