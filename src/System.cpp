
#include "System.h"
#include "misc_v15.h"

using namespace std;

//------------------------------------------------
// load parameters and functions
void System::load(Rcpp::NumericVector x, Rcpp::List args_params) {
  
  // data
  this->x = x;
  n_loci = x.size();
  
  // allele frequencies
  p = rcpp_to_vector_double(args_params["p"]);
  
}
