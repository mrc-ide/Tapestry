
#include "System.h"
#include "misc_v15.h"

using namespace std;

//------------------------------------------------
// load parameters and functions
void System::load(Rcpp::List args_data, Rcpp::List args_params) {
  
  // data
  a = rcpp_to_vector_int(args_data["a"]);
  r = rcpp_to_vector_int(args_data["r"]);
  n_loci = a.size();
  
  // other parameters
  p = rcpp_to_vector_double(args_params["p"]);
  c = args_params["c"];
  
}

//------------------------------------------------
// initialise Beta-binomial likelihood lookup table for each locus
void System::init_betabinom_lookup() {
  
  betabinom_lookup = vector<vector<double>>(n_loci, vector<double>(1001));
  for (int i = 0; i < n_loci; ++i) {
    double tmp0 = lgamma(a[i] + r[i] + 1) - lgamma(a[i] + 1) - lgamma(r[i] + 1);
    double tmp1 = lgamma(c) - lgamma(a[i] + r[i] + c);
    for (int j = 0; j < 1001; ++j) {
      double pi_ = j / double(1000);
      double tmp2 = lgamma(a[i] + pi_*c) - lgamma(pi_*c);
      double tmp3 = lgamma(r[i] + (1.0 - pi_)*c) - lgamma((1.0 - pi_)*c);
      betabinom_lookup[i][j] = exp(tmp0 + tmp1 + tmp2 + tmp3);
    }
  }
  
}

//------------------------------------------------
// get Beta-binomial likelihood for locus i at WSAF pi_
double System::get_betabinom(int i, double pi_) {
  
  // convert pi_ to corresponding integer index for lookup table, and return
  // value from table
  int index = round(pi_*1000);
  return betabinom_lookup[i][index];
}
