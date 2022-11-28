#include <Rcpp.h>

#include "misc_v15.h"
#include "probability_v17.h"

//------------------------------------------------
// run basic example mcmc
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(int K,
                        double e_0,
                        double e_1,
                        double v,
                        Rcpp::StringVector chroms,
                        Rcpp::NumericVector pos,
                        Rcpp::NumericVector refs,
                        Rcpp::NumericVector alts,
                        Rcpp::NumericVector plafs,
                        Rcpp::NumericVector wsafs,
                        Rcpp::List args_MCMC, Rcpp::List args_progress,
                        Rcpp::List args_functions);
