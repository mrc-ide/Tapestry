
#include <chrono>
#include "main.h"
#include "MCMC.h"
#include "System.h"

using namespace std;

//------------------------------------------------
// run basic example mcmc
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
                        Rcpp::List args_MCMC, 
                        Rcpp::List args_progress,
                        Rcpp::List args_functions) 
{   
  // start timer
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
  // extract R functions
  Rcpp::Function update_progress = args_functions["update_progress"];

  // instantiate system
  cout << "Building system..." << endl;
  SystemVCF system;
  system.load(chroms, pos, refs, alts, plafs, wsafs);
  system.set_parameters(K, e_0, e_1, v);
  cout << "Data:" << endl;
  cout << "  No. loci: " << system.n_loci << endl;
  cout << "Parameters:" << endl;
  cout << "  COI: " << system.K << endl;
  cout << "  e_0: " << system.e_0 << endl;
  cout << "  e_1: " << system.e_1 << endl;
  cout << "  v: " << system.v << endl;
  cout << "Done." << endl;

  // create MCMC object and load arguments
  MCMC mcmc(system, args_MCMC, args_progress);
  
  // run burn-in and sampling phases of MCMC
  mcmc.run_mcmc_burnin(update_progress);
  mcmc.run_mcmc_sampling(update_progress);
  
  // end timer
  double t_diff = chrono_timer(t0, "chain completed in ", "\n", !mcmc.silent);
  
  // return outputs in list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("mu_burnin") = mcmc.mu_burnin,
                                      Rcpp::Named("sigma_burnin") = mcmc.sigma_burnin,
                                      Rcpp::Named("w_burnin") = mcmc.w_burnin,
                                      Rcpp::Named("mu_sampling") = mcmc.mu_sampling,
                                      Rcpp::Named("sigma_sampling") = mcmc.sigma_sampling,
                                      Rcpp::Named("w_sampling") = mcmc.w_sampling,
                                      Rcpp::Named("MC_accept_burnin") = mcmc.MC_accept_burnin,
                                      Rcpp::Named("MC_accept_sampling") = mcmc.MC_accept_sampling,
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}
