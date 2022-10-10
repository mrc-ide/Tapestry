
#include <chrono>
#include "main.h"
#include "MCMC.h"

using namespace std;

//------------------------------------------------
// run basic example mcmc
Rcpp::List run_mcmc_cpp(Rcpp::NumericVector x, Rcpp::List args_params,
                        Rcpp::List args_MCMC, Rcpp::List args_progress,
                        Rcpp::List args_functions) {
  
  //-----------------------------
  // SCRATCH PAD
  // This part of the code is where you can play around with different C++ functions
  /*
   FYI, you can write a mult-line
   comment like this
  */
  
  /*
  // this is how you define a basic C++ object:
  int my_integer = 5;
  double my_double = 12.5;
  bool by_logical = true;
  
  // you can define a vector like this:
  vector<int> my_vec(3, 2);   // the "3" means length 3, the "2" fills the vector with the value 2
  
  // we index the vector with square brackets and 0-index notation:
  print(my_vec[0]);
  //print(my_vec[3]);   \\ this woud go off the end of the vector, and so would result in crash
  
  // the print() function is one that I wrote, and that is loaded in from
  // misc. You can list as many print values as you want:
  print("first vector value:", my_vec[0], "\n");
  
  // an matrix is just a vector of vectors. The initialisation syntax here is the same as for a
  // vector, i.e. we initialise a vector of length 2 containing vectors of length 5:
  vector<vector<double>> my_mat(2, vector<double>(5, 9.9));
  
  // this is a rectangular matrix because every vector contains a vector of the
  // same length, but in fact there is no contraint, you can have a "ragged"
  // matrix no problem
  
  // I also wrote functions for printing vectors and matrices
  print_vector(my_vec);
  print_matrix(my_mat);
  
  // there are also a wide range of Rcpp types that we could use. These have
  // slightly different syntax. It can get a bit tricky when you start merging
  // base types and Rcpp types, as you can easily get unexpected behaviour
  Rcpp::NumericVector my_rcpp_vector(4);
  
  // if you want to exit Rcpp back to R at any point you can use this function:
  //Rcpp::stop("error message");
  */
  //-----------------------------
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
  // extract R functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // create MCMC object and load arguments
  MCMC mcmc(x, args_params, args_MCMC, args_progress);
  
  // run burn-in and sampling phases of MCMC
  mcmc.run_mcmc_burnin(update_progress);
  mcmc.run_mcmc_sampling(update_progress);
  
  // end timer
  double t_diff = chrono_timer(t0, "chain completed in ", "\n", !mcmc.silent);
  
  // return outputs in list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("mu_burnin") = mcmc.mu_burnin,
                                      Rcpp::Named("sigma_burnin") = mcmc.sigma_burnin,
                                      Rcpp::Named("mu_sampling") = mcmc.mu_sampling,
                                      Rcpp::Named("sigma_sampling") = mcmc.sigma_sampling,
                                      Rcpp::Named("MC_accept_burnin") = mcmc.MC_accept_burnin,
                                      Rcpp::Named("MC_accept_sampling") = mcmc.MC_accept_sampling,
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}
