
#pragma once

#include "System.h"
#include "misc_v15.h"
#include "Particle.h"

#include <Rcpp.h>


//------------------------------------------------
// class defining MCMC
class MCMC {
  
public:
  // PUBLIC OBJECTS
  
  // system object
  System s;
  
  // MCMC parameters
  int burnin;
  int samples;
  
  // acceptance rates
  std::vector<double> MC_accept_burnin;
  std::vector<double> MC_accept_sampling;
  
  // rung parameters
  int rungs;
  std::vector<int> rung_order;
  std::vector<double> beta;
  
  // vector of particles
  std::vector<Particle> particle_vec;
  
  // misc parameters
  bool pb_markdown;
  bool silent;
  
  // objects for storing results
  std::vector<std::vector<double>> mu_burnin;
  std::vector<double> sigma_burnin;
  std::vector<double> w_burnin;
  std::vector<std::vector<double>> mu_sampling;
  std::vector<double> sigma_sampling;
  std::vector<double> w_sampling;
  
  // progress bars
  Rcpp::List args_progress;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  MCMC(Rcpp::List args_data, Rcpp::List args_params, Rcpp::List args_MCMC,
       Rcpp::List args_progress);
  
  // member functions
  void run_mcmc_burnin(Rcpp::Function update_progress);
  void run_mcmc_sampling(Rcpp::Function update_progress);
  void coupling(std::vector<double> &MC_accept);
  
};
