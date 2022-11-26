
#pragma once

#include "System.h"
#include "misc_v15.h"
#include "probability_v17.h"

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  SystemVCF * s_ptr;
  
  // model parameters
  std::vector<double> mu;
  double sigma;
  double w;
  
  // likelihood and prior
  double loglike;
  double logprior;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialisers
  void init(SystemVCF &system);
  
  // likelihood and prior
  double get_loglike(std::vector<double> mu, double sigma, double w);
  double get_logprior(std::vector<double> mu, double sigma, double w);
  
  // update functions
  void update(double beta);
  void update_mu(double beta);
  void update_sigma(double beta);
  void update_w(double beta);
  
};
