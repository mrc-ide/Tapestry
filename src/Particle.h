
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
  System * s_ptr;
  
  // model parameters
  std::vector<double> mu;
  double sigma;
  
  // likelihood and prior
  double loglike;
  double logprior;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialisers
  void init(System &s);
  
  // likelihood and prior
  double get_loglike(std::vector<double> mu, double sigma);
  double get_logprior(std::vector<double> mu, double sigma);
  
  // update functions
  void update(double beta);
  void update_mu(double beta);
  void update_sigma(double beta);
  
};
