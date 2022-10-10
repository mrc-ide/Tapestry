
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s) {
  
  // pointer to system object
  s_ptr = &s;
  
  // initialise model parameters
  mu = {0.1, 0.9};
  sigma = 0.1;
  
  // initialise likelihood and prior
  loglike = get_loglike(mu, sigma);
  logprior = get_logprior(mu, sigma);
  
}

//------------------------------------------------
// calculate loglikelihood given parameter inputs
double Particle::get_loglike(vector<double> mu, double sigma) {
  
  // calculate log-like over all data
  double ret = 0.0;
  for (int i = 0; i < s_ptr->n_loci; ++i) {
    double tmp = 0.5*R::dnorm4(s_ptr->x[i], mu[0], sigma, false) + 0.5*R::dnorm4(s_ptr->x[i], mu[1], sigma, false);
    ret += log(tmp);
  }
  
  return ret;
}

//------------------------------------------------
// calculate logprior given parameter inputs
double Particle::get_logprior(vector<double> mu, double sigma) {
  
  // calculate log-prior
  double ret = R::dnorm4(mu[0], 0.0, 1.0, true) +
    R::dnorm4(mu[1], 0.0, 1.0, true) +
    R::dexp(sigma, 1.0, true);
  
  return ret;
}

//------------------------------------------------
// propose new parameter values and accept/reject
void Particle::update(double beta) {
  
  // distinct update steps for each free parameter
  update_mu(beta);
  update_sigma(beta);
  
}

//------------------------------------------------
// update component means
void Particle::update_mu(double beta) {
  
  // propose a new value by drawing from normal around current value
  vector<double> mu_prop = mu;
  for (int i = 0; i < 2; ++i) {
    
    // propose new value of mu
    mu_prop[i] = rnorm1(mu[i], 0.1);
    
    // calculate likelihood and prior of proposed value
    double loglike_prop = get_loglike(mu_prop, sigma);
    double logprior_prop = get_logprior(mu_prop, sigma);
    
    // calculate Metropolis-Hastings ratio
    double MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    
    // accept or reject move
    if (log(R::runif(0, 1)) < MH) {
      mu[i] = mu_prop[i];
      loglike = loglike_prop;
      logprior = logprior_prop;
    } else {
      mu_prop[i] = mu[i];
    }
    
  }
  
}

//------------------------------------------------
// update component SDs
void Particle::update_sigma(double beta) {
  
  // propose a new value by drawing from reflected normal around current value
  double sigma_prop = rnorm1(sigma, 0.1);
  if (sigma_prop < 0) {
    sigma_prop = -sigma_prop;
  }
  
  // calculate likelihood and prior of proposed value
  double loglike_prop = get_loglike(mu, sigma_prop);
  double logprior_prop = get_logprior(mu, sigma_prop);
  
  // calculate Metropolis-Hastings ratio
  double MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  
  // accept or reject move
  if (log(R::runif(0, 1)) < MH) {
    sigma = sigma_prop;
    loglike = loglike_prop;
    logprior = logprior_prop;
  }
  
}
