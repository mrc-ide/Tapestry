#include "MCMC.h"
//#include "Particle.h"

using namespace std;


/**
 * Oustanding questions:
 * - Passing the system object instead of instantiating will reduce coupling
 *  - If you change System, you do *not* need to change MCMC
 * - However, we want to make sure we are passing by reference; i.e. no
 * memory copying
 * - I am not 100% sure I am doing that correctly
*/


//------------------------------------------------
// constructor
MCMC::MCMC(SystemVCF &system, Rcpp::List args_MCMC, Rcpp::List args_progress)
  : system(system)
{
  
  // // load data and parameters into system object
  // // s.load(x, args_params);
  // this->system = system;
  
  // extract MCMC parameters
  burnin = args_MCMC["burnin"];
  samples = args_MCMC["samples"];
  
  // extract rung parameters
  beta = rcpp_to_vector_double(args_MCMC["beta"]);
  rungs = beta.size();
  
  // extract misc parameters
  pb_markdown = args_MCMC["pb_markdown"];
  silent = args_MCMC["silent"];
  
  // extract progress bars
  this->args_progress = args_progress;
  
  // initialise acceptance rates
  MC_accept_burnin = vector<double>(rungs - 1);
  MC_accept_sampling = vector<double>(rungs - 1);
  
}

//------------------------------------------------
// run burn-in phase of PT MCMC
void MCMC::run_mcmc_burnin(Rcpp::Function update_progress) {
  
  // initialise rung order
  rung_order = seq_int(0, rungs - 1);
  
  // initialise particles
  particle_vec = vector<Particle>(rungs);
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].init(system);
  }
  
  // initialise objects for storing results
  mu_burnin = vector<vector<double>>(burnin, vector<double>(2));
  sigma_burnin = vector<double>(burnin);
  w_burnin = vector<double>(burnin);
  
  // load initial values into store objects
  int r_cold = rungs - 1;
  mu_burnin[0] = particle_vec[r_cold].mu;
  sigma_burnin[0] = particle_vec[r_cold].sigma;
  w_burnin[0] = particle_vec[r_cold].w;
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!silent) {
    print("burn-in");
  }
  
  // loop through burn-in iterations
  for (int rep = 1; rep < burnin; ++rep) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through rungs and update particles
    for (int i = 0; i < rungs; ++i) {
      int r = rung_order[i];
      particle_vec[r].update(beta[i]);
    }
    
    // perform Metropolis coupling
    coupling(MC_accept_burnin);
    
    // store results
    int r_cold = rung_order[rungs - 1];
    mu_burnin[rep] = particle_vec[r_cold].mu;
    sigma_burnin[rep] = particle_vec[r_cold].sigma;
    w_burnin[rep] = particle_vec[r_cold].w;
    
    // update progress bars
    if (!silent) {
      update_progress_cpp(args_progress, update_progress, "pb_burnin", rep + 1, burnin, !pb_markdown);
    }
    
  }  // end MCMC iterations
  
}

//------------------------------------------------
// run sampling phase of PT MCMC
void MCMC::run_mcmc_sampling(Rcpp::Function update_progress) {
  
  // initialise objects for storing results
  mu_sampling = vector<vector<double>>(samples, vector<double>(2));
  sigma_sampling = vector<double>(samples);
  w_sampling = vector<double>(samples);
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!silent) {
    print("sampling");
  }
  
  // loop through sampling iterations
  for (int rep = 0; rep < samples; ++rep) {
    
    // loop through rungs and update particles
    for (int i = 0; i < rungs; ++i) {
      int r = rung_order[i];
      particle_vec[r].update(beta[i]);
    }
    
    // perform Metropolis coupling
    coupling(MC_accept_sampling);
    
    // store results
    int r_cold = rung_order[rungs - 1];
    mu_sampling[rep] = particle_vec[r_cold].mu;
    sigma_sampling[rep] = particle_vec[r_cold].sigma;
    w_sampling[rep] = particle_vec[r_cold].w;
    
    // update progress bars
    if (!silent) {
      update_progress_cpp(args_progress, update_progress, "pb_samples", rep + 1, samples, !pb_markdown);
    }
    
  }  // end MCMC iterations
  
}

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void MCMC::coupling(std::vector<double> &MC_accept) {
  
  // loop over rungs, each time proposing a swap with the next rung along
  for (int i = 0; i < (rungs - 1); ++i) {
    
    // define rungs of interest
    int r1 = rung_order[i];
    int r2 = rung_order[i + 1];
    
    // get log-likelihoods and beta values of the two rungs
    double loglike1 = particle_vec[r1].loglike;
    double loglike2 = particle_vec[r2].loglike;
    
    double beta1 = beta[i];
    double beta2 = beta[i + 1];
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2 - loglike1)*(beta1 - beta2);
    
    // get acceptance rate in linear space (i.e. out of log space). Truncate at
    // 1.0 for values greater than 1.0, and draw whether to accept the move.
    // Only do the draw if the acceptance ratio is < 1
    double acceptance_linear = 1.0;
    bool accept_move = true;
    if (acceptance < 0) {
      acceptance_linear = exp(acceptance);
      accept_move = (R::runif(0, 1) < acceptance_linear);
    }
    
    // update acceptance rates
    MC_accept[i] += acceptance_linear;
    
    // implement swap
    if (accept_move) {
      rung_order[i] = r2;
      rung_order[i + 1] = r1;
    }
    
  } // end loop over rungs
  
}
