
#include "System.h"
#include "misc_v15.h"

using namespace std;

//------------------------------------------------
// load parameters and functions
void System::load(Rcpp::NumericVector x, Rcpp::List args_params)
{

  // data
  this->x = x;
  n_loci = x.size();

  // allele frequencies
  p = rcpp_to_vector_double(args_params["p"]);
}

/**
 * Load VCF data from as Rcpp data types into native C++ data types
 *
 * At present the system holds:
 * - all relevant *data*
 * - all model hyper-parameters (not updated in MCMC)
 *    - K
 *    - e0, e1
 *    - v
 *    - rho
 * - all pre-computed arrays required for the model, including:
 *    - Haplotype sampling probabilities
 *    - Relevant beta-binomial probabilities
 *    - IBD initiation probabilities
 *    - IBD transition probability array (not scaled by distance)
 *
 */
void SystemVCF::load(
    Rcpp::StringVector chroms,
    Rcpp::NumericVector pos,
    Rcpp::NumericVector refs,
    Rcpp::NumericVector alts,
    Rcpp::NumericVector plafs,
    Rcpp::NumericVector wsafs)
{
  this->chroms = rcpp_to_vector_string(chroms);
  this->pos = rcpp_to_vector_int(pos);
  this->refs = rcpp_to_vector_int(refs);
  this->alts = rcpp_to_vector_int(alts);
  this->plafs = rcpp_to_vector_double(plafs);
  this->wsafs = rcpp_to_vector_double(wsafs);
  n_loci = this->pos.size();
}

void SystemVCF::set_parameters(
    int K,
    double e_0,
    double e_1,
    double v)
{
  this->K = K;
  this->e_0 = e_0;
  this->e_1 = e_1;
  this->v = v;
}