#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "model_fitting.hpp"
#include "model_compare.hpp"
using namespace std;


// --------------------------------------------------------------------------------
// Model comparison through heuristics
// --------------------------------------------------------------------------------


ModelCompare::ModelCompare(vector<ModelFit> model_fits, int n_data, vector<int> Ks)
    : model_fits(model_fits),
    n_data(n_data),
    Ks(Ks),
    n_Ks(Ks.size()),
    n_params(Ks), // TODO: For now, maybe we need to derive something else here
    // TODO: below are all correct size, but have 0s as values
    map_logvalues(n_Ks),
    AIC_values(n_Ks),
    BIC_values(n_Ks),
    Keffs(n_Ks),
    f_ibds(n_Ks)
{
    for (int j = 0; j < n_Ks; ++j) {
        map_logvalues[j] = model_fits[j].logposterior;
        AIC_values[j] = calc_AIC(n_params[j], map_logvalues[j]);
        BIC_values[j] = calc_BIC(n_params[j], map_logvalues[j], n_data);
        Keffs[j] = model_fits[j].Keff;
        f_ibds[j] = model_fits[j].f_ibd;
    }
}

double ModelCompare::calc_AIC(double n_params, double log_map)
{
  return 2 * n_params - 2 * log_map;
}


double ModelCompare::calc_BIC(double n_params, double log_map, double n_data)
{
  return n_params * std::log(n_data) - 2 * log_map;
}


void ModelCompare::write_output(std::string output_csv)
{
  
  // Create file
  std::ofstream csv_file(output_csv);
  if (!csv_file.is_open()) {
    throw std::invalid_argument("Can't open file.");
  }

  csv_file << "K,log_map,AIC,BIC,Keff,f_ibd\n";
  for (int j = 0; j < n_Ks; ++j) {  // order is assumed
    csv_file << j + 1 << ",";
    csv_file << map_logvalues[j] << ",";
    csv_file << AIC_values[j] << ",";
    csv_file << BIC_values[j] << ",";
    csv_file << Keffs[j] << ",";
    csv_file << f_ibds[j] << "\n";
  }

  csv_file.close();
}


// --------------------------------------------------------------------------------
// Compute model evidence using Thermodynamic Evidence
// --------------------------------------------------------------------------------


ModelEvidence::ModelEvidence(const std::vector<std::unique_ptr<ParallelTempering>>& mcmc_ptrs)
  : n_mcmcs(mcmc_ptrs.size()),
  mcmc_ptrs(mcmc_ptrs),
  logevidences(ArrayXd::Constant(n_mcmcs, 0.0)),
  posterior(ArrayXd::Constant(n_mcmcs, -9999.0))
{}


double ModelEvidence::integrate_numerically(ArrayXd xs, ArrayXd ys) const
{
  if (xs.size() != ys.size()) {
    throw std::invalid_argument("For numeric integration x- and y-values must be same length.");
  }

  // TODO: could also implement more explicitly using Eigen features
  // TODO: note we are also assuming they are sorted, either ascending or descending
  double area = 0;
  for (int i = 1; i < xs.size(); ++i) {
    double height = 0.5 * (ys(i) + ys(i-1));
    double width = xs(i) - xs(i-1);
    area += (height * width);
  }

  return area;
}


ArrayXd ModelEvidence::calc_meanloglikelihoods(const ParallelTempering& mcmc) const
{
  // TODO: would be cleaner using Eigen views
  VectorXd m = mcmc.loglikelihoods(Eigen::seq(mcmc.n_burn_iters, Eigen::last), Eigen::all).colwise().mean();
  return m.array();
}


void ModelEvidence::calc_logevidences()
{
  for (int j = 0; j < n_mcmcs; ++j) {
    ArrayXd meanloglikelihoods = calc_meanloglikelihoods(*mcmc_ptrs[j]);

    // TODO: This is quite clumsy, but it is also hard to do otherwise
    ArrayXd betas = ArrayXd::Constant(mcmc_ptrs[j]->n_temps, -1.0);
    for (int k = 0; k < betas.size(); ++k) {
      betas[k] = mcmc_ptrs[j]->temps[k].beta;
    }

    // Integrate the mean-loglikelihood from beta [0, 1] for the MCMC
    logevidences[j] = integrate_numerically(betas, meanloglikelihoods);
  }
}


void ModelEvidence::calc_posterior()
{
  // TODO: verify this is best way to handle underflow
  ArrayXd temp_vals = logevidences - logevidences.maxCoeff(); // preventing underflow
  posterior = temp_vals.exp();
  posterior /= posterior.sum();
}


void ModelEvidence::calc_summary()
{
  calc_logevidences();
  calc_posterior();
}


void ModelEvidence::write_output(std::string output_csv)
{
  // Create file
  std::ofstream csv_file(output_csv);
  if (!csv_file.is_open()) {
    throw std::invalid_argument("Error, can't open file.");
  }

  // TODO: Order assumed. Improve this output.
  csv_file << "K,logevidence,posterior\n";
  for (int j = 0; j < n_mcmcs; ++j) {
    csv_file << j + 1 << ",";
    csv_file << logevidences[j] << ",";
    csv_file << posterior[j] << "\n";
  }

  csv_file.close();
}


