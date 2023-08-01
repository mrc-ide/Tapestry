#include <cmath>
#include <iostream>
#include <fstream>
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

