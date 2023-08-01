#pragma once
#include <vector>
#include "model_fitting.hpp"


// --------------------------------------------------------------------------------
// Model comparison through heuristics
// --------------------------------------------------------------------------------


class ModelCompare
{
private:
    double calc_AIC(double n_params, double log_map);
    double calc_BIC(double n_params, double log_map, double n_data);

    // Model fits (i.e. MAP from MCMC)
    vector<ModelFit> model_fits;
    
    // Data associateed
    int n_data;
    
    // Model associated
    std::vector<int> Ks;
    int n_Ks;
    std::vector<int> n_params;

    // Maximum aposteriori associated
    // All of the below should go into a struct
    std::vector<double> map_logvalues;
    std::vector<double> AIC_values;
    std::vector<double> BIC_values;

    // We will include these to in the summary, for now
    std::vector<double> Keffs;
    std::vector<double> f_ibds;

public:
    // Constructor
    ModelCompare(vector<ModelFit> model_fits, int n_data, vector<int> Ks);

    void write_output(std::string output_csv);
};