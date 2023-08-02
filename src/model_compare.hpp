#pragma once
#include <memory>
#include <vector>
#include "mcmcs.hpp"
#include "model_fitting.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::ArrayXd;


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


// --------------------------------------------------------------------------------
// Compute model evidence using Thermodynamic Evidence
// --------------------------------------------------------------------------------


/* Compute model evidence using Thermodynamic Integration
*
* Note that another implementation would work with just a single MCMC chain,
* and then the comparison part could occur with outputs from processing
* a single chain
*
*/
class ModelEvidence
{
private:
    // Pointer to vector of ParallelTempering MCMCs
    int n_mcmcs;
    const std::vector<std::unique_ptr<ParallelTempering>>& mcmc_ptrs;

    // Computed
    ArrayXd logevidences;       // The log-evidence P(model|data) for each MCMC
    ArrayXd posterior;          // 

    /* Numerical integration using the Trapezoidal Method
    *
    *  Pre: An equal number of x-values `xs` and y-values `ys` must
    *  be passed. These (x, y) pairs should be sorted by x-value.
    *  Post: Approximate area under the curve of y.
    *
    *  TODO: not doing extrapolation yet.
    */
    double integrate_numerically(ArrayXd xs, ArrayXd ys) const;
    ArrayXd calc_meanloglikelihoods(const ParallelTempering& mcmc) const;
    void calc_logevidences();

    /* Compute the posterior distribution over the MCMC models
    * 
    * This is deriveed from the log-evidence values that are
    * calculated using thermodynamic integration.
    *
    */
    void calc_posterior();

public:
    ModelEvidence(const std::vector<std::unique_ptr<ParallelTempering>>& mcmc_ptrs);

    /* Interface method; compute the log-evidences annd
    *  the posterior
    */
    void calc_summary();

    // Write
    void write_output(std::string output_csv);
};

