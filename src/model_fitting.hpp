#include <vector>
#include "bed.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;




/*
* The idea here is pretty simple;
* you pass your fit *proportions*
* and your fit *IBD* states
* and this class can automatically compute
* the summary statistics for you
*
* Probably can also take a write method
* Compared to before, what is the difference?
* Could fit models in different ways and re-use this class
*
*/
class ModelFit
{
private:
    // Proportions
    double static calc_effective_coi(VectorXd ws);

    // IBD-related
    void calc_ibd_summary_stats();


public:
    // MEMBERS
    // General
    const Parameters& params;
    const VCFData& data;
    const double logposterior;

    // Proportions
    const VectorXd ws;
    double Keff;

    // IBD
    const bool has_ibd;
    const VectorXi ibd_states;
    IBDContainer ibd;
    MatrixXi ibd_pairwise;
    vector<BEDRecord> ibd_segments;
    int n_ibd;      // Number of IBD segments
    double total_ibd, f_ibd, l_ibd, n50_ibd;  // could also be for each pair, how F? Need G length

    
    // METHODS
    // Constructor without IBD
    ModelFit(
        const Parameters& params,
        const VCFData& data,
        const double logposterior,
        const VectorXd ws  // TODO: need to be copied if created in fuction
    );
    // Constructor with IBD
    ModelFit(
        const Parameters& params,
        const VCFData& data,
        const double logposterior,
        const VectorXd ws,  // TODO: need to be copied if created in fuction
        const VectorXi ibd_states
    );

    void write_output(const string& output_dir);
};


// ModelFit compute_MAP_model_fit(    
//     const Parameters& params,
//     const NaiveIBDModel& model, 
//     const MCMC& mcmc
// );