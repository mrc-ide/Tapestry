#pragma once
#include <vector>
#include "data.hpp"
#include "libs/eigen-3.4.0/Dense"
#include "typedefs.hpp"
using namespace std;
using Eigen::RowVectorXd;  // TODO: I probably want RowVector


class PanelFreePhaser
{
private:
    // Input data from `tapestry infer`
    int K;
    RowVectorXd& ws; // but this must not be an array
    const InferredIBDPathData& data;
    const double e_0;
    const double e_1;
    const double v;
    
    // Required intermediate arrays for phasing
    const MatrixXi allele_configs;
    IBDContainer ibd;
    ArrayXd wsaf_adj; // the adjusted WSAF values

    // Directly used in computing phase
    std::vector<MatrixXd> sampling_probs;
    MatrixXd betabin_probs;

    // Phasing Outputs
    MatrixXi phased_alleles;
    MatrixXd phased_allele_probs;

public:
    PanelFreePhaser(
        int K, 
        RowVectorXd& ws,
        const InferredIBDPathData& data,
        double e_0,
        double e_1,
        double v);
    void phase();
    void write_output(const std::string& output_dir);
};

