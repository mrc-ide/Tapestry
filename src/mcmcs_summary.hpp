#pragma once
#include <string>
#include <vector>
#include "bed.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "particles.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;


class MCMCPointEstimator
{
private:
    // Passed
    const MCMC& mcmc;
    const VCFData& data;
    const NaiveIBDModel& model;
    const IBDContainer& ibd;  // TODO: access from model?
    
    // Computed
    Particle map_particle;
    ViterbiResult viterbi;
    MatrixXi pairwise_ibd;
    vector<BEDRecord> ibd_segments;

public:
    MCMCPointEstimator(
        const MCMC& mcmc,
        const VCFData& data,
        const NaiveIBDModel& model,
        const IBDContainer& ibd
    );

    // TODO: arguably in constructor
    void compute_point_estimate();
    void write_output(const string& output_dir) const;
};

