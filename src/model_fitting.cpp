#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "model_fitting.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;


// TODO: Lots of unnecessary copies in the arrangement below
// - First, proportions and IBD states are copied INTO ModelFit
// - Then, ModelFit is copied out of the `compute_MAP_model_fit()`


// ================================================================================
// Model fitting functions
//
// ================================================================================


/* 
* Can we return this as reference?
* Does that make sense?
* Idea would be to initialise outside (with params),
* and then compute the fit inside; by passing as reference
*/
// ModelFit compute_MAP_model_fit(
//     const Parameters& params,
//     const NaiveIBDModel& model,
//     const MCMC& mcmc
// )
// {

//     // // Get MAP proportions
//     // Particle map_particle = mcmc.get_map_particle();
//     // std::sort(map_particle.ws.begin(), map_particle.ws.end());

//     // // Get MAP IBD conditioned on proportions
//     // ViterbiResult viterbi = model.get_viterbi_path(map_particle);

//     // return ModelFit(
//     //     params,
//     //     viterbi.logposterior,
//     //     map_particle.ws,
//     //     viterbi.path
//     // );
// }


// ================================================================================
// Model fitting encapsulation
//
// ================================================================================

// --------------------------------------------------------------------------------
// Private
// --------------------------------------------------------------------------------


// TODO: Better to put this somewhere else, quite general
double ModelFit::calc_effective_coi(VectorXd ws)
{
    return 1.0 / ws.array().pow(2.0).sum();
}


void ModelFit::calc_ibd_summary_stats()
{}


// --------------------------------------------------------------------------------
// Public
// --------------------------------------------------------------------------------


ModelFit::ModelFit(
    const Parameters& params,
    const VCFData& data,
    const double logposterior,
    const VectorXd ws)
    : params(params),
    data(data),
    logposterior(logposterior),
    ws(ws),
    Keff(calc_effective_coi(ws)),
    has_ibd(false)
{}


ModelFit::ModelFit(     
    const Parameters& params,
    const VCFData& data,
    const double logposterior,
    const VectorXd ws,
    const VectorXi ibd_states)
    : params(params),
    data(data),
    logposterior(logposterior),
    ws(ws),
    Keff(calc_effective_coi(ws)),
    has_ibd(true),
    ibd(IBDContainer(params.K)),  // TODO: better to pass
    ibd_states(ibd_states),
    ibd_pairwise(convert_ibd_state_path_to_pairwise(ibd_states, ibd))
{
    // Get IBD segments
    vector<BEDRecord> pair_ibd_segments;
    for (int i = 0; i < ibd.column_index_to_pair.size(); ++i) {
        // Get IBD segments
        pair_ibd_segments = get_ibd_segments(
            data.chrom_names,
            data.pos,
            ibd_pairwise.col(i),
            ibd.pair_names[i]
        );

        // Insert in vector
        ibd_segments.insert(
            ibd_segments.end(), 
            pair_ibd_segments.begin(), 
            pair_ibd_segments.end()
        );
    }

    // IBD summary statistics
    n_ibd = 0;
    int total_l_ibd = 0;
    for (const BEDRecord& ibd_segment : ibd_segments) {
        total_l_ibd += (ibd_segment.end - ibd_segment.start);
        n_ibd += 1;
    }
    l_ibd = total_l_ibd / n_ibd;
}

// Behaviour will vary depending on whether or not it has IBD
void ModelFit::write_output(const string& output_dir)
{
    // Write proportions
    string prop_csv = output_dir + "/fit.proportions.csv";
    write_proportions(prop_csv, ws);

    // Write IBD, if provided
    if (!has_ibd) {
        return;
    }

    string pairwise_csv = output_dir + "/fit.ibd.pairwise.csv";
    write_data_with_annotation(
        pairwise_csv,
        data,
        ibd_pairwise,
        ibd.pair_names
    );

    string segment_bed = output_dir + "/fit.ibd.segments.bed";
    write_bed_records(segment_bed, ibd_segments);
}

