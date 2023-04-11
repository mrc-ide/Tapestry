#include <algorithm>
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
    ibd(params.K),  // TODO: better to pass
    ibd_states(ibd_states),
    ibd_pairwise(convert_ibd_state_path_to_pairwise(ibd_states, ibd)),
    n_ibd(0), total_ibd(0.0), f_ibd(0.0), l_ibd(0.0), n50_ibd(0.0)
{
    // Get IBD segments
    // - I could encapsulate
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

    // Compute IBD summary statistics
    n_ibd = ibd_segments.size();
    if (n_ibd == 0) {
        return;
    }

    // Get a vector of lengths
    double length;
    vector<double> ibd_segment_lengths;
    ibd_segment_lengths.reserve(n_ibd);
    for (const BEDRecord& ibd_segment : ibd_segments) {
        length = ibd_segment.end - ibd_segment.start;
        ibd_segment_lengths.emplace_back(length);
        total_ibd += length;
    }
    l_ibd = total_ibd / n_ibd;
    int G = data.genome_length * ibd.column_index_to_pair.size();
    f_ibd = total_ibd / G;

    // N50
    std::sort(ibd_segment_lengths.begin(), ibd_segment_lengths.end());
    double frac_cumsum = 0;
    for (double l : ibd_segment_lengths) {
        frac_cumsum += l / total_ibd;
        if (frac_cumsum > 0.5) {
            n50_ibd = l;
            break;
        }
    }
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
    
    string viterbi_csv = output_dir + "/fit.ibd.path.csv";
    write_data_with_annotation(
        viterbi_csv,
        data,
        ibd_states,
        vector<string>{"ibd_viterbi"}
    );

    string pairwise_csv = output_dir + "/fit.ibd.pairwise.csv";
    write_data_with_annotation(
        pairwise_csv,
        data,
        ibd_pairwise,
        ibd.pair_names
    );

    string segment_bed = output_dir + "/fit.ibd.segments.bed";
    write_bed_records(segment_bed, ibd_segments);

    // TODO: Encapsulate
    string stats_csv = output_dir + "/fit.ibd.stats.csv"; // could be json
    ofstream csv_file(stats_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open file.");
    }
    csv_file << "stat, value\n";
    csv_file << "n_ibd," << std::to_string(n_ibd) << "\n";
    csv_file << "f_ibd," << std::to_string(f_ibd) << "\n";
    csv_file << "l_ibd," << std::to_string(l_ibd) << "\n";
    csv_file << "n50_ibd," << std::to_string(n50_ibd) << "\n";
    csv_file.close();
}




/* Calculate the N50 value from a vector of lengths
*  N50 is a *weighted* median, here weights are
*  the lengths; half of the total length is in
*  segments of < N50.
*
*  Pre: Lengths should be positive; we copy since
*  we intend to sort in place. Should be at least one
*  length.
*  Post: 
*/
// double calc_n50(const VectorXi& lengths)
// {

// }

// def calc_ibd_n50(ibd_state, pos, chroms, from_segs=False):
//     """
//     Calculate the N50 IBD Segment Length
    
//     N50 is a version of the median where each
//     data-point (IBD segment) is weighted (here,
//     by it's length). In this case, we return
//     the shortest IBD segment such that the sum of
//     all smaller IBD segments sums to *greater than half
//     of the total IBD* observed.
    
//     This metric is robust to distributions with many
//     small-valued data points.
    
//     params
//         ibd_state: ndarray, int (n_sites/n_segments)
//             `ibd_state` is a one-dimensional numpy
//             array holding either the IBD segment lengths
//             themselves (if from_segs == True); or from
//             the IBD state at each position along the genome
//             (if from_segs == False).

//     returns
//         N50_length: int
//             The length of the N50 IBD segment.
            
//     """
//     if not from_segs:
//         ibd_segs = get_ibd_segments(ibd_state, pos, chroms)
//     else:
//         ibd_segs = ibd_state
    
//     if len(ibd_segs) > 0:
//         ibd_segs.sort()
//         total_ibd = ibd_segs.sum()
//         cumulative_ibd_frac = (ibd_segs/float(total_ibd)).cumsum()
//         ix = np.argmax(cumulative_ibd_frac > 0.5)
//     else:
//         return np.nan
//     return ibd_segs[ix]