#include <algorithm>
#include <iostream>
#include "bed.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "mcmcs_summary.hpp"
#include "models.hpp"
#include "particles.hpp"
using namespace std;


MCMCPointEstimator::MCMCPointEstimator(
    const MCMC& mcmc,
    const VCFData& data,
    const NaiveIBDModel& model,
    const IBDContainer& ibd)
    : mcmc(mcmc), 
    data(data),
    model(model),
    ibd(ibd)
{}


void MCMCPointEstimator::compute_point_estimate()
{
    
    // First, we get the MAP particle
    this->map_particle = mcmc.get_map_particle();
    std::sort(map_particle.ws.begin(), map_particle.ws.end());  // ascending sort
    cout << map_particle.ws << endl;

    // Highest probability IBD path
    // TODO: is this how I want to assign?
    this->viterbi = model.get_viterbi_path(map_particle);

    // Convert to a pairwise IBD matrix
    this->pairwise_ibd = convert_ibd_state_path_to_pairwise(viterbi.path, ibd);

    // IBD segments, a bit laboured
    std::vector<BEDRecord> pair_ibd_segments;
    for (int i = 0; i < ibd.column_index_to_pair.size(); ++i) {
        // Get IBD segments
        pair_ibd_segments = get_ibd_segments(
            data.chrom_names,
            data.pos,
            pairwise_ibd.col(i),
            ibd.pair_names[i]
        );

        // Insert in vector
        this->ibd_segments.insert(
            ibd_segments.end(), 
            pair_ibd_segments.begin(), 
            pair_ibd_segments.end()
        );
    }
}


// TODO: Probably should also write likelihood of point estimate
void MCMCPointEstimator::write_output(const string& output_dir) const
{
    std::string viterbi_csv = output_dir + "/ibd.viterbi_path.csv";
    write_data_with_annotation(
        viterbi_csv,
        data,
        viterbi.path,
        vector<string>{"ibd_viterbi"}
    );

    std::string pairwise_csv = output_dir + "/ibd.viterbi_path.pairwise.csv";
    write_data_with_annotation(
        pairwise_csv,
        data,
        pairwise_ibd,
        ibd.pair_names
    );

    string segment_bed = output_dir + "/ibd.viterbi_path.segments.bed";
    write_bed_records(segment_bed, ibd_segments);
}

