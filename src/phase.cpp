#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include "betabin.hpp"
#include "combinatorics.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "phase.hpp"
#include "sampling.hpp"
using Eigen::RowVectorXd;  // TODO: I probably want RowVector


PanelFreePhaser::PanelFreePhaser(
    int K,
    RowVectorXd& ws,
    const InferredIBDPathData& data,
    const double e_0,
    const double e_1,
    const double v)
    : K(K),
    ws(ws), // needs to come in as a vector
    data(data),
    e_0(e_0), e_1(e_1), v(v),
    allele_configs(create_powerset(K)),
    ibd(K)
{
    // Compute adjusted WSAF
    // Note: I believe this comes out as a COLUMN array?
    ArrayXd wsaf = (allele_configs.cast<double>() * ws.transpose()).array();
    wsaf_adj = (1 - wsaf) * e_0 + (1 - e_1) * wsaf;

    // Compute sampling probabilities
    sampling_probs = create_sampling_probs(data.plafs, allele_configs, ibd.states);
    betabin_probs = Betabinomial::calc_prob_matrix(
        data.refs,
        data.alts,
        wsaf_adj,
        v
    );

    // Prepare space
    phased_allele_probs = MatrixXd::Constant(data.n_sites, K, -1.0);
    phased_alleles = MatrixXi::Constant(data.n_sites, K, -1);
}


void PanelFreePhaser::phase()
{
    double total;
    for (int i = 0; i < data.n_sites; ++i) {

        // Get the MAP IBD state
        int j = data.ibd_path[i];
        
        // Compute the probability of each allele configuration
        ArrayXd allele_probs = betabin_probs.row(i).array() * sampling_probs[i].transpose().row(j).array();
        total = allele_probs.sum();
        allele_probs /= total;

        // Summarise...
        phased_allele_probs.row(i) = (allele_configs.cast<double>().array().colwise() * allele_probs).colwise().sum();
    }
}


void PanelFreePhaser::write_output(const std::string& output_dir)
{
    std::string phase_csv = output_dir + "/fit.ibd.phase.csv";
    std::vector<std::string> strain_names;
    for (int k = 1; k <= K; ++k) {
        ostringstream strain_name;
        strain_name << "strain" << k;
        strain_names.push_back(strain_name.str());
    }
    write_data_with_annotation(
        phase_csv,
        data,
        phased_allele_probs,
        strain_names
    );
}

