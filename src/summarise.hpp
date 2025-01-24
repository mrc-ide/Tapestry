#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include "bed.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "parameters.hpp"
#include "libs/eigen-3.4.0/Dense"
#include "libs/json.hpp"
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;


namespace ModelFit {
    

    // ================================================================================
    // Compute sample-level summary statistics based on a model fit
    // 
    // ================================================================================


    struct SampleStatistics
    {
        //std::string sample_name;
        int n_sites;
        int K;
        double Keff;
        std::vector<double> ws;
        double f_ibd;
        double l_ibd;
        double n50_ibd;
        double n_ibd;
        double logposterior;
        double aic;
        double bic;
        double logevidence = 0.0;
        double model_posterior = -1.0;

        SampleStatistics();
    };
    NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(
        SampleStatistics, 
        n_sites, K, Keff, ws, f_ibd, l_ibd, n50_ibd, n_ibd, logposterior, aic, bic, logevidence, model_posterior
    );


    bool compare_by_coi(
        const SampleStatistics& a,
        const SampleStatistics& b
    );

    SampleStatistics load_statistics_from_json(const std::string& json_path);
    void write_statistics_to_json(const SampleStatistics& sample_stats, const std::string& output_dir);


    // ================================================================================
    // Compute pairwise IBD statistics based on a model fit
    // 
    // ================================================================================


    struct PairwiseIBDStatistics
    {
    public:
        // Inputs 
		const std::vector<BEDRecord>& ibd_segments;
		const double genome_length;
		
		// Calculated
		int n_ibd;
		
		std::vector<int> ibd_segment_lengths; // can reserve in itialiser
		double total_ibd = 0;		
		double f_ibd = 0;
		double l_ibd = 0;
		double n50_ibd = 0;

        PairwiseIBDStatistics(
            const std::vector<BEDRecord>& ibd_segments,
            const double genome_length
        );

    private:
        void calc_ibd_segment_lengths_and_total();
		void calc_f_ibd();
		void calc_l_ibd();
		void calc_n50_ibd();
    };


    // ================================================================================
    // Summarise and write various outputs from a model fit
    // 
    // ================================================================================

    
    double calc_aic(double n_params, double log_map);
    double calc_bic(double n_params, double log_map, double n_data);

    class Summariser
    {
    private:
        double calc_effective_coi(VectorXd ws);
        void create_ibd_segments();
        void calc_ibd_summary_stats();
        void write_ibd_profiles(const std::string& output_dir);
        void write_ibd_pairwise_stats(const std::string& output_dir);
        void write_sample_stats(const std::string& output_dir);
    public:
        const Parameters& params;
        const VCFData& data;
        SampleStatistics& sample_stats;                             // Sample summary statistics
        const VectorXi ibd_states;
        IBDContainer ibd;
        MatrixXi ibd_pairwise;                                      // Pairwise IBD profiles
        std::vector<PairwiseIBDStatistics> ibd_pairwise_stats;      // Pairwise IBD summary statistics
        vector<BEDRecord> ibd_segments;                             // Complete list of IBD segments

        Summariser(
            const Parameters& params,
            const VCFData& data,
            SampleStatistics& sample_stats,
            const VectorXi ibd_states,
            const VectorXd ws,
            double logposterior
        );

        void write_output(const std::string& output_dir);
    };
}


// ================================================================================
// Compute model evidence using thermodynamic intergration
// 
// ================================================================================


//TODO: Could probably just be a function
class ModelEvidenceCalculator
{
    private:
        const MCMC& mcmc;
        ModelFit::SampleStatistics& sample_stats;
        double numerically_integrate(ArrayXd xs, ArrayXd ys) const;
        ArrayXd calc_meanloglikelihoods() const;

    public:
        ModelEvidenceCalculator(const MCMC& mcmc, ModelFit::SampleStatistics& sample_stats);
        void calc_logevidence();
};


// ================================================================================
// Merge and compare model fits across COI levels
// 
// ================================================================================


class CompareAcrossCOI
{
    private:
        const std::string& output_dir;
        std::vector<ModelFit::SampleStatistics> sample_stats_vec;
        std::vector<double> posterior;
    public:
        CompareAcrossCOI(const std::string& output_dir);
        void load_sample_statistics();
        void calc_posterior();
        void write_comparison_csv();
};

