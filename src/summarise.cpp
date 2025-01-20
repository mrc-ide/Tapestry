#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <regex>
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "summarise.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::VectorXi;
using Eigen::MatrixXi;
using namespace std;
namespace fs = std::filesystem;


namespace ModelFit {
    SampleStatistics::SampleStatistics() {};

    bool compare_by_coi(
        const SampleStatistics& a,
        const SampleStatistics& b)
    {
        return a.K < b.K;
    }

    SampleStatistics load_statistics_from_json(const std::string& json_path) {
        if (!fs::exists(json_path)) {
            throw std::invalid_argument("File '" + json_path + "' does not exist.");
        }
        nlohmann::json stats_json;
        std::ifstream file(json_path);
        file >> stats_json;
        auto sample_stats = stats_json.template get<SampleStatistics>();
        return sample_stats;
    }

    void write_statistics_to_json(const SampleStatistics& sample_stats, const std::string& output_dir)
    {
        nlohmann::json j = sample_stats;
        std::string sample_json = output_dir + "/fit.sample_stats.json";
        std::ofstream output_json(sample_json);
        output_json << std::setw(4) << j << std::endl;
    }

    double calc_aic(double n_params, double log_map)
    {
          return 2 * n_params - 2 * log_map;
    }
    double calc_bic(double n_params, double log_map, double n_data) 
    {
      return n_params * std::log(n_data) - 2 * log_map;
    }

    double Summariser::calc_effective_coi(VectorXd ws)
    {
        return 1.0 / ws.array().pow(2.0).sum();
    }

    void Summariser::create_ibd_segments() {
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
    }

    void Summariser::calc_ibd_summary_stats() {
        // Check if there is any IBD; if not default and exit.
        sample_stats.n_ibd = ibd_segments.size();
        if (sample_stats.n_ibd == 0) {
            sample_stats.f_ibd = 0.0;
            sample_stats.l_ibd = 0.0;
            sample_stats.n50_ibd = 0.0;
            return;
        }

        // Get a vector of IBD segment lengths
        double length;
        double total_ibd = 0;
        vector<double> ibd_segment_lengths;
        ibd_segment_lengths.reserve(sample_stats.n_ibd);
        for (const BEDRecord& ibd_segment : ibd_segments) {
            length = ibd_segment.end - ibd_segment.start;
            ibd_segment_lengths.emplace_back(length);
            total_ibd += length;
        }

        // Compute length and fraction IBD
        sample_stats.l_ibd = total_ibd / sample_stats.n_ibd;
        int G = data.genome_length * ibd.column_index_to_pair.size();
        sample_stats.f_ibd = total_ibd / G;

        // Compute N50 segment length
        std::sort(ibd_segment_lengths.begin(), ibd_segment_lengths.end());
        double frac_cumsum = 0;
        for (double l : ibd_segment_lengths) {
            frac_cumsum += l / total_ibd;
            if (frac_cumsum > 0.5) {
                sample_stats.n50_ibd = l;
                break;
            }
        }
    }

    Summariser::Summariser(
        const Parameters& params,
        const VCFData& data,
        SampleStatistics& sample_stats,
        const VectorXi ibd_states,
        const VectorXd ws,
        double logposterior) 
        : params(params),
        data(data),
        sample_stats(sample_stats),
        ibd_states(ibd_states),
        ibd(params.K),
        ibd_pairwise(convert_ibd_state_path_to_pairwise(ibd_states, ibd))
    {
        
        sample_stats.n_sites = data.n_sites;
        sample_stats.K = params.K;
        sample_stats.Keff = calc_effective_coi(ws);
        std::vector<double> v(ws.data(), ws.data() + ws.size());
        sample_stats.ws = v;

        create_ibd_segments();
        calc_ibd_summary_stats();

        sample_stats.logposterior = logposterior;
        sample_stats.aic = calc_aic(params.K - 1, logposterior);
        sample_stats.bic = calc_bic(params.K - 1, logposterior, data.n_sites);
    }

    void Summariser::write_ibd_profiles(const std::string& output_dir)
    {
        // std::string viterbi_csv = output_dir + "/fit.ibd.path.csv";
        // write_data_with_annotation(
        //     viterbi_csv,
        //     data,
        //     ibd_states,
        //     vector<string>{"ibd_viterbi"}
        // );

        std::string pairwise_csv = output_dir + "/fit.ibd.pairwise.csv";
        write_data_with_annotation(
            pairwise_csv,
            data,
            ibd_pairwise,
            ibd.pair_names
        );

        string segment_bed = output_dir + "/fit.ibd.segments.bed";
        write_bed_records(segment_bed, ibd_segments);
    }

    void Summariser::write_sample_stats(const std::string& output_dir)
    {
        nlohmann::json j = sample_stats;
        std::string sample_json = output_dir + "/fit.sample_stats.json";
        std::ofstream output_json(sample_json);
        output_json << std::setw(4) << j << std::endl;
    }

    void Summariser::write_output(const std::string& output_dir)
    {
        write_ibd_profiles(output_dir);
        //write_sample_stats(output_dir);
    }
}


double ModelEvidenceCalculator::numerically_integrate(ArrayXd xs, ArrayXd ys) const
{
  if (xs.size() != ys.size()) {
    throw std::invalid_argument("For numeric integration x- and y-values must be same length.");
  }

  double area = 0;
  for (int i = 1; i < xs.size(); ++i) {
    double height = 0.5 * (ys(i) + ys(i-1));
    double width = xs(i) - xs(i-1);
    area += (height * width);
  }

  return area;
}

ArrayXd ModelEvidenceCalculator::calc_meanloglikelihoods() const
{
  VectorXd m = mcmc.loglikelihoods(Eigen::seq(mcmc.n_burn_iters, Eigen::last), Eigen::all).colwise().mean();
  return m.array();
}

ModelEvidenceCalculator::ModelEvidenceCalculator(const MCMC& mcmc,  ModelFit::SampleStatistics& sample_stats) 
    : mcmc(mcmc), sample_stats(sample_stats) {}

void ModelEvidenceCalculator::calc_logevidence()
{
    ArrayXd meanloglikelihoods = calc_meanloglikelihoods();
    ArrayXd betas = ArrayXd::Constant(mcmc.n_temps, -1.0);
    for (int k = 0; k < betas.size(); ++k) {
      betas[k] = mcmc.temps[k].beta;
    }

    // Calculate new x- and y-values using transformed TI method, assuming beta raised to the power 2
    ArrayXd betas_trans = pow(betas, 0.5);
    ArrayXd meanloglikelihoods_trans = 2 * betas_trans * meanloglikelihoods;
    sample_stats.logevidence = numerically_integrate(betas_trans, meanloglikelihoods_trans);
}






CompareAcrossCOI::CompareAcrossCOI(const std::string& output_dir) : output_dir(output_dir) 
{
    sample_stats_vec.reserve(10); // TODO: arbitrarily set to enough space
    posterior.reserve(10);
}

void CompareAcrossCOI::load_sample_statistics()
{
    // Find the files and store them
    std::regex pattern("K[1-9]{1}");
    for (const auto& entry : fs::directory_iterator(output_dir)) {
        
        // Skip if not a directory
        if (!entry.is_directory()) {
            continue;
        }

        // Skip if doesn't match the regex
        if (!std::regex_match(entry.path().filename().string(), pattern)) {
            continue;
        }

        std::string json_path = entry.path().string() + "/fit.sample_stats.json";
        ModelFit::SampleStatistics sample_stats = ModelFit::load_statistics_from_json(json_path);
        sample_stats_vec.emplace_back(sample_stats);
    }
    std::sort(sample_stats_vec.begin(), sample_stats_vec.end(), ModelFit::compare_by_coi);
}

//TODO: Clean up
void CompareAcrossCOI::calc_posterior()
{
    int n = sample_stats_vec.size();
    ArrayXd logevidences(n);
    for (int j = 0; j < n; ++j) {
        logevidences[j] = sample_stats_vec[j].logevidence;
    }
    ArrayXd scaled_logevidence = logevidences - logevidences.maxCoeff(); // prevent underflow
    ArrayXd posterior = scaled_logevidence.exp();
    posterior /= posterior.sum();
    for (int j = 0; j < n; ++j) {
        sample_stats_vec[j].model_posterior = posterior[j];
    }
}

void CompareAcrossCOI::write_comparison_csv()
{
    std::string csv = output_dir + "/compare.sample_summary.csv";
    std::ofstream csv_file(csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open file.");
    }
    csv_file << "n_sites,K,Keff,ws,f_ibd,l_ibd,n50_ibd,n_ibd,logposterior,aic,bic,logevidence,model_posterior" << "\n";
    for (ModelFit::SampleStatistics& s : sample_stats_vec) {
        csv_file << s.n_sites << ",";
        csv_file << s.K << ",";
        csv_file << s.Keff << ",";
        for (int k = 0; k < s.K - 1; ++k) {
            csv_file << s.ws[k] << ";";
        }
        csv_file << s.ws[s.K - 1] << ",";
        csv_file << s.f_ibd << ",";
        csv_file << s.l_ibd << ",";
        csv_file << s.n50_ibd << ",";
        csv_file << s.n_ibd << ",";
        csv_file << s.logposterior << ",";
        csv_file << s.aic << ",";
        csv_file << s.bic << ",";
        csv_file << s.logevidence << ",";
        csv_file << s.model_posterior << "\n";
    }
    csv_file.close();
}
