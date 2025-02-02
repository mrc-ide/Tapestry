#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>
#include "libs/cli11/CLI11.hpp"
#include "libs/json.hpp"
#include "betabin.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "summarise.hpp"
#include "parameters.hpp"
#include "particle_writers.hpp"
#include "particles.hpp"
#include "phase.hpp"
#include "proposals.hpp"
#include "timer.hpp"
using Eigen::RowVectorXd;  // TODO: I probably want RowVector


int main(int argc, char* argv[])
{
    // PARSE CLI
    CLI::App app{"Infer COI and pairwise IBD for P. falciparum malaria"};

    // Subcommands
    // TODO: write a function to produce nice subcommand strings
    app.require_subcommand(1);
    CLI::App* cmd_filter = app.add_subcommand("filter", "Filter an input VCF prior to inference.");
    CLI::App* cmd_infer = app.add_subcommand("infer", "Run inference from an filtered VCF.");
    CLI::App* cmd_phase = app.add_subcommand("phase", "Phase genomes in a mixed infection, after inference is complete.");

    // DEFAULTS
    // IO
    std::string input_vcf;
    std::string sample_name;
    std::string output_dir = "tapestry_output";

    // Model hyperparameters
    int K;
    int minK = 1;
    int maxK = 4;
    double e_0 = 0.0001;            // REF -> ALT error probability
    double e_1 = 0.005;             // ALT -> REF error probability
    double v = 500;                 // WSAF dispersion
    double rho = 13.5;              // Recombination rate; kbp per cM
    double G = 5.0;                 // Generation parameter applied to IBD detection.
    int n_pi_bins = 1000;           // No. of WSAF bins in betabinomial lookup
    
    // MCMC parameters
    int n_burn_iters = 100;              // Number of iterations in burn-in phase
    int n_sample_iters = 900;            // Number of iterations in sampling phase
    int n_temps = 5;                     // Number of temperature levels for PT-MCMC
    double target_acceptance = 0.44;     // Target acceptance rate per MCMC rung
    int swap_freq = 1;                   // Number of iterations between proposing Metropolis-coupling swaps

    // OPTIONS
    // Filter
    cmd_filter->add_option("-i,--input_vcf", input_vcf, "Path to input VCF file.")
                ->group("Input and output")
                ->check(CLI::ExistingFile)
                ->required();

    // Infer
    cmd_infer->add_option("-i,--input_vcf", input_vcf, "Path to input VCF file.")
                ->group("Input and output")
                ->check(CLI::ExistingFile)
                ->required();
    cmd_infer->add_option("-s,--target_sample", sample_name, "Target sample in VCF.")
                ->group("Input and output")
                ->required();
    cmd_infer->add_option("-o,--output_dir", output_dir, "Output directory.")
                ->group("Input and output");
    cmd_infer->add_option("-k, --minK", minK, "Minimum COI.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(1, 9));
    cmd_infer->add_option("-K, --maxK", maxK, "Maximum COI.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(minK, 9));
    cmd_infer->add_option("-e, --error_ref", e_0, "Probability of REF->ALT error.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(0, 1));
    cmd_infer->add_option("-E, --error_alt", e_1, "Probability of ALT->REF error.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(0, 1));
    cmd_infer->add_option("-v, --var_wsaf", v, "Controls dispersion in WSAF. Larger is less dispersed.")
                ->group("Model Hyperparameters")
                ->check(CLI::PositiveNumber);
    cmd_infer->add_option("-r, --recomb_rate", rho, "Recombination rate in kbp/cM.")
                ->group("Model Hyperparameters")
                ->check(CLI::PositiveNumber);
    cmd_infer->add_option("-g, --gens", G, "Expected number of ancestral generations for IBD segment length.")
                ->group("Model Hyperparameters")
                ->check(CLI::PositiveNumber);
    cmd_infer->add_option("-b, --n_wsaf_bins", n_pi_bins, "Number of WSAF bins in Betabin lookup table.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(100, 10'000));
    cmd_infer->add_option("-t, --temps", n_temps, "Number of temperature levels in PT-MCMC.")
                ->group("MCMC Parameters")
                ->check(CLI::Range(5, 100));
    cmd_infer->add_option("-B, --burnin", n_burn_iters, "Number of burn-in iterations.")
                ->group("MCMC Parameters")
                ->check(CLI::PositiveNumber);
    cmd_infer->add_option("-S, --sampling", n_sample_iters, "Number of sampling iterations.")
                ->group("MCMC Parameters")
                ->check(CLI::PositiveNumber);

    // Phase
    cmd_phase->add_option("-o,--output_dir", output_dir, "Output directory from `tapestry infer`.")
                ->group("Input and output")
                ->required();
    cmd_phase->add_option("-K, --coi", K, "Inferred COI value for which to phase.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(1, 9))
                ->required();
    cmd_phase->add_option("-e, --error_ref", e_0, "Probability of REF->ALT error.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(0, 1));
    cmd_phase->add_option("-E, --error_alt", e_1, "Probability of ALT->REF error.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(0, 1));
    cmd_phase->add_option("-v, --var_wsaf", v, "Controls dispersion in WSAF. Larger is less dispersed.")
                ->group("Model Hyperparameters")
                ->check(CLI::PositiveNumber);


    // Parse
    CLI11_PARSE(app, argc, argv);

    // RUN
    Timer timer_total;
    nlohmann::json runtime;
    
    // Filter
    if (app.got_subcommand("filter")) {
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "Tapestry: Filtering VCF to informative bi-allelic SNPs" << std::endl;
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "Not yet implemented!" << std::endl;
    
    // Infer
    } else if (app.got_subcommand("infer")) {
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "Tapestry: Inferring COI, proportions and IBD" << std::endl;
        std::cout << std::string(80, '-') << std::endl;

        // Load data
        VCFData data(input_vcf, sample_name);
        data.print();

        // Make beta-binomial matrix
        Betabinomial::LookupMatrix betabin_lookup(data, n_pi_bins, e_0, e_1, v, false);

        for (int k = minK; k <= maxK; ++k) {
            std::cout << "Inferring under K = " << k << std::endl;
            std::string K_string = "K" + std::to_string(k);
            std::string K_output_dir = output_dir + "/" + K_string;
            Timer timer_coi;

            // Creation
            Parameters params(k, e_0, e_1, v, rho, G, n_pi_bins, target_acceptance, swap_freq);
            ProposalEngine proposal_engine(params);
            Model model(params, data, betabin_lookup);

            // Run MCMC
            std::cout << "  Runnning MCMC..." << std::endl;
            MCMC mcmc(params, model, proposal_engine, n_burn_iters, n_sample_iters, n_temps);
            mcmc.run();

            // Write MCMC outputs
            std::cout << "  Writing MCMC outputs..." << std::endl;
            ProportionParticleWriter particle_writer;
            mcmc.write_output(K_output_dir, particle_writer);

            // Fit
            std::cout << "  Fitting..." << std::endl;
            Particle map_particle = mcmc.get_map_particle();
            ViterbiResult viterbi = model.get_viterbi_path(map_particle);

            // Summarise
            ModelFit::SampleStatistics sample_stats;
            ModelFit::Summariser summariser(
                params, 
                data,
                sample_stats,
                viterbi.path, 
                map_particle.ws,
                viterbi.logposterior
            );
            summariser.write_output(K_output_dir);
            ModelEvidenceCalculator evidence_calculator(mcmc, sample_stats);
            evidence_calculator.calc_logevidence();
            ModelFit::write_statistics_to_json(sample_stats, K_output_dir);
            std::cout << "Done." << std::endl;

            // Store
            runtime[K_string] = timer_coi.elapsed<chrono::milliseconds>();
        }

        std::cout << "Comparing across COIs..." << std::endl;
        CompareAcrossCOI comparer(output_dir);
        comparer.load_sample_statistics();
        comparer.calc_posterior();
        comparer.write_comparison_csv();
        std::cout << "Done. Results are here: " << output_dir << std::endl;
    
    } else if (app.got_subcommand("phase")) {
        std::cout << std::string(80, '-') << std::endl;
        std::cout << "Tapestry: Phasing genomes after inference" << std::endl;
        std::cout << std::string(80, '-') << std::endl;

        std::cout << "Phasing under K = " << K << std::endl;
        std::string K_string = "K" + std::to_string(K);
        std::string K_output_dir = output_dir + "/" + K_string;

        InferredIBDPathData data(K_output_dir + "/fit.ibd.path.csv");
        data.print();

        // Load sample statistics
        ModelFit::SampleStatistics sample_stats = ModelFit::load_statistics_from_json(
            K_output_dir + "/fit.sample_stats.json"
        );
        RowVectorXd ws = Eigen::Map<ArrayXd>(sample_stats.ws.data(), sample_stats.ws.size());

        PanelFreePhaser phaser(
            K, ws, data,
            e_0, e_1, v
        );
        phaser.phase();
        phaser.write_output(K_output_dir);
        
        std::cout << "Done. Results are here: " << K_output_dir << std::endl;
    }

    double total_elapsed = timer_total.elapsed<chrono::milliseconds>();
    runtime["TOTAL"] = total_elapsed;

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Time elapsed (ms): " << total_elapsed << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    // Write runtime
    // TODO: this will get overwritten as different commands are run
    std::ofstream o(output_dir + "/runtime.json"); 
    o << std::setw(4) << runtime << std::endl;
}