#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>
#include "libs/cli11/CLI11.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "models.hpp"
#include "model_compare.hpp"
#include "model_fitting.hpp"
#include "parameters.hpp"
#include "particle_writers.hpp"
#include "particles.hpp"
#include "proposals.hpp"
#include "timer.hpp"
using namespace std; // TODO: let's remove this soon


int main(int argc, char* argv[])
{

    // PARSE CLI
    CLI::App app{"Infer COI and pairwise IBD for P. falciparum malaria"};

    // Subcommands
    // TODO: write a function to produce nice subcommand strings
    app.require_subcommand(1);
    CLI::App* cmd_filter = app.add_subcommand("filter", "Filter an input VCF prior to inference.");
    CLI::App* cmd_infer = app.add_subcommand("infer", "Run inference from an filtered VCF.");

    // DEFAULTS
    // IO
    string input_vcf;
    string sample_name;
    string output_dir = "tapestry_output";

    // Model hyperparameters
    // TODO: 
    // - Could I initialise Parameters here?
    // - Then options take params-><var>; &c
    // - However, would not be immutable
    int min_K = 1;
    int max_K = 4;
    int K = -1;
    double e_0 = 0.0001;            // REF -> ALT error probability
    double e_1 = 0.005;             // ALT -> REF error probability
    double v = 500;                 // WSAF dispersion
    double rho = 13.5;              // Recombination rate; kbp per cM
    double G = 5.0;                 // Generation parameter applied to IBD detection.
    int n_pi_bins = 1000;           // No. of WSAF bins in betabinomial lookup
    
    // MCMC parameters
    double w_proposal_sd = 0.1;     // Titres sampled from ~N(0, w_propsal_sd)
    int n_temps = 5;                // Number of temperature levels for PT-MCMC

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
    cmd_infer->add_option("-K, --COI", K, "Complexity of infection.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(min_K, max_K));
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
    cmd_infer->add_option("-w, --w_proposal", w_proposal_sd, "Controls variance in proportion proposals.")
                ->group("MCMC Parameters")
                ->check(CLI::PositiveNumber);
    cmd_infer->add_option("-t, --temps", n_temps, "Number of temperature levels in PT-MCMC.")
                ->group("MCMC Parameters")
                ->check(CLI::Range(5, 20));

    // Parse
    CLI11_PARSE(app, argc, argv);

    // RUN
    // Start timing...
    Timer timer;
    
    // Filter
    if (app.got_subcommand("filter")) {
        cout << string(80, '-') << endl;
        cout << "Tapestry: Filtering VCF to informative bi-allelic SNPs" << endl;
        cout << string(80, '-') << endl;
        cout << "Not yet implemented!" << endl;
    
    // Infer
    } else if (app.got_subcommand("infer")) {
        cout << string(80, '-') << endl;
        cout << "Tapestry: Inferring COI, proportions and IBD" << endl;
        cout << string(80, '-') << endl;

        // Load data
        VCFData data(input_vcf, sample_name);
        data.print();

        // Setup of Ks vector
        std::vector<int> Ks(1, K); 
        if (K == -1) { // iterate over multiple K values
            Ks.resize(max_K - min_K + 1);
            std::iota(Ks.begin(), Ks.end(), min_K);
        }

        // We will store the model fits for comparison later
        // std::vector<std::unique_ptr<ModelFits>> model_fits; TODO: figure out how to do this properly
        // TODO: Alternate implementation we store just a struct of the key things
        std::vector<ModelFit> model_fits;
        model_fits.reserve(Ks.size());
        std::vector<std::unique_ptr<ParallelTempering>> mcmc_ptrs;
        mcmc_ptrs.reserve(Ks.size());
        for (int k : Ks) {

            // Define output directory
            string K_output_dir = output_dir + "/K" + std::to_string(k);

            // Create objects for this COI
            Parameters params(k, e_0, e_1, v, rho, G, w_proposal_sd, n_pi_bins);
            ProposalEngine proposal_engine(params);
            NaiveIBDModel model(params, data); // TODO: Stop recreating BetabinArray 
            model.print();

            // Create MCMC on the heap
            cout << "Runnning MCMC..." << endl;
            mcmc_ptrs.emplace_back(std::make_unique<ParallelTempering>(params, model, proposal_engine, n_temps));
            mcmc_ptrs.back()->run();

            // Write MCMC outputs
            cout << "Writing MCMC outputs..." << endl;
            ProportionParticleWriter particle_writer;
            mcmc_ptrs.back()->write_output(K_output_dir, particle_writer);

            // Fitting
            cout << "Fitting..." << endl;
            Particle map_particle = mcmc_ptrs.back()->get_map_particle();
            std::sort(map_particle.ws.begin(), map_particle.ws.end());
            ViterbiResult viterbi = model.get_viterbi_path(map_particle);

            ModelFit model_fit(
                params,
                data,
                viterbi.logposterior, 
                map_particle.ws,
                viterbi.path
            );
            cout << "Writing fit outputs..." << endl;
            model_fit.write_output(K_output_dir);
            cout << "Done." << endl;

            // Store the fits
            model_fits.push_back(model_fit);
        }

        // Heuristic model comparison
        cout << "Comparing fits across models..." << endl;
        ModelCompare model_compare(
            model_fits,
            data.n_sites,
            Ks
        );
        string output_csv = output_dir + "/compare.heuristics.csv";
        model_compare.write_output(output_csv);
        cout << "Done." << endl;

        // Compute model evidence
        std::cout << "Computing model evidence..." << std::endl;
        ModelEvidence model_evidence(mcmc_ptrs);
        model_evidence.calc_summary();
        std::string evidence_csv = output_dir + "/compare.evidence.csv";
        model_evidence.write_output(evidence_csv);
        std::cout << "Done." << std::endl;

    } else {
        // Throw an exception 
    }
    cout << string(80, '-') << endl;
    cout << "Time elapsed (ms): " << timer.elapsed<chrono::milliseconds>() << endl;
    cout << string(80, '-') << endl;
}