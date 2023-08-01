#include <iostream>
#include <iomanip>
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
using namespace std; // let's remove this soon


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
    double e_0 = 0.01;          // REF -> ALT error probability
    double e_1 = 0.05;          // ALT -> REF error probability
    double v = 100;             // WSAF dispersion
    double rho = 13.5;         // Recombination rate; kbp per cM
    int n_pi_bins = 1000;       // No. of WSAF bins in betabinomial lookup
    
    // MCMC parameters
    double w_proposal_sd = 0.5;  // Titres sampled from ~N(0, w_propsal_sd)

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
    cmd_infer->add_option("-b, --n_wsaf_bins", n_pi_bins, "Number of WSAF bins in Betabin lookup table.")
                ->group("Model Hyperparameters")
                ->check(CLI::Range(100, 10'000));
    cmd_infer->add_option("-w, --w_proposal", w_proposal_sd, "Controls variance in proportion proposals.")
                ->group("MCMC Parameters")
                ->check(CLI::PositiveNumber);

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
        vector<int> Ks(1, K); 
        if (K == -1) { // iterate over multiple K values
            Ks.resize(max_K - min_K + 1);
            std::iota(Ks.begin(), Ks.end(), min_K);
        }

        // We will store the model fits for comparison later
        // std::vector<std::unique_ptr<ModelFits>> model_fits; TODO: figure out how to do this properly
        vector<ModelFit> model_fits;
        model_fits.reserve(max_K - min_K + 1);
        //model_fits.reserve(Ks.size());
        for (int k : Ks) {
            // Define output directory
            string K_output_dir = output_dir + "/K" + std::to_string(k);

            // Create key objects
            Parameters params(k, e_0, e_1, v, rho, w_proposal_sd, n_pi_bins);
            ProposalEngine proposal_engine(params);
            NaiveIBDModel model(params, data); // TODO: Stop recreating BetabinArray 
            MetropolisHastings mcmc(params, model, proposal_engine);

            // Some information to usre
            params.print();

            // Run MCMC
            cout << "Runnning MCMC..." << endl;
            mcmc.run();

            // Write MCMC outputs
            cout << "Writing MCMC outputs..." << endl;
            ProportionParticleWriter particle_writer;
            mcmc.write_output(K_output_dir, particle_writer);

            // Fitting
            cout << "Fitting..." << endl;
            Particle map_particle = mcmc.get_map_particle();
            std::sort(map_particle.ws.begin(), map_particle.ws.end());
            ViterbiResult viterbi = model.get_viterbi_path(map_particle);
            // The above gets MAP information, but still need to create
            // output files that are human-readable, and compute interesting
            // summary statistics
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

        // Model comparison
        cout << "Comparing fits across models..." << endl;
        ModelCompare model_compare(
            model_fits,
            data.n_sites,
            Ks
        );
        string output_csv = output_dir + "/compare.heuristics.csv";
        model_compare.write_output(output_csv);
        cout << "Done." << endl;
    } else {
        // Throw an exception 
    }
    cout << string(80, '-') << endl;
    cout << "Time elapsed (ms): " << timer.elapsed<chrono::milliseconds>() << endl;
    cout << string(80, '-') << endl;
}