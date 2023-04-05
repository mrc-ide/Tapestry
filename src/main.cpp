#include <iostream>
#include <iomanip>
#include <string>
#include "libs/cli11/CLI11.hpp"
#include "data.hpp"
#include "ibd.hpp"
#include "io.hpp"
#include "mcmcs.hpp"
#include "mcmcs_summary.hpp"
#include "models.hpp"
#include "model_fitting.hpp"
#include "parameters.hpp"
#include "particle_writers.hpp"
#include "particles.hpp"
#include "proposals.hpp"
#include "timer.hpp"
using namespace std;


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
    int K = 2;                  // COI
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
                ->check(CLI::Range(1, 6));
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

        // Create Parameters object
        Parameters params(K, e_0, e_1, v, rho, w_proposal_sd, n_pi_bins);
        params.print();

        // Create Data object
        VCFData data(input_vcf, sample_name);
        data.print();

        // Create a proposal engine
        cout << "Creating a proposal engine..." << endl;
        ProposalEngine proposal_engine(params);

        // Create a Model
        cout << "Creating a model..." << endl;
        NaiveIBDModel model(params, data);  // switch NoIBD model
        model.print();

        // Create MCMC
        cout << "Creating an MCMC..." << endl;
        MetropolisHastings mcmc(params, model, proposal_engine);
        
        cout << "Running..." << endl;
        mcmc.run();
        
        cout << "Writing output..." << endl;
        ProportionParticleWriter particle_writer;
        mcmc.write_output(output_dir, particle_writer);

        // Fit
        // Procedure
        cout << "Fitting..." << endl;
        Particle map_particle = mcmc.get_map_particle();
        std::sort(map_particle.ws.begin(), map_particle.ws.end());
        ViterbiResult viterbi = model.get_viterbi_path(map_particle);

        // Object
        ModelFit model_fit(
            params,
            data,
            viterbi.logposterior, 
            map_particle.ws,
            viterbi.path
        );
        cout << "  Log-posterior: " << model_fit.logposterior << endl;
        cout << "  Proportions: " << model_fit.ws.transpose() << endl;
        cout << "  Effective COI: " << model_fit.Keff << endl;
        if (model_fit.has_ibd) {
            cout << "  No. IBD Segments: " << model_fit.n_ibd << endl;
            cout << "  Mean IBD Segment Length (kbp): " << model_fit.l_ibd/1000 << endl;
        }


        model_fit.write_output(output_dir);

        // Particle
        // IBDContainer ibd(params.K);
        // MCMCPointEstimator mcmc_estimator(mcmc, data, model, ibd);
        // mcmc_estimator.compute_point_estimate();
        // mcmc_estimator.write_output(output_dir);

    } else {
        // Throw an exception 
    }
    cout << string(80, '-') << endl;
    cout << "Time elapsed (ms): " << timer.elapsed<chrono::milliseconds>() << endl;
    cout << string(80, '-') << endl;
}