#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "particles.hpp"
#include "particle_writers.hpp"
#include "libs/eigen-3.4.0/Dense"


// ================================================================================
// Particle writer interface
//
// ================================================================================


ParticleWriter::ParticleWriter() {};
ParticleWriter::~ParticleWriter() {};


// ================================================================================
// Concrete implementations
//
// ================================================================================

// --------------------------------------------------------------------------------
// Write proportions
// --------------------------------------------------------------------------------


ProportionParticleWriter::ProportionParticleWriter() {};

void ProportionParticleWriter::write_particle_trace(
    const std::string& output_csv,
    const std::vector<Particle>& particle_trace
) const
{
    // Get key statistics
    int n_iters = particle_trace.size();
    int K = particle_trace[0].ws.size();

    // Open the output file
    std::ofstream csv_file(output_csv);
    if (!csv_file.is_open()) {
        throw std::invalid_argument("Could not open output file.");
    }

    // Write column names
    csv_file << "iter";
    for (int k = 1; k <= K; ++k) {
        csv_file << ",prop" << k;
    }
    csv_file << "\n";

    // Format proportions
    const static Eigen::IOFormat CSVFormat(
        5,
        Eigen::DontAlignCols, 
        ", ", 
        "\n",
        "",
        "\n"
        );

    // Write proportion values, including iteration
    for (int i = 0; i < n_iters; ++i) {
        csv_file << i << "," << particle_trace[i].ws.format(CSVFormat);
    }

    csv_file.close();
}