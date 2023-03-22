#pragma once
#include "parameters.hpp"
#include "particles.hpp"
#include "random.hpp"


class ProposalEngine
{
private:
    RNG rng;

    const Parameters& params;

    // For particle creation
    const double alpha;
    const double beta;
    std::gamma_distribution<double> gamma_dist;  // TODO: can I const prob dists?

    // For particle updating
    std::uniform_int_distribution<int> unif_dist;
    std::normal_distribution<double> norm_dist;

public:
    // Constructors
    ProposalEngine(const Parameters& params);

    /*
    * Create a new particle initialised with random values
    */
    Particle create_particle();

    /*
    * Propose a new particle based on the values in an existing
    * particle
    * TODO:
    * - Best would be for this to take a *pointer* to the proposed particle
    * - Then we don't copy, but update in place, I think
    */
    Particle propose_particle(const Particle& particle);
};