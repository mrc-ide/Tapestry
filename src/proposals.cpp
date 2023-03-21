#include <cmath>
#include <random>
#include "proposals.hpp"
#include "particles.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::RowVectorXd;


ProposalEngine::ProposalEngine(const Parameters& params)
    : params(params),
    alpha(1.0),
    beta(1.0),
    gamma_dist(alpha, beta),
    unif_dist(0, params.K - 1),
    norm_dist(0, params.w_proposal_sd)
{};


Particle ProposalEngine::create_particle()
{
    // Sample from gamma
    RowVectorXd ws(params.K);
    for (int i = 0; i < params.K; ++i) {
        ws(i) = gamma_dist(rng.engine);
    }

    // Normalise to produce dirichlet
    double sum = ws.sum();
    ws /= sum;
    
    // Return a Particle
    Particle particle(ws);
    return particle;
}


Particle ProposalEngine::propose_particle(const Particle& particle)
{
    // Sample a strain index to update
    int ix = unif_dist(rng.engine);

    // Convert to real numbers by log transforming
    // and update by adding a normal variate ~N(0, w_proposal_sd)
    RowVectorXd proposed_ws = particle.ws;
    double titre = log(proposed_ws(ix));  // TODO: what base should I use?
    titre += norm_dist(rng.engine);
    proposed_ws(ix) = exp(titre);

    // Normalise
    proposed_ws /= proposed_ws.sum();

    // Return a particle
    Particle proposed_particle(proposed_ws);
    return proposed_particle;
}

