#include <cmath>
#include <random>
#include "proposals.hpp"
#include "particles.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::RowVectorXd;

// --------------------------------------------------------------------------------
// Propose w from reflected normal
// --------------------------------------------------------------------------------

ProposalEngine::ProposalEngine(const Parameters& params)
    : params(params),
    gamma_dist(1.0, 1.0),
    unif_dist(0, params.K - 1),
    unif_dist2(0, params.K - 2),
    norm_dist(0, 1)
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


Particle ProposalEngine::propose_particle(const Particle& particle, const double w_prop_sd)
{
    
    // Return if K = 1
    if (params.K == 1) {
        return particle;
    }

    // Sample the index of two strains without replacement
    int i1 = unif_dist(rng.engine);
    int i2 = unif_dist2(rng.engine);
    if (i2 >= i1) {
        i2++;
    }

    // Draw a new value for first weight from reflected normal, and calculate implied value of second weight
    double w1 = particle.ws[i1];
    double w2 = particle.ws[i2];
    double w_sum = w1 + w2;
    double w1_prop = rnorm_interval(w1, w_prop_sd, 0.0, w_sum);
    double w2_prop = w_sum - w1_prop;

    // Save values and return a particle
    RowVectorXd proposed_ws = particle.ws;
    proposed_ws[i1] = w1_prop;
    proposed_ws[i2] = w2_prop;

    Particle proposed_particle(proposed_ws);
    return proposed_particle;
}

double ProposalEngine::rnorm_interval(const double mean, const double sd, const double a, const double b) {

    // draw raw value relative to a
    double ret = sd*norm_dist(rng.engine) + mean - a;

    // reflect off boundries at 0 and (b-a)
    if (ret < 0 || ret > (b-a)) {

        // use multiple reflections to bring into range [-(b-a), 2(b-a)]
        if (ret < -(b - a)) {
            int n_double_intervals = floor(-ret / (b - a)) / 2;
            ret += 2 * (b - a) * (n_double_intervals + 1);
        } else if (ret > 2*(b - a)) {
            int n_double_intervals = floor(ret / (b - a) - 1) / 2;
            ret -= 2 * (b - a) * (n_double_intervals + 1);
        }

        // use one more reflection to bring into range [0, (b-a)]
        if (ret < 0) {
            ret = -ret;
        }
        if (ret > (b-a)) {
            ret = 2*(b-a) - ret;
        }
    }

    // no longer relative to a
    ret += a;

    return ret;
}

// --------------------------------------------------------------------------------
// Propose w from titre model
// --------------------------------------------------------------------------------

ProposalEngine_titre::ProposalEngine_titre(const Parameters& params)
    : params(params),
    alpha(1.0),
    beta(1.0),
    gamma_dist(alpha, beta),
    unif_dist(0, params.K - 1),
    norm_dist(0, 1.0)
{};


Particle ProposalEngine_titre::create_particle()
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


Particle ProposalEngine_titre::propose_particle(const Particle& particle, const double w_prop_sd)
{
    // Sample a strain index to update
    int ix = unif_dist(rng.engine);

    // Convert to real numbers by log transforming
    // and update by adding a normal variate ~N(0, w_proposal_sd)
    RowVectorXd proposed_ws = particle.ws;
    double titre = log(proposed_ws(ix));
    titre += w_prop_sd * norm_dist(rng.engine);
    proposed_ws(ix) = exp(titre);

    // Normalise
    proposed_ws /= proposed_ws.sum();

    // Return a particle
    Particle proposed_particle(proposed_ws);
    return proposed_particle;
}