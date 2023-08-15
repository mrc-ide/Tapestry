#include "particles.hpp"
#include "libs/eigen-3.4.0/Dense"
using Eigen::RowVectorXd;


Particle::Particle()
{};

Particle::Particle(int K)
    : ws(RowVectorXd::Constant(K, 1.0/K))
{};

Particle::Particle(RowVectorXd ws)
    : ws(ws)
{};

