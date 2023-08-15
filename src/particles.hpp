#pragma once
#include "libs/eigen-3.4.0/Dense"
using Eigen::RowVectorXd;  // TODO: I probably want RowVector


struct Particle
{
    // Members
    // double loglikelihood;           // TODO: REMOVE?
    // double logprior;                // TODO: REMOVE?
    RowVectorXd ws;                 // Strain proportions

    // Constructors
    Particle();                     // TODO: Init. empty; helpful when preallocating space, I think?
    Particle(int K);                // Init. with all equal
    Particle(RowVectorXd ws);       // Init. with passed values

    // TODO:
    // - Should often be referenced via pointer -- do I need a ~Particle();?   
};