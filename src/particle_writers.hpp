#pragma once
#include <string>
#include <vector>
#include "particles.hpp"


// ================================================================================
// Particle writer interface
//
// ================================================================================


/*  Write relevant information from a vector of Particle objects
*   to an output_csv
*
*   Pre: Must recieve Particle objects of the approporiate type
*   Post: Output CSV will exist
*   NB: This is basically just a function, for now.
*/
class ParticleWriter
{
public:
	ParticleWriter();

	virtual void write_particle_trace(
        const std::string& output_csv,
        const std::vector<Particle>& particle_trace
    ) const = 0;

    ~ParticleWriter();
};


// ================================================================================
// Concrete implementations
//
// ================================================================================

// --------------------------------------------------------------------------------
// Write proportions
// --------------------------------------------------------------------------------


class ProportionParticleWriter : public ParticleWriter
{
public:
    ProportionParticleWriter();
    void write_particle_trace(
        const std::string& output_csv,
        const std::vector<Particle>& particle_trace
    ) const override;
};

