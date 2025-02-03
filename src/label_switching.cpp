#include <iostream>
#include <vector>
#include "particles.hpp"
#include "label_switching.hpp"

// --------------------------------------------------------------------------------
// Solve label-switching problem
// --------------------------------------------------------------------------------

// Fix order of all parameters in particles trace
void label_switching::fix_labels(std::vector<Particle>& particle_trace)
{
    int k = particle_trace[0].ws.size();
    for (size_t i = 0; i < particle_trace.size(); i++) {
        std::sort(particle_trace[i].ws.data(), particle_trace[i].ws.data() + k);
    }
}
