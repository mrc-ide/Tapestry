#pragma once
#include <vector>

// --------------------------------------------------------------------------------
// Solve label-switching problem
// --------------------------------------------------------------------------------

class label_switching
{
public:
    
    // FUNCTIONS
    // Constructor
    label_switching() {};

    void fix_labels(std::vector<Particle>& particle_trace);

    // Destructor
    ~label_switching() {}
};

