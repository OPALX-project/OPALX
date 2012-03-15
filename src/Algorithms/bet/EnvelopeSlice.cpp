#include <Physics/Physics.h>
#include "Algorithms/bet/EnvelopeSlice.h"

EnvelopeSlice::EnvelopeSlice() {
    for(int i = 0; i < SLNPAR; i++) {
        p_old[i] = 0.0;
        p[i] = 0.0;
    }

    double g  = 1.0 + (1.0e6 * Physics::q_e / (Physics::EMASS * Physics::c * Physics::c));
    p[SLI_beta] = sqrt(1.0 / (1.0 - (1.0 / (g * g))));
    p[SLI_x]    = 1.0e-3;
    p[SLI_y]    = 1.0e-3;
    valid       = 1;

    backup();
}

void EnvelopeSlice::backup() {
    was_valid = valid;
    for(int i = 0; i < SLNPAR; i++)
        p_old[i] = p[i];
}

void EnvelopeSlice::restore() {
    valid = was_valid;
    for(int i = 0; i < SLNPAR; i++)
        p[i] = p_old[i];
}

int EnvelopeSlice::check() {
    int changed = 0;

    if(valid) {
        valid = (p[SLI_beta] > 0.0);
        if(!valid)
            changed = 1;
    }
    return changed;
}

