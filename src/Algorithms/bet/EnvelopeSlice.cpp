#include "Algorithms/bet/EnvelopeSlice.h"
#include <Physics/Physics.h>

EnvelopeSlice::EnvelopeSlice() {

    for(auto &slice_param : p) slice_param = 0.0;
    for(auto &slice_param : p_old) slice_param = 0.0;

    double initial_gamma  = 1.0 + (1.0e6 * Physics::q_e /
                           (Physics::EMASS * Physics::c * Physics::c));
    p[SLI_beta] = sqrt(1.0 / (1.0 - (1.0 / (initial_gamma * initial_gamma))));
    p[SLI_x]    = 1.0e-3;
    p[SLI_y]    = 1.0e-3;
    isValid_m   = 1;

    hasSliceBeenEmitted_m = false;
    backup();
}

void EnvelopeSlice::backup() {

    backupedValid_m= isValid_m;
    for(auto i = 0; i < SLNPAR; i++)
        p_old[i] = p[i];
}

void EnvelopeSlice::restore() {

    isValid_m = backupedValid_m;
    for(auto i = 0; i < SLNPAR; i++)
        p[i] = p_old[i];
}

int EnvelopeSlice::check() {

    int changed = 0;
    if(isValid_m) {
        isValid_m = (p[SLI_beta] > 0.0);
        if(!isValid_m)
            changed = 1;
    }
    return changed;
}

