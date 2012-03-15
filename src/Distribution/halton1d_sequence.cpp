#include <cmath>
#include <iostream>
#include "halton1d_sequence.hh"

#define DIM 1

typedef struct
{
    unsigned int sequence_count;
    unsigned int call_count;
    double x[DIM];
} 
halton_state_t;

static void
halton1d_set(void * vstate, unsigned long int s) {
    halton_state_t *h_state = (halton_state_t *) vstate;
    
    h_state->sequence_count = 0;
    h_state->call_count = 0;
}

static double 
halton1d_get_double(void *vstate) {
    
    // double x[DIM];
    // gsl_qrng_halton->get(vstate, DIM, x);
    // return x[DIM-1];

    halton_state_t *h_state = (halton_state_t *) vstate;
    if (h_state->call_count == 0) {
        gsl_qrng_halton->get(vstate, DIM, h_state->x);
    }
    double x = h_state->x[h_state->call_count ++];
    if (h_state->call_count == DIM) h_state->call_count = 0;

    return x;
}

static unsigned long int 
halton1d_get(void *vstate) {
    return static_cast<unsigned long int>(halton1d_get_double(vstate) * 16777216.0);
}

static const gsl_rng_type halton1d_type = 
    {"halton1d",
     0x00ffffffUL,
     0,
     sizeof(halton_state_t),
     &halton1d_set,
     &halton1d_get,
     &halton1d_get_double};

const gsl_rng_type *gsl_rng_1dhalton = &halton1d_type;

