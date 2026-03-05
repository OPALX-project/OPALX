//
// Random number generator to replace GSL RNG
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#ifndef OPAL_RANDOM_HH
#define OPAL_RANDOM_HH

#include <random>
#include <memory>
#include <cstddef>

// GSL-compatible random number generator wrapper
class gsl_rng {
public:
    gsl_rng() : engine_(std::mt19937_64(std::random_device{}())) {}
    
    double uniform() {
        return uniform_dist_(engine_);
    }
    
    double gaussian(double sigma) {
        std::normal_distribution<double> dist(0.0, sigma);
        return dist(engine_);
    }
    
    void seed(unsigned long seed) {
        engine_.seed(seed);
    }
    
    std::mt19937_64& engine() { return engine_; }
    const std::mt19937_64& engine() const { return engine_; }
    
private:
    std::mt19937_64 engine_;
    std::uniform_real_distribution<double> uniform_dist_{0.0, 1.0};
};

// GSL RNG type (not used, but kept for compatibility)
struct gsl_rng_type {
    const char* name;
};

// Default RNG type
extern const gsl_rng_type* gsl_rng_default;

// GSL-compatible functions
inline void gsl_rng_env_setup() {
    // No-op: environment setup not needed with std::random
}

inline gsl_rng* gsl_rng_alloc(const gsl_rng_type* /* type */) {
    return new gsl_rng();
}

inline void gsl_rng_free(gsl_rng* rng) {
    delete rng;
}

inline double gsl_rng_uniform(gsl_rng* rng) {
    return rng->uniform();
}

inline double gsl_ran_gaussian(gsl_rng* rng, double sigma) {
    return rng->gaussian(sigma);
}

// Dummy implementation (not used in practice)
inline const gsl_rng_type* get_gsl_rng_default() {
    static const gsl_rng_type dummy = {"mt19937"};
    return &dummy;
}
#define gsl_rng_default (get_gsl_rng_default())

#endif // OPAL_RANDOM_HH

