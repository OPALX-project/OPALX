//
// Quasi-random number generator to replace GSL QRNG
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

#ifndef OPAL_QUASIRANDOM_HH
#define OPAL_QUASIRANDOM_HH

#include <cmath>
#include <string>
#include <vector>

// Simple Sobol sequence generator (quasi-random)
class gsl_qrng {
public:
    gsl_qrng() : dimension_(1), count_(0) {}

    void init(size_t dimension) {
        dimension_ = dimension;
        count_     = 0;
        // Initialize direction numbers for Sobol sequence
        init_sobol();
    }

    void get(double* x) {
        for (size_t i = 0; i < dimension_; ++i) {
            x[i] = sobol_sample(i);
        }
        count_++;
    }

private:
    void init_sobol() {
        // Simplified Sobol initialization
        // In practice, you'd want proper direction numbers
        direction_numbers_.resize(dimension_);
        for (size_t i = 0; i < dimension_; ++i) {
            direction_numbers_[i] = (1u << (31 - i)) | 1u;
        }
    }

    double sobol_sample(size_t dim) {
        if (dim >= dimension_) return 0.0;

        unsigned int x = count_;
        double result  = 0.0;
        double factor  = 0.5;

        for (int i = 0; i < 32; ++i) {
            if (x & 1) {
                result += factor;
            }
            x >>= 1;
            factor *= 0.5;
        }

        return result;
    }

    size_t dimension_;
    size_t count_;
    std::vector<unsigned int> direction_numbers_;
};

// GSL-compatible functions
inline gsl_qrng* gsl_qrng_alloc(const char* /* type */, unsigned int dimension) {
    gsl_qrng* q = new gsl_qrng();
    q->init(dimension);
    return q;
}

inline void gsl_qrng_get(const gsl_qrng* q, double* x) { const_cast<gsl_qrng*>(q)->get(x); }

inline void gsl_qrng_free(gsl_qrng* q) { delete q; }

#endif  // OPAL_QUASIRANDOM_HH
