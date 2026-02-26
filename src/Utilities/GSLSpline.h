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

#ifndef OPALX_GSL_SPLINE_H
#define OPALX_GSL_SPLINE_H

#include "CubicSpline.h"
#include "LinearSpline.h"

// Compatibility layer for GSL-like interface

/// \brief GSL linear interpolation type identifiers.
constexpr int gsl_interp_linear  = 1;
constexpr int gsl_interp_cspline = 0;

/// \brief GSL-compatible spline type alias.
using gsl_spline = AbstractSpline;

/// \brief GSL-compatible accelerator type alias.
using gsl_interp_accel = AbstractSpline::Accelerator;

/// \brief Allocate a spline instance (size ignored).
/// \param type Input: interpolation type
/// \param size Input: number of points (unused).
/// \return Output: spline pointer.
inline gsl_spline* gsl_spline_alloc(const int type, size_t /*size*/) {
    gsl_spline* result{};
    if (type == gsl_interp_cspline) {
        result = new CubicSpline();
    } else if (type == gsl_interp_linear) {
        result = new LinearSpline();
    }
    return result;
}

/// \brief Allocate an interpolation accelerator.
/// \return Output: accelerator pointer.
inline gsl_interp_accel* gsl_interp_accel_alloc() { return new AbstractSpline::Accelerator(); }

/// \brief Initialize a spline with tabulated data.
/// \param spline Input/Output: spline to initialize.
/// \param x Input: strictly increasing x-coordinates.
/// \param y Input: y-values corresponding to \p x.
/// \param n Input: number of samples.
inline void gsl_spline_init(gsl_spline* spline, const double* x, const double* y, const size_t n) {
    spline->init(std::vector<double>{x, x + n}, std::vector<double>{y, y + n});
}

/// \brief Evaluate a spline at \p x using an accelerator.
/// \param spline Input: spline to evaluate.
/// \param x Input: x-coordinate.
/// \param accel Input/Output: accelerator cache.
/// \return Output: interpolated value.
inline double gsl_spline_eval(const gsl_spline* spline, const double x, gsl_interp_accel* accel) {
    return spline->eval(x, *accel);
}

/// \brief Evaluate the integral of a spline at \p x using an accelerator.
/// \param spline Input: spline to evaluate.
/// \param xa Input: lower bound x-coordinate.
/// \param xb Input: upper bound x-coordinate.
/// \param accel Input/Output: accelerator cache.
/// \return Output: integrated value.
inline double gsl_spline_eval_integ(
    const gsl_spline* spline, const double xa, const double xb, gsl_interp_accel* accel
) {
    return spline->evalIntegral(xa, xb, *accel);
}

/// \brief Free a spline instance.
/// \param spline Input: spline to release (can be null).
inline void gsl_spline_free(const gsl_spline* spline) { delete spline; }

/// \brief Free an accelerator instance.
/// \param accel Input: accelerator to release (can be null).
inline void gsl_interp_accel_free(const gsl_interp_accel* accel) { delete accel; }

/// \brief Reset an accelerator to the initial state.
/// \param accel Input/Output: accelerator to reset.
inline void gsl_interp_accel_reset(gsl_interp_accel* accel) {
    if (accel) {
        accel->reset();
    }
}

#endif  // OPALX_GSL_SPLINE_H