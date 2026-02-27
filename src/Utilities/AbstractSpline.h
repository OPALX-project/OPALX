//
// Cubic Spline Interpolation to replace GSL spline
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

#ifndef OPALX_ABSTRACTSPLINE_H
#define OPALX_ABSTRACTSPLINE_H

#include "OPALTypes.h"

/// \brief Base class for the linear and cubic interpolation spline classes
class AbstractSpline {
public:
    /// \brief Default constructor.
    AbstractSpline() = default;

    /// \brief Polymorphic destructor
    virtual ~AbstractSpline() = default;

    // Accelerator for repeated evaluations (compatible with GSL interface)
    /// \brief Accelerator caching last interval indices.
    class Accelerator {
    public:
        Accelerator() = default;
        size_t last_index_{};
        size_t last_upper_index_{};
        void reset() { last_index_ = 0; last_upper_index_ = 0; }
    };

    /// \brief Initialize from tabulated data (natural spline).
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    virtual void init(const std::vector<double>& x, const std::vector<double>& y);

    /// \brief Evaluate the spline at x using an accelerator.
    /// \param x Input: x-coordinate.
    /// \param accel Input/Output: cached interval index.
    /// \return Output: interpolated value.
    virtual double eval(double x, Accelerator& accel) const = 0;

    /// \brief Evaluate the integral of the spline at \p x.
    /// \param xa Input: x-coordinate lower bound.
    /// \param xb Input: x-coordinate upper bound.
    /// \param accel Input/output: accelerator for interval caching.
    /// \return Output: interpolated integral value (with linear extrapolation).
    virtual double evalIntegral(double xa, double xb, Accelerator& accel) const = 0;

protected:

    /// \brief  Return the interval index for the given x-coordinate, using the
    /// supplied cached value if possible
    /// \param x Input: x-coordinate.
    /// \param intervalCache Input/output: The cached interval index.
    /// \return Output: The index of the interval containing x.
    size_t findInterval(double x, size_t& intervalCache) const;

    // Members
    std::vector<double> x_{};
    std::vector<double> y_{};
};


#endif //OPALX_ABSTRACTSPLINE_H
