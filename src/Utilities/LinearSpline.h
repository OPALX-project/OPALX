//
// Linear Spline Interpolation to replace GSL spline
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

#ifndef OPALX_LINEAR_SPLINE_H
#define OPALX_LINEAR_SPLINE_H

#include "AbstractSpline.h"
#include <vector>

// \brief Linear spline interpolation class that implements an equivelent of
// the gsl_interp_linear mode of the GSL spline class
class LinearSpline final : public AbstractSpline {
public:
    /// \brief Default constructor.
    LinearSpline() = default;

    /// \brief Construct and initialize from tabulated data.
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    LinearSpline(const std::vector<double>& x, const std::vector<double>& y);

    /// \brief Initialize from tabulated data (natural spline).
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    void init(const std::vector<double>& x, const std::vector<double>& y) override;

    /// \brief Evaluate the spline at x using an accelerator.
    /// \param x Input: x-coordinate.
    /// \param accel Input/Output: cached interval index.
    /// \return Output: interpolated value.
    double eval(double x, Accelerator& accel) const override;

    /// \brief Evaluate the integral of the spline at \p x.
    /// \param xa Input: x-coordinate lower bound.
    /// \param xb Input: x-coordinate upper bound.
    /// \param accel Input/output: accelerator for interval caching.
    /// \return Output: interpolated integral value (with linear extrapolation).
    double evalIntegral(double xa, double xb, Accelerator& accel) const override;

private:
    /// \brief Compute linear spline coefficients.
    void computeCoefficients();

    /// \brief Compute the integrals of the intervals.
    void computeIntegrals();

    /// \brief Calculate the integral from x_[i] to x_[i]+dx
    /// \param i Input: the interval index
    /// \param dx Input: the distance in the interval over which to calculate
    /// \return Output: The integral.
    [[nodiscard]] double integral(size_t i, double dx) const;

    /// \brief Linear extrapolation to the left of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    [[nodiscard]] double extrapolateLeft(double x) const;

    /// \brief Linear extrapolation to the right of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    [[nodiscard]] double extrapolateRight(double x) const;

    // Members
    std::vector<double> m_{};  // Gradients
    std::vector<double> c_{};  // Constants
    std::vector<double> integrals_{};
};

#endif  // OPALX_LINEAR_SPLINE_H
