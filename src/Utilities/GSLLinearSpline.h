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

#ifndef OPALX_LINEAR_SPLINE_H
#define OPALX_LINEAR_SPLINE_H

#include <vector>
#include <algorithm>

// \brief Linear spline interpolation class that implements an equivelent of
// the gsl_interp_linear mode of the GSL spline class
class LinearSpline {
public:
    /// \brief Default constructor.
    LinearSpline() = default;

    /// \brief Construct and initialize from tabulated data.
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    LinearSpline(const std::vector<double>& x, const std::vector<double>& y) {
        init(x, y);
    }

    /// \brief Initialize from tabulated data (natural spline).
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    void init(const std::vector<double>& x, const std::vector<double>& y) {
        // Validate and save the parameters
        if (x.size() != y.size()) {
            throw std::invalid_argument("LinearSpline: x and y must be the same size");
        }
        if (x.size() < 2) {
            throw std::invalid_argument("LinearSpline: need at least 2 points");
        }
        for (size_t i = 1; i < x.size(); ++i) {
            if (x[i] <= x[i - 1]) {
                throw std::invalid_argument("LinearSpline: x must be strictly increasing");
            }
        }
        x_ = x;
        y_ = y;
        // The gradients of each interval
        computeCoefficients();
        // And the integrals over each interval
        computeIntegrals();
    }

    // Accelerator for repeated evaluations (compatible with GSL interface)
    /// \brief Accelerator caching last interval index.
    class Accelerator {
    public:
        Accelerator() = default;
        size_t last_index_{};
        size_t last_upper_index_{};
    };

    /// \brief Evaluate the spline at x using an accelerator.
    /// \param x Input: x-coordinate.
    /// \param accel Input/Output: cached interval index.
    /// \return Output: interpolated value.
    double eval(const double x, Accelerator& accel) const {
        double result{};
        if (x_.empty()) {
            throw std::runtime_error("LinearSpline: not initialized");
        }
        if (x <= x_.front()) {
            result = extrapolateLeft(x);
        } else if (x >= x_.back()) {
            result = extrapolateRight(x);
        } else {
            const size_t i = findInterval(x, accel.last_index_);
            result = m_[i] * x + c_[i];
        }
        return result;
    }

    /// \brief Evaluate the integral of the spline at \p x.
    /// \param xa Input: x-coordinate lower bound.
    /// \param xb Input: x-coordinate upper bound.
    /// \param accel Input/output: accelerator for interval caching.
    /// \return Output: interpolated integral value (with linear extrapolation).
    double evalIntegral(double xa, double xb, Accelerator& accel) const {
        double result{};
        // Verify parameters
        if (x_.empty()) {
            throw std::runtime_error("LinearSpline: not initialized");
        }
        if (xb < xa) {
            std::swap(xa, xb);
        }
        // Linear extrapolation below the defined intervals
        if (xa < x_.front()) {
            const auto y  = extrapolateLeft(xa);
            const auto dx = x_.front() - xa;
            const auto dy = y_.front() - y;
            result += dx * y_.front() - dx * dy / 2.0;
            xa = x_.front();
        }
        // Linear extrapolation above the defined intervals
        if (xb > x_.back()) {
            const auto y  = extrapolateRight(xb);
            const auto dx = xb - x_.back();
            const auto dy = y - y_.back();
            result += dx * y_.back() + dx * dy / 2.0;
            xb = x_.back();
        }
        // Find the intervals
        const size_t ia = findInterval(xa, accel.last_index_);
        const size_t ib = findInterval(xb, accel.last_upper_index_);
        // The first partial interval
        result -= integral(ia, xa - x_[ia]);
        // All the full interval integrals
        for (size_t m = ia; m < ib; ++m) {
            result += integrals_[m];
        }
        // The last partial interval
        result += integral(ib, xb - x_[ib]);
        return result;
    }

private:
    /// \brief Compute linear spline coefficients.
    void computeCoefficients() {
        m_.resize(x_.size() - 1);
        c_.resize(x_.size() - 1);
        for (size_t i = 0; i < x_.size() - 1; ++i) {
            m_[i] = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
            c_[i] = y_[i] - m_[i] * x_[i];
        }
    }

    /// \brief Compute the integrals of the intervals.
    void computeIntegrals() {
        integrals_.resize(x_.size() - 1);
        for (size_t i = 0; i < x_.size() - 1; ++i) {
            integrals_[i] = integral(i, x_[i + 1] - x_[i]);
        }
    }

    /// \brief  Return the interval index for the given x-coordinate, using the
    /// accelerator if possible
    /// \param x Input: x-coordinate.
    /// \param intervalCache Input/output: The cached interval index.
    /// \return Output: The index of the interval containing x.
    size_t findInterval(const double x, size_t& intervalCache) const {
        if (intervalCache < x_.size() - 1 && x >= x_[intervalCache] && x < x_[intervalCache + 1]) {
            // x is in the cached interval
        } else {
            // Find the new interval
            const auto it = std::ranges::lower_bound(x_, x);
            intervalCache = it == x_.begin() ? 0 : std::distance(x_.begin(), it) - 1;
        }
        return intervalCache;
    }

    /// \brief Calculate the integral from x_[i] to x_[i]+dx
    /// \param i Input: the interval index
    /// \param dx Input: the distance in the interval over which to calculate
    /// \return Output: The integral.
    double integral(const size_t i, const double dx) const {
        const auto x2 = x_[i] + dx;
        return m_[i] * (x2 * x2 - x_[i] * x_[i]) / 2.0 + c_[i] * (x2 - x_[i]);
    }

    /// \brief Linear extrapolation to the left of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    double extrapolateLeft(const double x) const {
        return m_.front() * x + c_.front();
    }

    /// \brief Linear extrapolation to the right of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    double extrapolateRight(const double x) const {
        return m_.back() * x + c_.back();
    }

    // Members
    std::vector<double> x_{};
    std::vector<double> y_{};
    std::vector<double> m_{};  // Gradients
    std::vector<double> c_{};  // Constants
    std::vector<double> integrals_{};
};

#endif  // OPALX_LINEAR_SPLINE_H
