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

#ifndef OPAL_CUBIC_SPLINE_HH
#define OPAL_CUBIC_SPLINE_HH

#include <vector>
#include <algorithm>
#include <stdexcept>

/// \see https://www.gnu.org/software/gsl/doc/html/interp.html
/// \see https://www.gnu.org/software/gsl/
// Simple cubic spline interpolation class to replace GSL spline
/// \brief Natural cubic spline interpolator.
class CubicSpline {
public:
    CubicSpline() : x_(), y_(), b_(), c_(), d_() {}
    
    /// \brief Construct and initialize from tabulated data.
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    /// \param n Input: number of samples (>= 2).
    CubicSpline(const double* x, const double* y, size_t n) {
        init(x, y, n);
    }
    
    /// \brief Initialize from tabulated data (natural spline).
    /// \param x Input: strictly increasing x-coordinates.
    /// \param y Input: y-values corresponding to \p x.
    /// \param n Input: number of samples (>= 2).
    void init(const double* x, const double* y, size_t n) {
        if (n < 2) {
            throw std::invalid_argument("CubicSpline: need at least 2 points");
        }
        
        x_.assign(x, x + n);
        y_.assign(y, y + n);
        
        // Check that x is strictly increasing
        for (size_t i = 1; i < n; ++i) {
            if (x_[i] <= x_[i-1]) {
                throw std::invalid_argument("CubicSpline: x must be strictly increasing");
            }
        }
        
        // Natural cubic spline: second derivatives at endpoints are zero
        computeCoefficients();
    }
    
    /// \brief Evaluate the spline at \p x.
    /// \param x Input: x-coordinate.
    /// \return Output: interpolated value (with linear extrapolation).
    double eval(double x) const {
        if (x_.empty()) {
            throw std::runtime_error("CubicSpline: not initialized");
        }
        
        // Handle extrapolation
        if (x <= x_[0]) {
            return extrapolateLeft(x);
        }
        if (x >= x_.back()) {
            return extrapolateRight(x);
        }
        
        // Find the interval containing x
        auto it = std::lower_bound(x_.begin(), x_.end(), x);
        size_t i = std::distance(x_.begin(), it) - 1;
        
        if (i >= x_.size() - 1) {
            i = x_.size() - 2;
        }
        
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }
    
    // Accelerator for repeated evaluations (compatible with GSL interface)
    /// \brief Accelerator caching last interval index.
    class Accelerator {
    public:
        Accelerator() : last_index_(0) {}
        size_t last_index_;
    };
    
    /// \brief Evaluate the spline at \p x using an accelerator.
    /// \param x Input: x-coordinate.
    /// \param accel Input/Output: cached interval index.
    /// \return Output: interpolated value.
    double eval(double x, Accelerator& accel) const {
        // Use cached index if possible
        if (accel.last_index_ < x_.size() - 1 && 
            x >= x_[accel.last_index_] && 
            x < x_[accel.last_index_ + 1]) {
            // x is in the cached interval
        } else {
            // Find new interval
            auto it = std::lower_bound(x_.begin(), x_.end(), x);
            accel.last_index_ = (it == x_.begin()) ? 0 : std::distance(x_.begin(), it) - 1;
            if (accel.last_index_ >= x_.size() - 1) {
                accel.last_index_ = x_.size() - 2;
            }
        }
        
        size_t i = accel.last_index_;
        double dx = x - x_[i];
        return y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }
    
private:
    /// \brief Compute natural spline coefficients.
    void computeCoefficients() {
        size_t n = x_.size();
        b_.resize(n);
        c_.resize(n);
        d_.resize(n);
        
        std::vector<double> h(n - 1);
        std::vector<double> alpha(n - 1);
        
        // Compute differences
        for (size_t i = 0; i < n - 1; ++i) {
            h[i] = x_[i + 1] - x_[i];
            alpha[i] = 3.0 * ((y_[i + 1] - y_[i]) / h[i] - (i > 0 ? (y_[i] - y_[i - 1]) / h[i - 1] : 0.0));
        }
        
        // Natural spline: second derivatives at endpoints are zero
        std::vector<double> l(n, 1.0);
        std::vector<double> mu(n, 0.0);
        std::vector<double> z(n, 0.0);
        
        // Forward elimination
        for (size_t i = 1; i < n - 1; ++i) {
            l[i] = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        
        // Back substitution
        c_[n - 1] = 0.0;  // Natural spline
        for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
            c_[i] = z[i] - mu[i] * c_[i + 1];
            b_[i] = (y_[i + 1] - y_[i]) / h[i] - h[i] * (c_[i + 1] + 2.0 * c_[i]) / 3.0;
            d_[i] = (c_[i + 1] - c_[i]) / (3.0 * h[i]);
        }
        
        // Set last coefficients
        b_[n - 1] = 0.0;
        d_[n - 1] = 0.0;
    }
    
    /// \brief Linear extrapolation to the left of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    double extrapolateLeft(double x) const {
        double dx = x - x_[0];
        return y_[0] + b_[0] * dx;
    }
    
    /// \brief Linear extrapolation to the right of the data range.
    /// \param x Input: x-coordinate.
    /// \return Output: extrapolated value.
    double extrapolateRight(double x) const {
        size_t n = x_.size();
        double dx = x - x_[n - 1];
        return y_[n - 1] + b_[n - 2] * dx;
    }
    
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> b_;  // Linear coefficients
    std::vector<double> c_;  // Quadratic coefficients
    std::vector<double> d_;  // Cubic coefficients
};

// Compatibility aliases for GSL-like interface
/// \brief GSL-compatible spline type alias.
using gsl_spline = CubicSpline;
/// \brief GSL-compatible accelerator type alias.
using gsl_interp_accel = CubicSpline::Accelerator;

/// \brief Allocate a spline instance (type/size ignored).
/// \param type Input: interpolation type (unused).
/// \param size Input: number of points (unused).
/// \return Output: spline pointer.
inline CubicSpline* gsl_spline_alloc(int /* type */, size_t /* size */) {
    return new CubicSpline();
}

/// \brief Allocate an interpolation accelerator.
/// \return Output: accelerator pointer.
inline CubicSpline::Accelerator* gsl_interp_accel_alloc() {
    return new CubicSpline::Accelerator();
}

/// \brief Initialize a spline with tabulated data.
/// \param spline Input/Output: spline to initialize.
/// \param x Input: strictly increasing x-coordinates.
/// \param y Input: y-values corresponding to \p x.
/// \param n Input: number of samples.
inline void gsl_spline_init(CubicSpline* spline, const double* x, const double* y, size_t n) {
    spline->init(x, y, n);
}

/// \brief Evaluate a spline at \p x using an accelerator.
/// \param spline Input: spline to evaluate.
/// \param x Input: x-coordinate.
/// \param accel Input/Output: accelerator cache.
/// \return Output: interpolated value.
inline double gsl_spline_eval(const CubicSpline* spline, double x, CubicSpline::Accelerator* accel) {
    return spline->eval(x, *accel);
}

/// \brief Free a spline instance.
/// \param spline Input: spline to release (can be null).
inline void gsl_spline_free(CubicSpline* spline) {
    delete spline;
}

/// \brief Free an accelerator instance.
/// \param accel Input: accelerator to release (can be null).
inline void gsl_interp_accel_free(CubicSpline::Accelerator* accel) {
    delete accel;
}

/// \brief Reset an accelerator to the initial state.
/// \param accel Input/Output: accelerator to reset.
inline void gsl_interp_accel_reset(CubicSpline::Accelerator* accel) {
    if (accel) {
        accel->last_index_ = 0;
    }
}

// GSL interpolation type constants (not used in our implementation, but kept for compatibility)
/// \brief GSL linear interpolation type identifier.
constexpr int gsl_interp_linear = 1;
constexpr int gsl_interp_cspline = 0;

#endif // OPAL_CUBIC_SPLINE_HH

