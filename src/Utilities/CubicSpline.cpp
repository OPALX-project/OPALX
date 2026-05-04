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
#include "CubicSpline.h"
#include <stdexcept>

CubicSpline::CubicSpline(const std::vector<double>& x, const std::vector<double>& y) {
    CubicSpline::init(x, y);
}

void CubicSpline::init(const std::vector<double>& x, const std::vector<double>& y) {
    // Base class handles the points
    AbstractSpline::init(x, y);
    // Natural cubic spline: second derivatives at endpoints are zero
    computeCoefficients();
    // And the integrals over each interval
    computeIntegrals();
}

double CubicSpline::eval(const double x, Accelerator& accel) const {
    double result{};
    if (x_.empty()) {
        throw std::runtime_error("CubicSpline: not initialized");
    }
    if (x <= x_.front()) {
        result = extrapolateLeft(x);
    } else if (x >= x_.back()) {
        result = extrapolateRight(x);
    } else {
        const size_t i  = findInterval(x, accel.last_index_);
        const double dx = x - x_[i];
        result          = y_[i] + b_[i] * dx + c_[i] * dx * dx + d_[i] * dx * dx * dx;
    }
    return result;
}

double CubicSpline::evalIntegral(double xa, double xb, Accelerator& accel) const {
    double result = 0.0;
    if (x_.empty()) {
        throw std::runtime_error("CubicSpline: not initialized");
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

void CubicSpline::computeCoefficients() {
    // Size the coefficients
    const size_t n = x_.size();
    b_.resize(n);
    c_.resize(n);
    d_.resize(n);
    // Compute differences
    std::vector<double> h(n - 1);
    std::vector<double> alpha(n - 1);
    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = x_[i + 1] - x_[i];
        alpha[i] =
                3.0 * ((y_[i + 1] - y_[i]) / h[i] - (i > 0 ? (y_[i] - y_[i - 1]) / h[i - 1] : 0.0));
    }
    // Natural spline: second derivatives at endpoints are zero
    std::vector l(n, 1.0);
    std::vector mu(n, 0.0);
    std::vector z(n, 0.0);
    // Forward elimination
    for (size_t i = 1; i < n - 1; ++i) {
        l[i]  = 2.0 * (x_[i + 1] - x_[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
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

void CubicSpline::computeIntegrals() {
    integrals_.resize(x_.size() - 1);
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        integrals_[i] = integral(i, x_[i + 1] - x_[i]);
    }
}

double CubicSpline::integral(const size_t i, const double dx) const {
    return y_[i] * dx + b_[i] * dx * dx / 2.0 + c_[i] * dx * dx * dx / 3.0
           + d_[i] * dx * dx * dx * dx / 4.0;
}

double CubicSpline::extrapolateLeft(double x) const {
    double dx = x - x_[0];
    return y_[0] + b_[0] * dx;
}

double CubicSpline::extrapolateRight(double x) const {
    size_t n  = x_.size();
    double dx = x - x_[n - 1];
    return y_[n - 1] + b_[n - 2] * dx;
}
