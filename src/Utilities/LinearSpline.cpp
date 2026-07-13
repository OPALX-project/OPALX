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

#include "LinearSpline.h"
#include <stdexcept>

LinearSpline::LinearSpline(const std::vector<double>& x, const std::vector<double>& y) {
    LinearSpline::init(x, y);
}

void LinearSpline::init(const std::vector<double>& x, const std::vector<double>& y) {
    // Base class handles the x and y points
    AbstractSpline::init(x, y);
    // The gradients of each interval
    computeCoefficients();
    // And the integrals over each interval
    computeIntegrals();
}

double LinearSpline::eval(const double x, Accelerator& accel) const {
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
        result         = m_[i] * x + c_[i];
    }
    return result;
}

double LinearSpline::evalIntegral(double xa, double xb, Accelerator& accel) const {
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

void LinearSpline::computeCoefficients() {
    m_.resize(x_.size() - 1);
    c_.resize(x_.size() - 1);
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        m_[i] = (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
        c_[i] = y_[i] - m_[i] * x_[i];
    }
}

void LinearSpline::computeIntegrals() {
    integrals_.resize(x_.size() - 1);
    for (size_t i = 0; i < x_.size() - 1; ++i) {
        integrals_[i] = integral(i, x_[i + 1] - x_[i]);
    }
}

double LinearSpline::integral(const size_t i, const double dx) const {
    const auto x2 = x_[i] + dx;
    return m_[i] * (x2 * x2 - x_[i] * x_[i]) / 2.0 + c_[i] * (x2 - x_[i]);
}

double LinearSpline::extrapolateLeft(const double x) const { return m_.front() * x + c_.front(); }

double LinearSpline::extrapolateRight(const double x) const { return m_.back() * x + c_.back(); }
