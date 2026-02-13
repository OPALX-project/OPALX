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

#include "AbstractSpline.h"
#include <algorithm>
#include <stdexcept>

void AbstractSpline::init(const std::vector<double>& x, const std::vector<double>& y) {
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
}

size_t AbstractSpline::findInterval(const double x, size_t& intervalCache) const {
    if (intervalCache < x_.size() - 1 && x >= x_[intervalCache] && x < x_[intervalCache + 1]) {
        // x is in the cached interval
    } else {
        // Find the new interval
        const auto it = std::ranges::lower_bound(x_, x);
        intervalCache = it == x_.begin() ? 0 : std::distance(x_.begin(), it) - 1;
    }
    return intervalCache;
}
