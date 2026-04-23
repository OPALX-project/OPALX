//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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

#include "Algorithms/SplineTimeDependence.h"
#include <stdexcept>
#include "Utility/Inform.h"

SplineTimeDependence::SplineTimeDependence(
        const size_t splineOrder, const std::vector<double>& times,
        const std::vector<double>& values) {
    setSpline(splineOrder, times, values);
}

SplineTimeDependence* SplineTimeDependence::clone() {
    auto* timeDep = new SplineTimeDependence();
    timeDep->setSpline(splineOrder_m, times_m, values_m);
    return timeDep;
}

Inform& SplineTimeDependence::print(Inform& os) const {
    if (!spline_m) {
        os << "Uninitialised SplineTimeDependence" << endl;
    } else {
        os << "SplineTimeDependence of order " << splineOrder_m << " with " << times_m.size()
           << " entries" << endl;
    }
    return os;
}

void SplineTimeDependence::setSpline(
        const size_t splineOrder, const std::vector<double>& times,
        const std::vector<double>& values) {
    if (times.size() != values.size()) {
        throw std::invalid_argument(
                "SplineTimeDependence::SplineTimeDependence: "
                "Times and values should be of equal length");
    }
    if (times.size() <= splineOrder) {
        throw std::invalid_argument(
                "SplineTimeDependence::SplineTimeDependence: "
                "Times and values should be of length > splineOrder");
    }
    if (splineOrder != LinearInterpolation and splineOrder != CubicInterpolation) {
        throw std::invalid_argument(
                "SplineTimeDependence::SplineTimeDependence: "
                "Only linear or cubic interpolation is supported");
    }
    for (size_t i = 0; i < times.size() - 1; ++i) {
        if (times[i] >= times[i + 1]) {
            throw std::invalid_argument(
                    "SplineTimeDependence::SplineTimeDependence: "
                    "Times should increase monotonically");
        }
    }
    splineOrder_m = splineOrder;
    times_m       = times;
    values_m      = values;
    splineAcc_m   = std::make_unique<AbstractSpline::Accelerator>();
    if (splineOrder_m == LinearInterpolation) {
        spline_m = std::make_unique<LinearSpline>(times_m, values_m);
    } else {
        spline_m = std::make_unique<CubicSpline>(times_m, values_m);
    }
}

double SplineTimeDependence::getValue(const double time) {
    double result{};
    if (time < times_m[0] or time > times_m.back()) {
        std::stringstream ss;
        ss << "SplineTimeDependence::getValue: time out of spline range: " << time;
        throw std::invalid_argument(ss.str());
    }
    result = spline_m->eval(time, *splineAcc_m);
    return result;
}

double SplineTimeDependence::getIntegral(const double time) {
    double result{};
    if (time < times_m[0] or time > times_m.back()) {
        std::stringstream ss;
        ss << "SplineTimeDependence::getValue: time out of spline range: " << time;
        throw std::invalid_argument(ss.str());
    }
    result = spline_m->evalIntegral(0, time, *splineAcc_m);
    return result;
}
