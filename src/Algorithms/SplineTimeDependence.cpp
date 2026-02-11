/*
 *  Copyright (c) 2018, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "Algorithms/SplineTimeDependence.h"
#include "Utilities/GSLCubicSpline.h"
#include "Utilities/GeneralClassicException.h"
#include "Utility/Inform.h"

SplineTimeDependence::SplineTimeDependence(
    const size_t splineOrder, const std::vector<double>& times, const std::vector<double>& values) {
    setSpline(splineOrder, times, values);
}

SplineTimeDependence::SplineTimeDependence(const SplineTimeDependence& rhs) {
    setSpline(rhs.splineOrder_m, rhs.times_m, rhs.values_m);
}

SplineTimeDependence* SplineTimeDependence::clone() {
    auto* timeDep = new SplineTimeDependence();
    timeDep->setSpline(splineOrder_m, times_m, values_m);
    return timeDep;
}

Inform& SplineTimeDependence::print(Inform& os) const {
    if (!cubicSpline_m && !linearSpline_m) {
        os << "Uninitialised SplineTimeDependence" << endl;
    } else {
        os << "SplineTimeDependence of order " << splineOrder_m << " with " << times_m.size()
           << " entries" << endl;
    }
    return os;
}

void SplineTimeDependence::setSpline(
    const size_t splineOrder, const std::vector<double>& times, const std::vector<double>& values) {
    if (times.size() != values.size()) {
        throw GeneralClassicException(
            "SplineTimeDependence::SplineTimeDependence",
            "Times and values should be of equal length");
    }
    if (times.size() <= splineOrder) {
        throw GeneralClassicException(
            "SplineTimeDependence::SplineTimeDependence",
            "Times and values should be of length > splineOrder");
    }
    if (splineOrder != LinearInterpolation and splineOrder != CubicInterpolation) {
        throw GeneralClassicException(
            "SplineTimeDependence::SplineTimeDependence",
            "Only linear or cubic interpolation is supported");
    }
    for (size_t i = 0; i < times.size() - 1; ++i) {
        if (times[i] >= times[i + 1]) {
            throw GeneralClassicException(
                "SplineTimeDependence::SplineTimeDependence",
                "Times should increase monotonically");
        }
    }
    splineOrder_m = splineOrder;
    times_m  = times;
    values_m = values;
    if (splineOrder_m == LinearInterpolation) {
        linearSpline_m = std::make_unique<LinearSpline>(times_m, values_m);
        linearAcc_m = std::make_unique<LinearSpline::Accelerator>();
        cubicSpline_m.reset();
        cubicAcc_m.reset();
    } else {
        linearSpline_m.reset();
        linearAcc_m.reset();
        cubicSpline_m = std::make_unique<CubicSpline>(times_m.data(), values_m.data(),
                times_m.size());
        cubicAcc_m = std::make_unique<CubicSpline::Accelerator>();
    }
}

double SplineTimeDependence::getValue(const double time) {
    double result{};
    if (time < times_m[0] or time > times_m.back()) {
        std::stringstream ss;
        ss << "time out of spline range: " << time;
        throw GeneralClassicException("SplineTimeDependence::getValue", ss.str());
    }
    if (splineOrder_m == LinearInterpolation) {
        result = linearSpline_m->eval(time, *linearAcc_m);
    } else {
        result = cubicSpline_m->eval(time, *cubicAcc_m);
    }
    return result;
}

double SplineTimeDependence::getIntegral(const double time) {
    double result{};
    if (time < times_m[0] or time > times_m.back()) {
        throw GeneralClassicException("SplineTimeDependence::getValue", "time out of spline range");
    }
    if (splineOrder_m == LinearInterpolation) {
        result = linearSpline_m->evalIntegral(0, time, *linearAcc_m);
    } else {
        result = cubicSpline_m->evalIntegral(0, time, *cubicAcc_m);
    }
    return result;
}
