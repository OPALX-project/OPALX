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

#ifndef ALGORITHMS_SPLINETIMEDEPENDENCE_H_
#define ALGORITHMS_SPLINETIMEDEPENDENCE_H_

#include <vector>

#include "Algorithms/AbstractTimeDependence.h"
#include "Utilities/CubicSpline.h"
#include "Utilities/LinearSpline.h"

class Inform;

/** @class SplineTimeDependence
 *
 *  Time dependence that follows a spline. Interpolation is supported at
 *  linear or cubic order..
 *
 *  Interpolation is done using gsl_spline routines (1st or 3rd order).
 *
 *  Nb: I tried to use PolynomialPatch, as it is arbitrary order interpolation,
 *  but it turns out that only works for regular grids.
 */
class SplineTimeDependence : public AbstractTimeDependence {
public:
    /** Constructor
     *
     *  @param splineOrder the order of the spline used to fit the time
     *         dependence; either 1 (linear interpolation) or 3 (cubic spline
     *         with quadratic smoothing)
     *  @param times the times of successive elements in the time dependence
     *  @param values the values at each time step.
     *
     *  It is an error if times and values are not of equal length, times and
     *  values length < splineOrder + 1, or times do not increase strictly
     *  monotonically.
     */
    SplineTimeDependence(
            size_t splineOrder, const std::vector<double>& times,
            const std::vector<double>& values);

    /** Default Constructor makes a dependence of length 2 with values 0*/
    SplineTimeDependence() = default;

    /** Destructor cleans up the GSL spline stuff */
    ~SplineTimeDependence() override = default;

    /** Return the value of the spline at a given time
     *
     *  @param time: time at which the spline is referenced. If time is off
     *         either end of the spline, the last few values are used to
     *         extrapolate past the end of the spline.
     *
     */
    double getValue(double time) override;

    /** Return the integral of the spline from 0 to the given time
     *
     *  @param time: time at which the spline is referenced.
     *
     */
    double getIntegral(double time) override;

    /** Inheritable copy constructor
     *
     *  @returns new SplineTimeDependence.
     */
    SplineTimeDependence* clone() override;

    /** Print summary information about the time dependence
     *
     *  @param os "Inform" stream to which the information is printed.
     */
    Inform& print(Inform& os) const;

    /** Set the spline, deleting any existing spline data
     *
     * @param splineOrder 1 for linear, 3 for cubic, all other values invalid
     * @param times a list of times in seconds
     * @param values a list of values corresponding to the times
     */
    void setSpline(
            size_t splineOrder, const std::vector<double>& times,
            const std::vector<double>& values);

    /* Getters for the test case use */
    [[nodiscard]] const std::vector<double>& getTimes() const { return times_m; }
    [[nodiscard]] const std::vector<double>& getValues() const { return values_m; }
    [[nodiscard]] size_t getSplineOrder() const { return splineOrder_m; }

    // Spline order constants
    static constexpr size_t LinearInterpolation = 1;
    static constexpr size_t CubicInterpolation  = 3;

private:
    std::unique_ptr<AbstractSpline> spline_m;
    std::unique_ptr<AbstractSpline::Accelerator> splineAcc_m;
    size_t splineOrder_m{LinearInterpolation};
    std::vector<double> times_m;
    std::vector<double> values_m;
};

inline Inform& operator<<(Inform& os, const SplineTimeDependence& p) { return p.print(os); }

#endif
