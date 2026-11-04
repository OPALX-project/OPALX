//
// Class ReferenceTracer
//   Traces the reference particle of the primary beam once through the lattice before tracking.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_REFERENCE_TRACER_H
#define OPALX_REFERENCE_TRACER_H

#include "Algorithms/ExternalFieldRayTracker.h"
#include "Algorithms/GeneralElementMapping.h"
#include "Algorithms/StepSizeConfig.h"

#include <cstddef>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>

class OpalBeamline;
class PartData;

/**
 * @brief Traces the reference particle of the primary beam once through the lattice.
 *
 * The pass takes Boris steps of the smallest DT. At each step it applies the elements that
 * OpalBeamline::getElements() selects at the step midpoint. It first steps backwards over the
 * bunch tail, then forwards. On the way it
 * - autophases each RFCAVITY and TRAVELINGWAVE when the reference first reaches it,
 * - checks that the DT of the current segment resolves each element it reaches
 *   (ElementBase::getRequiredNumberOfTimeSteps()),
 * - checks the element reaches against the reference
 * (GeneralElementMapping::checkReferencePoint()),
 * - writes <basename>_DesignPath.dat, which the FFT2D5 solver reads.
 *
 * A LINE pass ends when the reference is outside every element and either has left the lab box
 * of all elements and moves away from it, or is two bunch lengths past the last ZSTOP. It also
 * ends when the step budget of the tracking schedule is used up. A RING pass ends after one
 * circumference.
 */
class ReferenceTracer {
public:
    using RayState = ExternalFieldRayTracker::State;

    /**
     * @param reference Particle data of the primary beam.
     * @param beamline Placed lattice whose element grid is built.
     * @param stepSizes Tracking schedule: DT, ZSTOP and MAXSTEPS per segment.
     * @param dt Time step of the pass [s], normally the smallest DT of the schedule.
     * @param ringPeriod Ring circumference [m]; nonpositive for a LINE.
     */
    ReferenceTracer(
            const PartData& reference, OpalBeamline& beamline, const StepSizeConfig& stepSizes,
            double dt, double ringPeriod);

    /**
     * @brief Trace the reference particle.
     * @param start Lab position [m], momentum (beta gamma), time [s] and path length [m].
     * @param bunchTail Distance [m] the bunch extends behind the reference. The pass first steps
     *                  back over it, but not before path length 0.
     * @return The state where the pass ended.
     * @throws OpalException If DT does not resolve an element or the reference would be deleted.
     */
    RayState execute(RayState start, double bunchTail);

    /// Whether a ray from @p position along @p direction hits @p box, including from inside.
    static bool rayHitsBox(
            const BoundingBox& box, const Vector_t<double, 3>& position,
            const Vector_t<double, 3>& direction);

private:
    using ElementSet = std::set<std::shared_ptr<ElementBase>>;

    /// One Boris step with the elements at the step midpoint; returns false on a material hit.
    bool step(double dt, ElementSet& fields);
    void checkElementLengths(const ElementSet& fields);
    void autophaseCavities(const ElementSet& fields);
    void openDesignPath();
    void writeDesignPath(
            const ExternalFieldRayTracker::Step& step, const ElementSet& fields, double dt);
    void checkReach(const Vector_t<double, 3>& position, double pathLength);
    void printSupportOutsideReach() const;
    bool reachedEnd(const ElementSet& fields) const;

    const PartData& reference_m;
    OpalBeamline& beamline_m;
    /// Advanced along the pass to find the DT of the current segment.
    StepSizeConfig stepSizes_m;
    double dt_m;
    const double ringPeriod_m;
    ExternalFieldRayTracker rayTracker_m;

    RayState ray_m;
    /// Path length [m] at the start of the forward pass.
    double startPathLength_m{0.0};
    /// Path length [m] past which a LINE pass may end: last ZSTOP plus two bunch lengths.
    double endPathLength_m{0.0};
    /// Signed step counter; negative during the forward pass after the backward steps.
    long step_m{0};
    /// Forward steps allowed: the schedule's steps at the finest DT.
    unsigned long long maxForwardSteps_m{0};
    BoundingBox labBounds_m;

    std::set<std::string> visitedElements_m;
    ElementSet previousFields_m;
    std::map<const ElementBase*, GeneralElementMapping::SupportOutsideBox> supportOutsideReach_m;

    std::ofstream designPath_m;
    std::size_t designPathFrequency_m;
};

#endif  // OPALX_REFERENCE_TRACER_H
