//
// Class OrbitThreader
//
// This class determines the design path by tracking the reference particle through
// the 3D lattice.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_ORBITTHREADER_H
#define OPAL_ORBITTHREADER_H

#include "Algorithms/IndexMap.h"
#include "Algorithms/LinearTransferMapBuilder.h"
#include "Algorithms/StepSizeConfig.h"

#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include "Elements/OpalBeamline.h"
#include "Steppers/BorisPusher.h"
#include "Structure/BoundingBox.h"
#include "Structure/ValueRange.h"

/**
 * @brief Trace the reference orbit and orchestrate optional element-attached linear maps.
 *
 * The first pass builds an IndexMap of encountered field supports. For the design
 * beam with Options::enableLinearTransferMaps, it also samples reference states and
 * transported frames and marks observed support overlaps. A second pass delegates
 * the numerical work to LinearTransferMapBuilder: nominal bodies define map owners,
 * but reference and private rays always see the full support-selected external field.
 *
 * This class attaches map copies to runtime owners and keeps the complete combined
 * product, including unowned intervals. Secondary species may build their own
 * support index but do not replace the shared design-beam maps or overlap flags.
 * The selected map integrator is shared by the map-enabled reference and shadow rays;
 * it does not change the production bunch pusher. This is not a closed-orbit finder.
 */
class OrbitThreader {
public:
    /**
     * @brief Thread a LINE interval or a RING return section.
     *
     * Positive @p period is the nominal design circumference in metres, never the sum of
     * field-support lengths. A ring ends at the launch-direction crossing of
     * \f$(\mathbf r-\mathbf r_0)\cdot\widehat{\mathbf p}_0=0\f$ after travelling at least
     * half the design circumference. Search is bounded by twice that distance (or four
     * nominal flight times). The measured return length sets the distance-index period;
     * the design circumference is unchanged. This simple-ring return search does not
     * solve for a fixed point or certify position/momentum closure.
     *
     * @param ref Particle rest energy and charge; must outlive the threader.
     * @param r Initial reference position in the lab frame [m].
     * @param p Initial lab mechanical momentum in beta-gamma units, p/(mc).
     * @param s Initial accumulated reference-path coordinate [m].
     * @param maxDiffZBunch Longitudinal bunch extent [m], used for pre-roll/end padding.
     * @param t Initial time [s].
     * @param dT Signed nominal reference time step [s].
     * @param stepSizes Production stepping limits; map/ring threading can continue
     * beyond the production step-count budget to cover its requested interval.
     * @param bl Placed runtime lattice; must outlive the threader.
     * @param isDesignBeam Whether this threader may update shared design diagnostics.
     * @param period Nominal circumference [m] for RING; nonpositive selects LINE behavior.
     */
    OrbitThreader(
            const PartData& ref, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p,
            double s, double maxDiffZBunch, double t, double dT, StepSizeConfig stepSizes,
            OpalBeamline& bl, bool isDesignBeam, double period = 0.0);

    /**
     * @brief Run reference threading and, when enabled for the design beam, generate maps.
     *
     * A map-enabled design pass first clears existing runtime maps and overlap flags.
     * After successful construction, each owner receives its segment copies with a
     * zero-based per-element attachment ordinal. The full product and diagnostics
     * remain on this threader and are printed on rank zero.
     * Disabled-map and secondary-species passes leave existing element maps untouched.
     * This operation is not transactional: a failure need not restore prior diagnostics.
     */
    void execute();

    /// Support candidates around path coordinate step with half-width length, both [m].
    /// For RING the lookup wraps with the measured return period; see IndexMap::query().
    IndexMap::value_t query(IndexMap::key_t::first_type step, IndexMap::key_t::second_type length);

    /// First contiguous support crossing of this runtime occurrence, not all ring visits.
    IndexMap::key_t getRange(const IndexMap::value_t::value_type& element) const;

    BoundingBox getBoundingBox() const;

    /**
     * @brief Combined map in reference-path order, available after execute() when enabled.
     *
     * If the ordered element/section maps are \f$M_1,\ldots,M_N\f$, this is
     * \f[
     * M_{\rm total}=M_N\cdots M_2M_1.
     * \f]
     * The product contains every unique nominal-owner segment, including unowned intervals.
     * For a periodic OrbitThreader this is a same-section return derivative; it becomes a
     * closed-orbit one-turn map only after reference position and momentum closure.
     * Absent before map generation or when no nondegenerate map interval was produced.
     * The returned reference belongs to this threader, not the elements; copy the value
     * if it must outlive the threader.
     */
    const std::optional<matrix6x6_t>& getCombinedLinearTransferMap() const;

    /**
     * @brief Absolute determinant error \f$|\det(M)-1|\f$ for the combined map.
     *
     * This checks volume preservation in the reported coordinates. Unit determinant is
     * necessary but not sufficient for canonical symplecticity; the chosen noncanonical
     * coordinates need not preserve it universally. Availability/lifetime match the map.
     */
    const std::optional<double>& getCombinedDeterminantResidual() const;

    /**
     * @brief Maximum component of \f$M^TJM-J\f$ for the combined map.
     *
     * This is a canonical-form diagnostic. The reported slopes and mechanical momenta are not
     * globally canonical, so a nonzero value does not by itself demonstrate nonsymplectic
     * physical tracking. Uses the J defined in LinearTransferMap; availability/lifetime
     * match getCombinedLinearTransferMap().
     */
    const std::optional<double>& getCombinedSymplecticResidual() const;

    /// Nominal circumference argument [m], supplied by Ring::getLength(); nonpositive for LINE.
    double getDesignCircumference() const { return period_m; }
    /// Measured travel distance [m], populated only after a RING return-plane crossing.
    /// This is not necessarily a closed-orbit length and is independent of map enablement.
    const std::optional<double>& getReferenceReturnLength() const { return referenceReturnLength_m; }
    /// Current displacement from the launch position [m] in lab coordinates.
    /// Interpret as a return-closure diagnostic only when getReferenceReturnLength() is present.
    Vector_t<double, 3> getReferenceReturnDisplacement() const { return r_m - ringOrigin_m.position; }

private:
    /// position of reference particle in lab coordinates
    Vector_t<double, 3> r_m;
    /// Lab mechanical momentum of the reference particle in beta-gamma units.
    Vector_t<double, 3> p_m;
    /// position of reference particle in path length
    double pathLength_m;
    /// distance to track back before tracking forward
    /// (length of bunch but not beyond cathode)
    double distTrackBack_m;
    /// the simulated time
    double time_m;
    /// the time step
    double dt_m;
    ValueRange<long> stepRange_m;
    long currentStep_m{0};

    /// final position in path length
    StepSizeConfig stepSizes_m;
    const double sStop_m;
    const double period_m;
    ValueRange<double> pathLengthRange_m;

    OpalBeamline& itsOpalBeamline_m;
    IndexMap imap_m;

    /// True for the single design beam (container 0's species): only this threader may
    /// autophase cavities, set per-element design energy, and write the geometry/design-path
    /// outputs. Secondary species build their own IndexMap but reuse that shared element
    /// state.
    bool isDesignBeam_m;

    unsigned int errorFlag_m;

    const PartData& reference_m;

    std::ofstream logger_m;
    size_t loggingFrequency_m;

    BoundingBox globalBoundingBox_m;

    using ReferenceSample = LinearTransferMapBuilder::ReferenceSample;
    using RayState = ExternalFieldRayTracker::State;
    Vector_t<double, 3> positionCorrection_m = Vector_t<double, 3>(0.0);
    double timeCorrection_m{0.0};
    double pathLengthCorrection_m{0.0};
    RayState currentRay() const;
    RayState ringOrigin_m;
    bool referencePass_m{false};
    std::optional<double> referenceReturnLength_m;
    void checkRingReturn(const RayState& before, double stepDt);
    void setCurrentRay(const RayState& ray);
    Vector_t<double, 3> nextMidpointPosition() const;
    /// Snapshot of map numerical OPTION settings; map-enabled reference and perturbed rays
    /// use the same integrator. Other reference threading retains BORIS.
    LinearTransferMapBuilder::Settings mapSettings_m;
    ExternalFieldRayTracker rayTracker_m;

    /// Collect only the forward reference pass, not pre-roll tracking or private rays.
    bool collectReferenceSamples_m{false};
    /// Requested map entrance [m], excluding any reference pre-roll.
    double transferMapStartPathLength_m{0.0};
    std::vector<ReferenceSample> referenceSamples_m;
    std::optional<matrix6x6_t> combinedLinearTransferMap_m;
    std::optional<double> combinedDeterminantResidual_m;
    std::optional<double> combinedSymplecticResidual_m;

    void trackBack();
    void integrate(const IndexMap::value_t& activeSet, double maxDrift = 10.0);
    bool containsCavity(const IndexMap::value_t& activeSet);
    void autophaseCavities(
            const IndexMap::value_t& activeSet, const std::set<std::string>& visitedElements);
    double getMaxDesignEnergy(const IndexMap::value_t& elementSet) const;

    void setDesignEnergy(ElementList& allElements, const std::set<std::string>& visitedElements);
    void computeBoundingBox();
    void updateBoundingBoxWithCurrentPosition();
    bool reachedThreadingEnd() const;
    double computeDriftLengthToBoundingBox(
            const std::set<std::shared_ptr<ElementBase>>& elements,
            const Vector_t<double, 3>& position, const Vector_t<double, 3>& direction) const;

    void checkElementLengths(const std::set<std::shared_ptr<ElementBase>>& elements);

    /// Capture the compensated ray and transported axes. At RING return, reuse launch
    /// axes while retaining the actual returned position/momentum; do not impose closure.
    void recordReferenceSample();
    void printCombinedLinearTransferMap() const;
};

inline IndexMap::value_t OrbitThreader::query(
        IndexMap::key_t::first_type pathLength, IndexMap::key_t::second_type length) {
    return imap_m.query(pathLength, length);
}

inline IndexMap::key_t OrbitThreader::getRange(const IndexMap::value_t::value_type& element) const {
    return imap_m.getRange(element);
}

inline BoundingBox OrbitThreader::getBoundingBox() const { return globalBoundingBox_m; }

inline const std::optional<matrix6x6_t>& OrbitThreader::getCombinedLinearTransferMap() const {
    return combinedLinearTransferMap_m;
}

inline const std::optional<double>& OrbitThreader::getCombinedDeterminantResidual() const {
    return combinedDeterminantResidual_m;
}

inline const std::optional<double>& OrbitThreader::getCombinedSymplecticResidual() const {
    return combinedSymplecticResidual_m;
}
#endif
