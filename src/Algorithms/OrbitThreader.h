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

class OrbitThreader {
public:
    OrbitThreader(
            const PartData& ref, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p,
            double s, double maxDiffZBunch, double t, double dT, StepSizeConfig stepSizes,
            OpalBeamline& bl, bool isDesignBeam, double period = 0.0);

    void execute();

    IndexMap::value_t query(IndexMap::key_t::first_type step, IndexMap::key_t::second_type length);

    IndexMap::key_t getRange(const IndexMap::value_t::value_type& element) const;

    BoundingBox getBoundingBox() const;

    /**
     * @brief Combined map in reference-path order, available after execute() when enabled.
     *
     * If the ordered element/section maps are \f$M_1,\ldots,M_N\f$, this is
     * \f[
     * M_{\rm total}=M_N\cdots M_2M_1.
     * \f]
     * The product contains every unique constant-active-set segment, including field-free
     * intervals. For a periodic OrbitThreader this is the one-turn map.
     */
    const std::optional<matrix6x6_t>& getCombinedLinearTransferMap() const;

    /**
     * @brief Absolute determinant error \f$|\det(M)-1|\f$ for the combined map.
     *
     * This checks phase-space volume preservation. It is necessary, but for a \f$6\times6\f$
     * map not sufficient, for symplecticity.
     */
    const std::optional<double>& getCombinedDeterminantResidual() const;

    /**
     * @brief Maximum component of \f$M^TJM-J\f$ for the combined map.
     *
     * This is a canonical-form diagnostic. The reported slopes and mechanical momenta are not
     * globally canonical, so a nonzero value does not by itself demonstrate nonsymplectic
     * physical tracking.
     */
    const std::optional<double>& getCombinedSymplecticResidual() const;

private:
    /// position of reference particle in lab coordinates
    Vector_t<double, 3> r_m;
    /// momentum of reference particle
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

    BorisPusher integrator_m;
    const PartData& reference_m;

    std::ofstream logger_m;
    size_t loggingFrequency_m;

    BoundingBox globalBoundingBox_m;

    using ReferenceSample = LinearTransferMapBuilder::ReferenceSample;
    using RayState = ExternalFieldRayTracker::State;
    ExternalFieldRayTracker rayTracker_m;

    bool collectReferenceSamples_m{false};
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
