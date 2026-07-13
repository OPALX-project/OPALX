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
#include "Algorithms/StepSizeConfig.h"

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include "Elements/OpalBeamline.h"
#include "Steppers/BorisPusher.h"
#include "Structure/BoundingBox.h"
#include "Structure/ValueRange.h"

class OrbitThreader {
public:
    OrbitThreader(
            const PartData& ref, const Vector_t<double, 3>& r, const Vector_t<double, 3>& p,
            double s, double maxDiffZBunch, double t, double dT, StepSizeConfig stepSizes,
            OpalBeamline& bl, bool isDesignBeam);

    void execute();

    IndexMap::value_t query(IndexMap::key_t::first_type step, IndexMap::key_t::second_type length);

    IndexMap::key_t getRange(const IndexMap::value_t::value_type& element) const;

    BoundingBox getBoundingBox() const;

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

    void trackBack();
    void integrate(const IndexMap::value_t& activeSet, double maxDrift = 10.0);
    bool containsCavity(const IndexMap::value_t& activeSet);
    void autophaseCavities(
            const IndexMap::value_t& activeSet, const std::set<std::string>& visitedElements);
    double getMaxDesignEnergy(const IndexMap::value_t& elementSet) const;

    void setDesignEnergy(ElementList& allElements, const std::set<std::string>& visitedElements);
    void computeBoundingBox();
    void updateBoundingBoxWithCurrentPosition();
    double computeDriftLengthToBoundingBox(
            const std::set<std::shared_ptr<ElementBase>>& elements,
            const Vector_t<double, 3>& position, const Vector_t<double, 3>& direction) const;

    void checkElementLengths(const std::set<std::shared_ptr<ElementBase>>& elements);
};

inline IndexMap::value_t OrbitThreader::query(
        IndexMap::key_t::first_type pathLength, IndexMap::key_t::second_type length) {
    return imap_m.query(pathLength, length);
}

inline IndexMap::key_t OrbitThreader::getRange(const IndexMap::value_t::value_type& element) const {
    return imap_m.getRange(element);
}

inline BoundingBox OrbitThreader::getBoundingBox() const { return globalBoundingBox_m; }
#endif
