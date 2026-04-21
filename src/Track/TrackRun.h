//
// Class TrackRun
//   The RUN command.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_TrackRun_HH
#define OPAL_TrackRun_HH

#include "OPALTypes.h"

#include "AbstractObjects/Action.h"
#include "PartBunch/PartBunch.h"
#include "Distribution/SamplingBase.hpp"

#include "Structure/FieldSolverCmd.h"

#include "Utilities/BiMap.h"

#include <string>
#include <vector>

class Beam;
class OpalData;
class DataSink;
class Distribution;
class EmissionSource;
class H5PartWrapper;
class Inform;
class Tracker;

class TrackRun : public Action {
    using emittingSamplers_t = std::vector<std::shared_ptr<SamplingBase>>;
public:
    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /// Make clone.
    virtual TrackRun* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

    // Bring base class print into scope to avoid hiding warning
    using Object::print;
    
    Inform& print(Inform& os) const;

private:
    enum class RunMethod : unsigned short { NONE, PARALLEL };

    // Not implemented.
    TrackRun(const TrackRun&);
    void operator=(const TrackRun&);

    // Clone constructor.
    TrackRun(const std::string& name, TrackRun* parent);

    void setRunMethod();
    std::string getRunMethodName() const;

    void initDataSink(size_t numParticleContainers);

    void setupBoundaryGeometry();

    /// Build samplers for all emission sources, perform initial sampling for t0 == 0
    /// sources, and populate emittingSamplers_m for time-dependent or delayed sources.
    /// Applied to particle container [index]
    void setupDistributionsAndSamplers(
        const std::vector<EmissionSource*>& sources, 
        Beam* beam, 
        emittingSamplers_t& emittingSamplers,
        size_t index=0
    );

    /// Compute total number of macroparticles for the bunch from BEAM::NALLOC and
    /// optional per-distribution NPARTDIST values on the emission sources.
    size_t computeTotalAllocationForBunch(
        Beam* beam,
        const std::vector<EmissionSource*>& sources) const;

    Tracker* itsTracker_m;

    /// Distributions referenced by all emission sources (non-owning raw pointers).
    std::vector<Distribution*> distrs_m;

    /// Samplers for time-dependent (emitting) sources; tracker calls emitParticles(t, dt) on each.
    //std::vector<std::shared_ptr<SamplingBase>> emittingSamplers_m;

    std::shared_ptr<FieldSolverCmd> fs_m;

    std::shared_ptr<DataSink> ds_m;

    std::vector<H5PartWrapper*> phaseSpaceSinks_m;

    OpalData* opal_m;

    /*

      this is the ippl bunch
    */

    using bunch_type = PartBunch_t;
    std::shared_ptr<bunch_type> bunch_m;

    bool isFollowupTrack_m;

    RunMethod method_m;
    static const BiMap<RunMethod, std::string> stringMethod_s;
    
};

inline Inform& operator<<(Inform& os, const TrackRun& b) {
    return b.print(os);
}

#endif  // OPAL_TrackRun_HH
