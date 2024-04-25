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
#include "PartBunch/PartBunch.hpp"

#include "Structure/FieldSolverCmd.h"

#include <boost/bimap.hpp>

#include <string>
#include <vector>

class Beam;
class OpalData;
class DataSink;
class Distribution;
class H5PartWrapper;
class Inform;
class Tracker;

class TrackRun : public Action {
public:
    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /// Make clone.
    virtual TrackRun* clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

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

    void initDataSink();

    void setupBoundaryGeometry();

    double setupDistribution(Beam* beam);

    Tracker* itsTracker_m;

    Distribution* dist_m;

    std::vector<Distribution*> distrs_m;

    FieldSolverCmd* fs_m;

    DataSink* ds_m;

    H5PartWrapper* phaseSpaceSink_m;

    OpalData* opal_m;

    /*

      this is the ippl bunch
    */

    using bunch_type = PartBunch_t;
    std::unique_ptr<bunch_type> bunch_m;

    bool isFollowupTrack_m;

    static const std::string defaultDistribution_m;

    RunMethod method_m;
    static const boost::bimap<RunMethod, std::string> stringMethod_s;

    /*
      \todo this is only for one bunch / beam
    */

    // macro mass / charge for simulation particles
    double macromass_m;
    double macrocharge_m;
};

inline Inform& operator<<(Inform& os, const TrackRun& b) {
    return b.print(os);
}

#endif  // OPAL_TrackRun_HH
