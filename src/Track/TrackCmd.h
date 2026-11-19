//
// Class TrackCmd
//   The class for the OPAL TRACK command.
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
#ifndef OPAL_TrackCmd_HH
#define OPAL_TrackCmd_HH

#include "AbstractObjects/Action.h"

#include <vector>

namespace opalx::run {
    struct TrackDescription;
}

/**
 * @brief The TRACK command. It stores its settings until ENDTRACK; each RUN in between reads
 *        them with openTrack().
 */
class TrackCmd : public Action {
public:
    TrackCmd();

    virtual ~TrackCmd();

    /// Store the TRACK settings. Throws if a TRACK is already open.
    virtual void execute();

    /// The settings of the TRACK that has no ENDTRACK yet; nullptr if there is none.
    static const opalx::run::TrackDescription* openTrack();

    /// Forget the open TRACK (ENDTRACK). Throws if there is none.
    static void closeTrack();

    /// Return the timestep in seconds
    std::vector<double> getDT() const;

    /// Return the elapsed time (sec) of the bunch
    double getT0() const;

    /// Return the maximum timsteps we integrate the system
    std::vector<unsigned long long> getMaxSteps() const;

    /// location at which the simulation starts
    double getZStart() const;

    /// location at which the simulation stops
    std::vector<double> getZStop() const;

private:
    // Not implemented.
    TrackCmd(const TrackCmd&);
    void operator=(const TrackCmd&);
};

#endif  // OPAL_TrackCmd_HH
