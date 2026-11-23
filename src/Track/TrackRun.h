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

#include "AbstractObjects/Action.h"

class TrackRun : public Action {
public:
    /// Exemplar constructor.
    TrackRun();

    virtual ~TrackRun();

    /**
     * @brief Run freeze(), validate() and buildElements(); with --check write the run
     *        description and stop, else build the tracking setup and track.
     * @throws OpalException with every input error of the RUN, before anything is allocated.
     */
    virtual void execute();

private:
    // Not implemented.
    TrackRun(const TrackRun&);
    void operator=(const TrackRun&);
};

#endif  // OPAL_TrackRun_HH
