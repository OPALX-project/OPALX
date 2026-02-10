//
// Namespace Configure
//   The OPAL configurator.
//   This class must be modified to configure the commands to be contained
//   in an executable OPAL program. For each command an exemplar object
//   is constructed and linked to the main directory. This exemplar is then
//   available to the OPAL parser for cloning.
//   This class could be part of the class OpalData.  It is separated from
//   that class and opale into a special module in order to reduce
//   dependencies between modules.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "OpalConfigure/Configure.h"
#include "AbstractObjects/OpalData.h"

#include "Distribution/Distribution.h"

// Basic action commands.
#include "BasicActions/Quit.h"

// Commands introducing a special mode.
#include "Track/TrackCmd.h"

// Table-related commands.
#include "Structure/Beam.h"
#include "Structure/FieldSolverCmd.h"

// Value definitions commands.
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalDrift.h"
#include "Elements/OpalMarker.h"

// Structure-related commands.
#include "Lines/Line.h"

#include "changes.h"

// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace {

    void makeActions() {
        OpalData* opal = OpalData::getInstance();
        opal->create(new Quit());
        opal->create(new TrackCmd());
    }

    void makeDefinitions() {
        OpalData* opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->create(new RealConstant());
        opal->create(new StringConstant());

        opal->create(new Beam());
        opal->create(new FieldSolverCmd());
        opal->create(new Distribution());
    }

    void makeElements() {
        OpalData* opal = OpalData::getInstance();
        opal->create(new OpalDrift());
        opal->create(new OpalMarker());
        opal->create(new Line());
    }
};  // namespace

namespace Configure {
    void configure() {
        makeDefinitions();
        makeElements();
        makeActions();
        Versions::fillChanges();
    }
};  // namespace Configure
