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
#include "BasicActions/Call.h"
#include "BasicActions/DumpFields.h"
#include "BasicActions/DumpEMFields.h"
#include "BasicActions/Echo.h"
#include "BasicActions/Help.h"
#include "BasicActions/Option.h"
#include "BasicActions/Select.h"
#include "BasicActions/Stop.h"
#include "BasicActions/Quit.h"
#include "BasicActions/System.h"
#include "BasicActions/PSystem.h"
#include "BasicActions/Title.h"
#include "BasicActions/Value.h"

// Macro command.                                                                                                                                                                                                                             
#include "OpalParser/MacroCmd.h"


// Commands introducing a special mode.
#include "Track/TrackCmd.h"    

// Table-related commands.
#include "Structure/FieldSolver.h"
#include "Structure/Beam.h"
//#include "Tables/List.h"

// Value definitions commands.
#include "ValueDefinitions/BoolConstant.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalMarker.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalRingDefinition.h"
#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"
#include "Elements/OpalVerticalFFAMagnet.h"
#include "Elements/OpalProbe.h"

// Structure-related commands.
#include "Lines/Line.h"

#include "changes.h"

// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace {

    void makeActions() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new Call());
        opal->create(new DumpFields());
        opal->create(new DumpEMFields());
        opal->create(new Echo());
        opal->create(new Help());
        opal->create(new Option());
        opal->create(new Select());
        opal->create(new Stop());
        opal->create(new Quit());
        opal->create(new PSystem());
        opal->create(new System());
        opal->create(new Title());
        opal->create(new TrackCmd());
        opal->create(new Value());
    }

    void makeDefinitions() {
        OpalData *opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->create(new BoolConstant());
        opal->create(new RealConstant());
        opal->create(new RealVariable());
        opal->create(new RealVector());
        opal->create(new StringConstant());

        opal->create(new MacroCmd());

        opal->create(new Beam());
        opal->create(new FieldSolver());
        opal->create(new Distribution());        
    }

    void makeElements() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new OpalDrift());
        opal->create(new OpalMarker());
        opal->create(new OpalProbe());
        opal->create(new OpalRingDefinition());
        opal->create(new Line());
        opal->create(new OpalOffset::OpalLocalCartesianOffset());
        opal->create(new OpalVerticalFFAMagnet());
    }
};

namespace Configure {
    void configure() {
        makeDefinitions();
        makeElements();
        makeActions();
        Versions::fillChanges();
    }
};
