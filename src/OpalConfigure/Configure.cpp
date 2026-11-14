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
#include "BasicActions/DumpEMFields.h"
#include "BasicActions/Echo.h"
#include "BasicActions/Help.h"
#include "BasicActions/Option.h"
#include "BasicActions/PSystem.h"
#include "BasicActions/Quit.h"
#include "BasicActions/Select.h"
#include "BasicActions/Stop.h"
#include "BasicActions/System.h"
#include "BasicActions/Title.h"
#include "BasicActions/Value.h"

// Macro command.
#include "OpalParser/MacroCmd.h"

// Commands introducing a special mode.
#include "Track/CofCmd.h"
#include "Track/TrackCmd.h"

// Table-related commands.
#include "Structure/Beam.h"
#include "Structure/BinningCmd.h"
#include "Structure/EmissionSource.h"
#include "Structure/FieldSolverCmd.h"
// #include "Tables/List.h"

// Value definitions commands.
#include "ValueDefinitions/BoolConstant.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalBox.h"
#include "Elements/OpalCavity.h"
#include "Elements/OpalCollimator.h"
#include "Elements/OpalConstantEFieldCavity.h"
#include "Elements/OpalConstantFocusing.h"
#include "Elements/OpalCyclotronSector.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalFieldmapElement.h"
#include "Elements/OpalLaser.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalMultipole.h"
#include "Elements/OpalMultipoleT.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalQuadrupole.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalSBend.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Elements/OpalSinusoidalTimeDependence.h"
#include "Elements/OpalSolenoid.h"
#include "Elements/OpalSplineTimeDependence.h"
#include "Elements/OpalTravelingWave.h"
#include "Elements/OpalVariableRFCavity.h"
#include "Elements/OpalVerticalFFAMagnet.h"
#include "Lines/Ring.h"

// EndFieldModels
#include "Elements/OpalAsymmetricEnge.h"
#include "Elements/OpalEnge.h"
#include "Elements/OpalTanh.h"

// Structure-related commands.
#include "Lines/EmissionSourceList.h"
#include "Lines/Line.h"

#include "changes.h"

// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace {

    void makeActions() {
        OpalData* opal = OpalData::getInstance();
        opal->addKeyword(new Call());
        opal->addKeyword(new DumpEMFields());
        opal->addKeyword(new Echo());
        opal->addKeyword(new Help());
        opal->addKeyword(new Option());
        opal->addKeyword(new Select());
        opal->addKeyword(new Stop());
        opal->addKeyword(new Quit());
        opal->addKeyword(new PSystem());
        opal->addKeyword(new System());
        opal->addKeyword(new Title());
        opal->addKeyword(new TrackCmd());
        opal->addKeyword(new CofCmd());
        opal->addKeyword(new Value());
    }

    void makeDefinitions() {
        OpalData* opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->addKeyword(new BoolConstant());
        opal->addKeyword(new RealConstant());
        RealConstant::defineBuiltins();
        opal->addKeyword(new RealVariable());
        RealVariable::defineBuiltins();
        opal->addKeyword(new RealVector());
        opal->addKeyword(new StringConstant());
        StringConstant::defineBuiltins();

        opal->addKeyword(new MacroCmd());

        opal->addKeyword(new Beam());
        Beam::defineBuiltins();
        opal->addKeyword(new EmissionSource());
        opal->addKeyword(new FieldSolverCmd());
        opal->addKeyword(new Distribution());
        opal->addKeyword(new EmissionSourceList());
        opal->addKeyword(new BinningCmd());
    }

    void makeElements() {
        OpalData* opal = OpalData::getInstance();
        opal->addKeyword(new OpalCavity());
        opal->addKeyword(new OpalTravelingWave());
        opal->addKeyword(new OpalConstantEFieldCavity());
        opal->addKeyword(new OpalConstantFocusing());
        opal->addKeyword(new OpalCyclotronSector());
        opal->addKeyword(new OpalTrimCoil());
        opal->addKeyword(new OpalBox());
        opal->addKeyword(new OpalCollimator());
        opal->addKeyword(new OpalDrift());
        opal->addKeyword(new OpalLaser());
        opal->addKeyword(new OpalMonitor());
        opal->addKeyword(new OpalMarker());
        OpalMarker::defineBuiltins();
        opal->addKeyword(new OpalProbe());
        opal->addKeyword(new OpalMultipole());
        opal->addKeyword(new OpalQuadrupole());
        opal->addKeyword(new OpalRBend());
        opal->addKeyword(new OpalSBend());
        opal->addKeyword(new OpalScalingFFAMagnet());
        opal->addKeyword(new OpalFieldmapElement());
        opal->addKeyword(new OpalSolenoid());
        opal->addKeyword(new Line());
        opal->addKeyword(new Ring());
        opal->addKeyword(new OpalVerticalFFAMagnet());
        opal->addKeyword(new OpalPolynomialTimeDependence());
        opal->addKeyword(new OpalSinusoidalTimeDependence());
        opal->addKeyword(new OpalSplineTimeDependence());
        opal->addKeyword(new OpalMultipoleT());
        opal->addKeyword(new OpalVariableRFCavity());

        opal->addKeyword(new OpalEnge());
        opal->addKeyword(new OpalTanh());
        opal->addKeyword(new OpalAsymmetricEnge());
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
