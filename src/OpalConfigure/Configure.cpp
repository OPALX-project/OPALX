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
        opal->addKeyword<Call>();
        opal->addKeyword<DumpEMFields>();
        opal->addKeyword<Echo>();
        opal->addKeyword<Help>();
        opal->addKeyword<Option>();
        opal->addKeyword<Select>();
        opal->addKeyword<Stop>();
        opal->addKeyword<Quit>();
        opal->addKeyword<PSystem>();
        opal->addKeyword<System>();
        opal->addKeyword<Title>();
        opal->addKeyword<TrackCmd>();
        opal->addKeyword<CofCmd>();
        opal->addKeyword<Value>();
    }

    void makeDefinitions() {
        OpalData* opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->addKeyword<BoolConstant>();
        opal->addKeyword<RealConstant>();
        RealConstant::defineBuiltins();
        opal->addKeyword<RealVariable>();
        RealVariable::defineBuiltins();
        opal->addKeyword<RealVector>();
        opal->addKeyword<StringConstant>();
        StringConstant::defineBuiltins();

        opal->addKeyword<MacroCmd>();

        opal->addKeyword<Beam>();
        Beam::defineBuiltins();
        opal->addKeyword<EmissionSource>();
        opal->addKeyword<FieldSolverCmd>();
        opal->addKeyword<Distribution>();
        opal->addKeyword<EmissionSourceList>();
        opal->addKeyword<BinningCmd>();
    }

    void makeElements() {
        OpalData* opal = OpalData::getInstance();
        opal->addKeyword<OpalCavity>();
        opal->addKeyword<OpalTravelingWave>();
        opal->addKeyword<OpalConstantEFieldCavity>();
        opal->addKeyword<OpalConstantFocusing>();
        opal->addKeyword<OpalCyclotronSector>();
        opal->addKeyword<OpalTrimCoil>();
        opal->addKeyword<OpalBox>();
        opal->addKeyword<OpalCollimator>();
        opal->addKeyword<OpalDrift>();
        opal->addKeyword<OpalLaser>();
        opal->addKeyword<OpalMonitor>();
        opal->addKeyword<OpalMarker>();
        OpalMarker::defineBuiltins();
        opal->addKeyword<OpalProbe>();
        opal->addKeyword<OpalMultipole>();
        opal->addKeyword<OpalQuadrupole>();
        opal->addKeyword<OpalRBend>();
        opal->addKeyword<OpalSBend>();
        opal->addKeyword<OpalScalingFFAMagnet>();
        opal->addKeyword<OpalFieldmapElement>();
        opal->addKeyword<OpalSolenoid>();
        opal->addKeyword<Line>();
        opal->addKeyword<Ring>();
        opal->addKeyword<OpalVerticalFFAMagnet>();
        opal->addKeyword<OpalPolynomialTimeDependence>();
        opal->addKeyword<OpalSinusoidalTimeDependence>();
        opal->addKeyword<OpalSplineTimeDependence>();
        opal->addKeyword<OpalMultipoleT>();
        opal->addKeyword<OpalVariableRFCavity>();

        opal->addKeyword<OpalEnge>();
        opal->addKeyword<OpalTanh>();
        opal->addKeyword<OpalAsymmetricEnge>();
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
