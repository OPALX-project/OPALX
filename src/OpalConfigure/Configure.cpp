// ------------------------------------------------------------------------
// $RCSfile: Configure.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Configure
//   Contains methods for configuring the OPAL-9 program.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 12:40:49 $
// $Author: opal $
//
// JMJ & JP adding Aperture and Split 18/4/2000
// ------------------------------------------------------------------------

#include "OpalConfigure/Configure.h"
#include "AbstractObjects/OpalData.h"

#include "Distribution/Distribution.h"

// Basic action commands.
#include "BasicActions/Call.h"
#include "BasicActions/Dump.h"
#include "BasicActions/Echo.h"
#include "BasicActions/Help.h"
#include "BasicActions/Option.h"
#include "BasicActions/Save.h"
#include "BasicActions/Select.h"
#include "BasicActions/Show.h"
#include "BasicActions/Stop.h"
#include "BasicActions/Quit.h"
#include "BasicActions/What.h"
#include "BasicActions/System.h"
#include "BasicActions/Title.h"
#include "BasicActions/Value.h"

// Macro command.
#include "OpalParser/MacroCmd.h"

// Physics action commands.
#include "PhysicsActions/Dynamic.h"
#include "PhysicsActions/MakeSequence.h"
#include "PhysicsActions/SetIntegrator.h"
#include "PhysicsActions/Static.h"

// Commands introducing a special mode.
#include "Editor/EditCmd.h"
#include "Errors/ErrorCmd.h"
#include "Match/MatchCmd.h"
#include "Track/TrackCmd.h"

// Table-related commands.
#include "Structure/Beam.h"
#include "Structure/FieldSolver.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/OpalWake.h"
#include "Structure/SurfacePhysics.h"
#include "Utilities/OpalFilter.h"
#include "Tables/AttList.h"
#include "Tables/Eigen.h"
#include "Tables/Envelope.h"
#include "Tables/Insertion.h"
#include "Tables/List.h"
#include "Tables/MatrixCmd.h"
#include "Tables/Micado.h"
#include "Tables/Period.h"
#include "Tables/Survey.h"
#include "Tables/ThreadAll.h"
#include "Tables/ThreadBpm.h"
#include "Tables/Twiss3.h"
#include "Aperture/Aperture.h"
#include "Aperture/Split.h"


// Value definitions commands.
#include "ValueDefinitions/BoolConstant.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalBeamBeam.h"
#include "Elements/OpalBeamBeam3D.h"
#include "Elements/OpalCavity.h"
#include "Elements/OpalCCollimator.h"
#include "Elements/OpalCyclotron.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalECollimator.h"
#include "Elements/OpalDegrader.h"
#include "Elements/OpalHKicker.h"
#include "Elements/OpalHMonitor.h"
#include "Elements/OpalInstrument.h"
#include "Elements/OpalKicker.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalMultipole.h"
#include "Elements/OpalOctupole.h"
#include "Elements/OpalPepperPot.h"
#include "Elements/OpalPatch.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalQuadrupole.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalRCollimator.h"
#include "Elements/OpalSBend.h"
#include "Elements/OpalSBend3D.h"
#include "Elements/OpalSeparator.h"
#include "Elements/OpalSeptum.h"
#include "Elements/OpalSextupole.h"
#include "Elements/OpalSlit.h"
#include "Elements/OpalSolenoid.h"
#include "Elements/OpalSRot.h"
#include "Elements/OpalTravelingWave.h"
#include "Elements/OpalVKicker.h"
#include "Elements/OpalVMonitor.h"
#include "Elements/OpalWire.h"
#include "Elements/OpalYRot.h"
#include "Elements/OpalParallelPlate.h"
#include "Elements/OpalCyclotronValley.h"
#include "Elements/OpalStripper.h"
#include "Elements/OpalRingDefinition.h"

// Structure-related commands.
#include "Lines/Line.h"
#include "Lines/Sequence.h"


// Namespace Configure
// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace Configure {

    void makeActions() {
        OpalData *OPAL = OpalData::getInstance();
        OPAL->create(new Call());
        OPAL->create(new Dump());
        OPAL->create(new Echo());
        OPAL->create(new Dynamic());
        OPAL->create(new Eigen());
        OPAL->create(new Envelope());
        OPAL->create(new Help());
        OPAL->create(new EditCmd());
        OPAL->create(new ErrorCmd());
        OPAL->create(new List());
        OPAL->create(new MakeSequence());
        OPAL->create(new MatchCmd());
        OPAL->create(new MatrixCmd());
        OPAL->create(new Micado());
        OPAL->create(new Option());
        OPAL->create(new Save());
        OPAL->create(new Select());
        OPAL->create(new Show());
        OPAL->create(new SetIntegrator());
        OPAL->create(new Static());
        OPAL->create(new Stop());
        OPAL->create(new Quit());
        OPAL->create(new System());
        OPAL->create(new ThreadAll());
        OPAL->create(new ThreadBpm());
        OPAL->create(new Title());
        OPAL->create(new TrackCmd());
        OPAL->create(new Twiss3());
        OPAL->create(new Aperture());
        OPAL->create(new MSplit());
        OPAL->create(new Value());
        OPAL->create(new What());
    }


    void makeDefinitions() {
        OpalData *OPAL = OpalData::getInstance();
        // Must create the value definitions first.
        OPAL->create(new BoolConstant());
        OPAL->create(new RealConstant());
        OPAL->create(new RealVariable());
        OPAL->create(new RealVector());
        OPAL->create(new StringConstant());

        OPAL->create(new AttList());
        OPAL->create(new Beam());
        OPAL->create(new FieldSolver());
        OPAL->create(new BoundaryGeometry());
        OPAL->create(new OpalWake());
        OPAL->create(new SurfacePhysics());

        OPAL->create(new OpalFilter());

        OPAL->create(new Distribution());

        OPAL->create(new MacroCmd());
        OPAL->create(new Period());
        OPAL->create(new Insertion());
        OPAL->create(new Survey());
    }


    void makeElements() {
        OpalData *OPAL = OpalData::getInstance();
        OPAL->create(new OpalBeamBeam());
        OPAL->create(new OpalBeamBeam3D());
        OPAL->create(new OpalCavity());
        OPAL->create(new OpalCCollimator());
        OPAL->create(new OpalCyclotron());
        OPAL->create(new OpalDrift());
        OPAL->create(new OpalECollimator());
        OPAL->create(new OpalDegrader());
        OPAL->create(new OpalHKicker());
        OPAL->create(new OpalHMonitor());
        OPAL->create(new OpalInstrument());
        OPAL->create(new OpalKicker());
        OPAL->create(new OpalMarker());
        OPAL->create(new OpalMonitor());
        OPAL->create(new OpalMultipole());
        OPAL->create(new OpalOctupole());
        OPAL->create(new OpalPatch());
        OPAL->create(new OpalProbe());
        OPAL->create(new OpalPepperPot());
        OPAL->create(new OpalQuadrupole());
        OPAL->create(new OpalRBend());
        OPAL->create(new OpalRCollimator());
        OPAL->create(new OpalSBend());
        OPAL->create(new OpalSBend3D());
        OPAL->create(new OpalSeparator());
        OPAL->create(new OpalSeptum());
        OPAL->create(new OpalSextupole());
        OPAL->create(new OpalSlit());
        OPAL->create(new OpalSolenoid());
        OPAL->create(new OpalSRot());
        OPAL->create(new OpalTravelingWave());
        OPAL->create(new OpalVKicker());
        OPAL->create(new OpalVMonitor());
        OPAL->create(new OpalWire());
        OPAL->create(new OpalYRot());
        OPAL->create(new OpalParallelPlate());
        OPAL->create(new OpalCyclotronValley());
        OPAL->create(new OpalStripper());
        OPAL->create(new Line());
        OPAL->create(new Sequence());
        OPAL->create(new OpalRingDefinition());
    }


    void configure() {
        makeDefinitions();
        makeElements();
        makeActions();
    }
};
