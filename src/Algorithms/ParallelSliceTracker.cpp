// ------------------------------------------------------------------------
// $RCSfile: ParallelSliceTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelSliceTracker
//   The visitor class for tracking particles with time as independent
//   variable.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include <cfloat>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Algorithms/ParallelSliceTracker.h"


#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"
#include "Lines/Sequence.h"

#include "Fields/BMultipoleField.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsMath.h"
#include "FixedAlgebra/FVps.h"

#include "Utilities/NumToStr.h"
#include "Distribution/Distribution.h"
#include "Utilities/Timer.h"

#define PSdim 6
class PartData;
using namespace OPALTimer;
using Physics::c;

extern Inform *gmsg;
extern Inform *gmsg2all;

// Class ParallelSliceTracker
// ------------------------------------------------------------------------

ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        const PartData &reference,
        bool revBeam,
        bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    itsOpalBeamline_m() {
}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        EnvelopeBunch &bunch,
        DataSink &ds,
        const PartData &reference,
        bool revBeam,
        bool revTrack,
        int maxSTEPS,
        double zstop):
    Tracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    maxSteps_m(maxSTEPS),
    zstop_m(zstop) {
    //itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    scaleFactor_m = itsBunch->getdT() * c;

    timeIntegrationTimer1_m  = IpplTimings::getTimer("Time integration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("Time integration2");
    timeFieldEvaluation_m  = IpplTimings::getTimer("Field evaluation");

    BinRepartTimer_m   = IpplTimings::getTimer("Time of Binary repart.");
    WakeFieldTimer_m   = IpplTimings::getTimer("Time of Wake Field calc.");

}


ParallelSliceTracker::~ParallelSliceTracker() {
}

void ParallelSliceTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelSliceTracker::applyExitFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}

void ParallelSliceTracker::execute() {

    Inform msg("ParallelSliceTracker");
    double tmp;
    double recpgamma, gamma;
    double t = itsBunch->getT();
    double dt = itsBunch->getdT();
    double tEmission = itsBunch->getTEmission();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);
    long long step = 0.0;

    if(OPAL.inRestartRun()) {
        int prevDumpFreq = OPAL.getRestartDumpFreq();
        step = OPAL.getRestartStep() * prevDumpFreq + 1;
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    else if(OPAL.hasBunchAllocated() && Options::scan) {
        step = 1;
        if(!itsBunch->doEmission())
            writePhaseSpace(step - 1, 0.0); // write initial phase space
        itsBunch->setT(0.0);
    } else {
        step = OPAL.getLastStep() + 1;
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    *gmsg << "Track start at: t= " << itsBunch->getT() << " zstop@ " << zstop_m << " [m]" << endl;

    Vector_t um, a, s;
    Vector_t externalE, externalB, KR, KT;
    BorisPusher pusher(itsReference);
    Vector_t rmin, rmax;

    bool global_EOL;
    bool bends;               // flag which indicates wheter any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch

    bool hasWake = false;     // flag which indicates wheter any particle is within the influence of a wake field

    //unsigned long bends = 0;            // flag which indicates whether any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch

    //unsigned long hasWake = 0;          // flag which indicates whether any particle is within the influence of a wake field

    int wfSection = -1;
    WakeFunction *wf;

    size_t totalParticles_i = itsBunch->getTotalNum();

    msg << "executing ParallelSliceTracker, initial DT " << itsBunch->getdT()
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM() * 1e-6 << " MeV, its charge: " << itsReference.getQ() << endl;

    itsBeamline_m.accept(*this);  // fill list allElements
    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(msg);

    for(int i = 0; i < itsBunch->getLocalNum(); i++) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));
        itsOpalBeamline_m.getSectionIndexAt(pos, l);
    }

    //itsBunch->calcBeamParameters();

    msg << "ParallelEnvelopeTracker starting tracking..." << endl;

    for(step; step < maxSteps_m; ++step) {

        //check if any particle hasn't reached the end of the field from the last element
        global_EOL = true;
        bends = false;
        hasWake = false;
        wfSection = -1;
        double margin = 1e-7;

        itsOpalBeamline_m.resetStatus();
        t = itsBunch->getT();
        itsOpalBeamline_m.switchElements(itsBunch->zTail() - margin, itsBunch->zHead() + margin);

        //external field for all local slices
        IpplTimings::startTimer(timeFieldEvaluation_m);
        for(int i = 0; i < itsBunch->getLocalNum(); i++) {

            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            KR        = Vector_t(0.0);
            KT        = Vector_t(0.0);

            Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));

            long &ls = itsBunch->LastSection[i];
            itsOpalBeamline_m.getSectionIndexAt(pos, ls);

            if(ls != itsBunch->LastSection[i])
                itsBunch->LastSection[i] = ls;

            //itsBunch->actT();
            //t = itsBunch->getT();

            const unsigned int &rtv = itsOpalBeamline_m.getFieldAt(i, pos, ls, t + itsBunch->dt[i] / 2.0, externalE, externalB);
            global_EOL = global_EOL && (rtv & BEAMLINE_EOL);

            //if((rtv & BEAMLINE_WAKE) && !hasWake) {
            //    wfSection = ls;
            //}

            //hasWake = hasWake || (rtv & BEAMLINE_WAKE);
            //bends = bends || (rtv & BEAMLINE_BEND);

            // Calculate factors for the envelope equation
            itsOpalBeamline_m.getKFactors(i, pos, ls, t + itsBunch->dt[i] / 2.0, KR, KT);

            // pass K-values and E/B fields to EnvelopeBunch
            itsBunch->setKR(KR, i);
            itsBunch->setKT(KT, i);
            itsBunch->setBF(externalB, i);
            itsBunch->setEF(externalE, i);
        }
        IpplTimings::stopTimer(timeFieldEvaluation_m);

        //TEST the non-len reduce: reduce(&global_EOL, &global_EOL, OpBitwiseOrAssign());
        reduce(&global_EOL, &global_EOL + 1, &global_EOL, OpBitwiseAndAssign());
        if(global_EOL)
            break;

        //FIXME: synchronizeSliceInformation() (currently done in updateFields())

        // make sure I-profile and space charge is updated
        itsBunch->updateFields();

        // do timestep for all slices
        // calls the function BetBunch->Run() and solves the envelope equation
        IpplTimings::startTimer(timeIntegrationTimer1_m);
        itsBunch->run(itsBunch->getdT());
        IpplTimings::stopTimer(timeIntegrationTimer1_m);

        // sets time for EnvelopeBunch
        itsBunch->actT();

        // FIXME: WHY DO I CALL THIS TWICE?
        // calculates new profile
        //itsBunch->updateFields();

        if(step % 1000 == 0)
            msg << " Step " << step << " at " << itsBunch->zAvg() << " [m] t= " << itsBunch->getT() << " [s] E=" << itsBunch->Eavg() * 1e-6 << " [MeV]" << endl;

        //t after a full global timestep with dT "synchronization point" for simulation time
        t += itsBunch->getdT();
        itsBunch->setT(t);

        double sposRef = itsBunch->get_sPos();
        if(step != 0 && (step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0))
            writePhaseSpace(step, sposRef);

        //  Stop simulation if beyond zstop_m
        if(sposRef > zstop_m) {
            maxSteps_m = step;
        }
    }

    OPAL.setLastStep(step);

    //itsOpalBeamline_m.switchElements(numeric_limits<double>::max(), numeric_limits<double>::min());
    itsOpalBeamline_m.switchElementsOff();

    msg << "done executing ParallelSliceTracker" << endl;
    //msg << "dumping stashed data to h5 file" << endl;
    //itsDataSink->dumpStashedPhaseSpaceEnvelope();
    //msg << "done dumping data to h5 file" << endl;
    writeLastStepPhaseSpace(step, itsBunch->get_sPos());
}

