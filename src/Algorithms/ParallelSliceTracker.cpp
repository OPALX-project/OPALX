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
    Tracker(beamline, reference, revBeam, revTrack) {
    itsOpalBeamline_m = new OpalBeamline();
}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        EnvelopeBunch &bunch,
        DataSink &ds,
        const PartData &reference,
        bool revBeam,
        bool revTrack,
        int maxSTEPS,
        double zstop,
        ParallelTTracker &mySlApTracker):
    Tracker(beamline, reference, revBeam, revTrack),
    maxSteps_m(maxSTEPS),
    zstop_m(zstop) {
    itsBunch = &bunch;
    itsOpalBeamline_m = new OpalBeamline();

    mySlApTracker_m = & mySlApTracker;

    itsDataSink = &ds;
    scaleFactor_m = itsBunch->getdT() * c;

    timeIntegrationTimer1_m  = IpplTimings::getTimer("Time integration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("Time integration2");
    timeFieldEvaluation_m  = IpplTimings::getTimer("Field evaluation");

    BinRepartTimer_m   = IpplTimings::getTimer("Time of Binary repart.");
    WakeFieldTimer_m   = IpplTimings::getTimer("Time of Wake Field calc.");

}


ParallelSliceTracker::~ParallelSliceTracker() {
    delete itsOpalBeamline_m;
}


void ParallelSliceTracker::updateRFElement(string elName, double maxPhi) {
    /**
       The maximum phase is added to the nominal phase of
       the element. This is done on all nodes except node 0 where
       the Autophase took place.
    */
    FieldList cl = itsOpalBeamline_m->getElementByType("RFCavity");
    FieldList twl = itsOpalBeamline_m->getElementByType("TravelingWave");
    cl.merge(twl, OpalField::SortAsc);
    double phi = 0.0;

    for(FieldList::iterator fit = cl.begin(); fit != cl.end(); ++ fit) {
        if((*fit).getElement()->getName() == elName) {
            if((*fit).getElement()->getType() == "TravelingWave") {
                phi  =  static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
                phi += maxPhi;
                static_cast<TravelingWave *>((*fit).getElement())->setPhasem(phi);
            } else {
                phi  = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
                phi += maxPhi;
                static_cast<RFCavity *>((*fit).getElement())->setPhasem(phi);
            }
        }
    }
}

void ParallelSliceTracker::updateAllRFElements() {
    /**
       All RF-Elements gets updated, where the phiShift is the
       global phase shift in units of seconds.
    */

    FieldList cl = itsOpalBeamline_m->getElementByType("RFCavity");
    FieldList twl = itsOpalBeamline_m->getElementByType("TravelingWave");
    cl.merge(twl, OpalField::SortAsc);

    double phi;

    for(FieldList::iterator it = cl.begin(); it != cl.end(); ++ it) {
        if((*it).getElement()->getType() == "TravelingWave") {
            phi = static_cast<TravelingWave *>((*it).getElement())->getPhasem();
            phi = getCavityPhase(cavities_m, (*it).getElement()->getName());
            static_cast<TravelingWave *>((*it).getElement())->setPhasem(phi);
        } else {
            phi = static_cast<RFCavity *>((*it).getElement())->getPhasem();
            phi = getCavityPhase(cavities_m, (*it).getElement()->getName());
            static_cast<RFCavity *>((*it).getElement())->setPhasem(phi);
        }
    }
}


double ParallelSliceTracker::getCavityPhase(FieldList cav, string name) {
    double phi = 0.0;
    for(FieldList::iterator fit = cav.begin(); fit != cav.end(); ++ fit) {
        if((*fit).getElement()->getName() == name) {
            if((*fit).getElement()->getType() == "TravelingWave")
                phi = static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
            else
                phi = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
        }
    }
    return phi;
}


void ParallelSliceTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelSliceTracker::applyExitFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}

void ParallelSliceTracker::execute() {

    Inform msg("ParallelSliceTracker");
    double t = itsBunch->getT();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);
    long long step = 0.0;
    OpalData *OPAL = OpalData::getInstance();

    if(OPAL->inRestartRun()) {
        int prevDumpFreq = OPAL->getRestartDumpFreq();
        step = OPAL->getRestartStep() * prevDumpFreq + 1;
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    else if(OPAL->hasBunchAllocated() && Options::scan) {
        step = 1;
        if(!itsBunch->doEmission())
            writePhaseSpace(step - 1, 0.0); // write initial phase space
        itsBunch->setT(0.0);
    } else {
        step = OPAL->getLastStep() + 1;
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    *gmsg << "Track start at: t= " << itsBunch->getT() << " zstop@ " << zstop_m << " [m]" << endl;

    Vector_t um, a, s;
    Vector_t externalE, externalB, KR, KT;
    BorisPusher pusher(itsReference);
    Vector_t rmin, rmax;

    bool global_EOL;

    msg << "executing ParallelSliceTracker, initial DT " << itsBunch->getdT()
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM() * 1e-6 << " MeV, its charge: " << itsReference.getQ() << endl;

    itsBeamline_m.accept(*this);  // fill list allElements
    itsOpalBeamline_m->prepareSections();
    itsOpalBeamline_m->print(msg);

    for(int i = 0; i < itsBunch->getLocalNum(); i++) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));
        itsOpalBeamline_m->getSectionIndexAt(pos, l);
    }


    /*
      AAA DO AUTOPHASING
    */
    if(Options::autoPhase > 0 && !OPAL->hasBunchAllocated()) {
        cavities_m = mySlApTracker_m->executeAutoPhaseForSliceTracker();
        updateAllRFElements();
    } else if(Options::autoPhase > 0 && OPAL->hasBunchAllocated()) {
        for(std::vector<MaxPhasesT>::iterator it = OPAL->getFirstMaxPhases(); it < OPAL->getLastMaxPhases(); it++) {
            updateRFElement((*it).first, (*it).second);
        }
    }

    /**
       save autophase information in order to skip
       autophase in a restart run
    */

    if((!OPAL->inRestartRun()) && (Options::autoPhase > 0))
        itsDataSink->storeCavityInformation();


    msg << " is starting tracking..." << endl;

    for(; step < maxSteps_m; ++step) {

        //check if any particle hasn't reached the end of the field from the last element
        global_EOL = true;
        double margin = 1e-7;

        itsOpalBeamline_m->resetStatus();
        t = itsBunch->getT();
        itsOpalBeamline_m->switchElements(itsBunch->zTail() - margin, itsBunch->zHead() + margin);

        //external field for all local slices
        IpplTimings::startTimer(timeFieldEvaluation_m);
        for(int i = 0; i < itsBunch->getLocalNum(); i++) {

            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            KR        = Vector_t(0.0);
            KT        = Vector_t(0.0);

            //FIXME: why not x=y=0.0?
            Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));

            long &ls = itsBunch->LastSection[i];
            itsOpalBeamline_m->getSectionIndexAt(pos, ls);

            if(ls != itsBunch->LastSection[i])
                itsBunch->LastSection[i] = ls;

            //XXX: disregard itsBunch->dt
            unsigned long rtv = itsOpalBeamline_m->getFieldAt(i, pos, ls, t , externalE, externalB);
            global_EOL = global_EOL && (rtv & BEAMLINE_EOL);

            //XXX: disregard itsBunch->dt
            // Calculate factors for the envelope equation
            itsOpalBeamline_m->getKFactors(i, pos, ls, t, KR, KT);

            itsBunch->setExternalFields(i, externalE, externalB, KR, KT);
        }
        IpplTimings::stopTimer(timeFieldEvaluation_m);

        //TEST the non-len reduce: reduce(&global_EOL, &global_EOL, OpBitwiseOrAssign());
        //reduce(&global_EOL, &global_EOL + 1, &global_EOL, OpBitwiseAndAssign());
        MPI_Allreduce(MPI_IN_PLACE, &global_EOL, 1, MPI_INT, MPI_LAND, Ippl::getComm());
        if((bool)global_EOL)
            break;

        itsBunch->computeSpaceCharge();

        IpplTimings::startTimer(timeIntegrationTimer1_m);
        // solve the envelope equation for all active slices
        itsBunch->timeStep(itsBunch->getdT());
        IpplTimings::stopTimer(timeIntegrationTimer1_m);

        // sets time for EnvelopeBunch
        itsBunch->actT();

        if(step % 1000 == 0)
            msg << " Step " << step << " at " << itsBunch->zAvg() << " [m] t= "
                << itsBunch->getT() << " [s] E=" << itsBunch->Eavg() * 1e-6
                << " [MeV]" << endl;

        //t after a full global timestep with dT "synchronization point" for simulation time
        t += itsBunch->getdT();
        itsBunch->setT(t);

        double sposRef = itsBunch->get_sPos();
        if(step != 0 && (step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0))
            writePhaseSpace(step, sposRef);

        //  Stop simulation if beyond zstop_m or bunch invalid
        if(sposRef > zstop_m || !(itsBunch->isValid_m)) {
            maxSteps_m = step;
        }
    }

    OPAL->setLastStep(step);

    //dump last phase space
    writeLastStepPhaseSpace(step, itsBunch->get_sPos());

    //itsOpalBeamline_m->switchElements(numeric_limits<double>::max(), numeric_limits<double>::min());
    itsOpalBeamline_m->switchElementsOff();

    msg << "done executing ParallelSliceTracker" << endl;


}

