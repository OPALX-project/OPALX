#ifdef HAVE_ENVELOPE_SOLVER
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

extern Inform* gmsg;
extern Inform* gmsg2all;

// Class ParallelSliceTracker
// ------------------------------------------------------------------------

ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        const PartData &reference,
        bool revBeam, 
        bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    itsOpalBeamline_m()
{
}


ParallelSliceTracker::ParallelSliceTracker(const Beamline &beamline,
        SLPartBunch &bunch,
        SLDataSink &ds,
        const PartData &reference,
        bool revBeam,
        bool revTrack,
        int maxSTEPS,
        double zstop):
    Tracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    maxSteps_m(maxSTEPS),
    zstop_m(zstop)
{ 
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


ParallelSliceTracker::~ParallelSliceTracker()
{
}

void ParallelSliceTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale)
{
}


void ParallelSliceTracker::applyExitFringe(double angle, double curve,
        const BMultipoleField &field, double scale)
{
}

// 2007/04/19 CKR
void ParallelSliceTracker::execute(FILE* BetOutputDat, FILE* BetOutputSli)
{
    //prepare output

    fstream fileZPOS;
    fileZPOS.open ("zpos.dat", fstream::out);

    double tmp;
    double recpgamma, gamma;
    double t = itsBunch->getT(); 
    double dt = itsBunch->getdT();
    double tEmission = itsBunch->getTEmission();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    long long step =  OPAL.inRestartRun() ? OPAL.getRestartStep() + 1 :lround(t/dt);

    Vector_t um, a, s;
    Vector_t externalE, externalB, K1, K2;
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

    *gmsg << "executing ParallelSliceTracker, initial DT " << itsBunch->getdT() 
        << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl
        << "the mass is: " << itsReference.getM()*1e-6 << " MeV, its charge: " << itsReference.getQ() << endl;

    itsBeamline_m.accept(*this);  // fill list allElements
    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(*gmsg);

    for(int i=0; i < itsBunch->getLocalNum(); i++) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));
        itsOpalBeamline_m.getSectionIndexAt(pos, l);
        itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l), itsOpalBeamline_m.getSectionStart(l));
    }

    itsBunch->calcBeamParameters();

    *gmsg << "ParallelEnvelopeTracker starting tracking..." << endl;

    for(step; step < maxSteps_m; ++step) {
        
        //check if any particle hasn't reached the end of the field from the last element
        global_EOL = true;
        bends = false;
        hasWake = false;
        wfSection = -1;
        double margin = 1e-6;

        itsOpalBeamline_m.resetStatus();
        itsBunch->calcBeamParameters();
        t = itsBunch->getT();

        //FIXME: this does not yet work
        //for (int i = 0; i < itsBunch->getLocalNumSlices(); ++i) {
        
        //external field for all slices
        for (int i = 0; i < itsBunch->getN(); i++) {

            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            K1        = Vector_t(0.0);
            K2        = Vector_t(0.0);

            Vector_t pos = Vector_t(itsBunch->getX(i), itsBunch->getY(i), itsBunch->getZ(i));

            long &ls = itsBunch->LastSection[i];
            itsOpalBeamline_m.getSectionIndexAt(pos, ls);
            
            //FIXME: find better solution for max(..., 0.0)
            //only check switch on/off elements for first and last slice
            if(i == 0 || i == itsBunch->getN()-1)
                itsOpalBeamline_m.switchElements(itsBunch->getZ(i) - margin, max(itsBunch->getZ(i),0.0) + margin);

            if(ls != itsBunch->LastSection[i]) {
                itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
                itsBunch->LastSection[i] = ls;
            }

            itsBunch->actT();

            const unsigned int& rtv = itsOpalBeamline_m.getFieldAt(i, pos, ls, t + itsBunch->dt[i]/2.0, externalE, externalB);
            global_EOL = global_EOL && (rtv & BEAMLINE_EOL);

            if((rtv & BEAMLINE_WAKE) && !hasWake) {
                wfSection = ls;
            }

            hasWake = hasWake || (rtv & BEAMLINE_WAKE);
            bends = bends || (rtv & BEAMLINE_BEND);

            // Calculate factors for the envelope equation
            const CompVec& elements = itsOpalBeamline_m.getSectionAt(pos,ls).getElements();
            for(CompVec::const_iterator elit = elements.begin(); elit != elements.end(); ++ elit) {
                (*elit)->addKR(i,t,K1);
                (*elit)->addKT(i,t,K2);
            }
     
            if(i==0) {
                FieldsOutput(itsBunch->getZ(i), externalE(0), externalE(1), externalE(2), externalB(0), externalB(1), externalB(2));
            }

            // pass K-values and E/B fields to EnvelopeBunch
            itsBunch->setKR(K1, i);
            itsBunch->setKT(K2, i);
            itsBunch->setBF(externalB, i);
            itsBunch->setEF(externalE, i);
        }
        
        //TEST the non-len reduce: reduce(&global_EOL, &global_EOL, OpBitwiseOrAssign());
        reduce(&global_EOL, &global_EOL + 1, &global_EOL, OpBitwiseAndAssign());
        if(global_EOL)
            break;

        /** do timestep for all slices
         *  calls the function BetBunch->Run() and solves the envelope equation 
         */
        // set cathode position to 0
        itsBunch->tstep(0, itsBunch->getdT(), step); 

        /** Write Output in .dat and .sli files */
        itsBunch->BetOut(BetOutputDat, BetOutputSli);

        //t after a full global timestep with dT "synchronization point" for simulation time
        t += itsBunch->getdT();

        itsBunch->setT(t);

        //double sposRef = 1.;
        //if (step % Options::psDumpFreq == 0 ) 
            //writePhaseSpace(sposRef);           
        
        /**
          Stop simulation if beyond zstop_m	   
          */
        double sposRef = itsBunch->getZ(itsBunch->getN()/2);
        if(sposRef > zstop_m) {
            maxSteps_m = step;
        }
    }


    fileZPOS.close();
    
    //itsOpalBeamline_m.switchElements(numeric_limits<double>::max(), numeric_limits<double>::min());
    itsOpalBeamline_m.switchElementsOff();

    *gmsg << "done executing ParallelSliceTracker" << endl;
}

void ParallelSliceTracker::FieldsOutput(double z, double Ex, double Ey, double Ez,  double Bx, double By, double Bz) {
    fstream output;
    output.open("testfields.dat",fstream::out|fstream::app);
    output << z << "\t" << Ex << "\t" << Ey << "\t" << Ez << "\t" <<  Bx << "\t"<<  By << "\t"<< Bz << endl;
    output.close();
}

#endif
