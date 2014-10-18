// ------------------------------------------------------------------------
// $RCSfile: ParallelTTracker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParallelTTracker
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
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/PartPusher.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Corrector.h"
#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/Marker.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Probe.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/Separator.h"
#include "AbsBeamline/Septum.h"
#include "AbsBeamline/Solenoid.h"
#include "AbsBeamline/ParallelPlate.h"
#include "AbsBeamline/CyclotronValley.h"
#include "Beamlines/Beamline.h"
#include "Lines/Sequence.h"
//--------- Added by Xiaoying Pang 04/22/2014 ---------------
#include "Solvers/CSRWakeFunction.hh"

#include "AbstractObjects/OpalData.h"

#include "BasicActions/Option.h"

#include "Distribution/Distribution.h"
#include "ValueDefinitions/RealVariable.h"
#include "Utilities/Timer.h"
#include "Utilities/OpalException.h"
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Structure/BoundaryGeometry.h"

class PartData;

using namespace std;
using namespace OPALTimer;

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
				   size_t N):
Tracker(beamline, reference, revBeam, revTrack),
itsBunch(NULL),
itsDataSink_m(NULL),
specifiedNPart_m(N),
bgf_m(NULL),
itsOpalBeamline_m(),
lineDensity_m(),
RefPartR_zxy_m(0.0),
RefPartP_zxy_m(0.0),
RefPartR_suv_m(0.0),
RefPartP_suv_m(0.0),
globalEOL_m(false),
wakeStatus_m(false),
//--------- Added by Xiaoying Pang 04/22/2014 ---------------
wakeFunction_m(NULL),
surfaceStatus_m(false),
secondaryFlg_m(false),
mpacflg_m(true),
nEmissionMode_m(false),
zStop_m(0.0),
scaleFactor_m(1.0),
vscaleFactor_m(scaleFactor_m),
recpGamma_m(1.0),
rescale_coeff_m(1.0),
dtTrack_m(0.0),
surfaceEmissionStop_m(-1),
minStepforReBin_m(-1),
minBinEmitted_m(std::numeric_limits<size_t>::max()),
repartFreq_m(-1),
lastVisited_m(-1),
numRefs_m(-1),
gunSubTimeSteps_m(-1),
emissionSteps_m(std::numeric_limits<unsigned int>::max()),
localTrackSteps_m(0),
maxNparts_m(0),
numberOfFieldEmittedParticles_m(std::numeric_limits<size_t>::max()),
bends_m(0),
numParticlesInSimulation_m(0),
space_orientation_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
timeFieldEvaluation_m(IpplTimings::getTimer("Fieldeval")),
BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
Nimpact_m(0),
SeyNum_m(0.0),
sphys_m(NULL){
}

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunch &bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   int maxSTEPS,
                                   double zstop,
                                   int timeIntegrator,
				   size_t N):
Tracker(beamline, reference, revBeam, revTrack),
itsBunch(&bunch),
specifiedNPart_m(N),
itsDataSink_m(&ds),
bgf_m(NULL),
itsOpalBeamline_m(),
lineDensity_m(),
RefPartR_zxy_m(0.0),
RefPartP_zxy_m(0.0),
RefPartR_suv_m(0.0),
RefPartP_suv_m(0.0),
globalEOL_m(false),
wakeStatus_m(false),
surfaceStatus_m(false),
secondaryFlg_m(false),
mpacflg_m(true),
nEmissionMode_m(false),
zStop_m(zstop),
scaleFactor_m(itsBunch->getdT() * Physics::c),
vscaleFactor_m(scaleFactor_m),
recpGamma_m(1.0),
rescale_coeff_m(1.0),
dtTrack_m(0.0),
surfaceEmissionStop_m(-1),
minStepforReBin_m(-1),
minBinEmitted_m(std::numeric_limits<size_t>::max()),
repartFreq_m(-1),
lastVisited_m(-1),
numRefs_m(-1),
gunSubTimeSteps_m(-1),
emissionSteps_m(numeric_limits<unsigned int>::max()),
localTrackSteps_m(maxSTEPS),
maxNparts_m(0),
numberOfFieldEmittedParticles_m(numeric_limits<size_t>::max()),
bends_m(0),
numParticlesInSimulation_m(0),
space_orientation_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
timeFieldEvaluation_m(IpplTimings::getTimer("Fieldeval")),
BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
timeIntegrator_m(timeIntegrator),
Nimpact_m(0),
SeyNum_m(0.0) {

    //    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());

}

#ifdef HAVE_AMR_SOLVER
ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunch &bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   int maxSTEPS,
                                   double zstop,
                                   int timeIntegrator,
				   size_t N,
				   Amr* amrptr_in):
Tracker(beamline, reference, revBeam, revTrack),
itsBunch(&bunch),
specifiedNPart_m(N),
itsDataSink_m(&ds),
bgf_m(NULL),
itsOpalBeamline_m(),
lineDensity_m(),
RefPartR_zxy_m(0.0),
RefPartP_zxy_m(0.0),
RefPartR_suv_m(0.0),
RefPartP_suv_m(0.0),
globalEOL_m(false),
wakeStatus_m(false),
surfaceStatus_m(false),
secondaryFlg_m(false),
mpacflg_m(true),
nEmissionMode_m(false),
zStop_m(zstop),
scaleFactor_m(itsBunch->getdT() * Physics::c),
vscaleFactor_m(scaleFactor_m),
recpGamma_m(1.0),
rescale_coeff_m(1.0),
dtTrack_m(0.0),
surfaceEmissionStop_m(-1),
minStepforReBin_m(-1),
minBinEmitted_m(std::numeric_limits<size_t>::max()),
repartFreq_m(-1),
lastVisited_m(-1),
numRefs_m(-1),
gunSubTimeSteps_m(-1),
emissionSteps_m(numeric_limits<unsigned int>::max()),
localTrackSteps_m(maxSTEPS),
maxNparts_m(0),
numberOfFieldEmittedParticles_m(numeric_limits<size_t>::max()),
bends_m(0),
numParticlesInSimulation_m(0),
space_orientation_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
timeFieldEvaluation_m(IpplTimings::getTimer("Fieldeval")),
BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
timeIntegrator_m(timeIntegrator),
Nimpact_m(0),
SeyNum_m(0.0),
amrptr(amrptr_in){
}
#endif

ParallelTTracker::~ParallelTTracker() {

}

void ParallelTTracker::applySchottkyCorrection(PartBunch &itsBunch, int ne, double t, double rescale_coeff) {

    Inform msg("ParallelTTracker ");
    const long ls = 0;
    /*
     Now I can calculate E_{rf} at each position
     of the newely generated particles and rescale Q

     Note:
     For now I only sample the field of the last emitted particles.
     Space charge is not yet included
     */


    double laser_erg = itsBunch.getLaserEnergy(); // 4.7322; energy of single photon of 262nm laser  [eV]
    double workFunction = itsBunch.getWorkFunctionRf(); // espace energy for copper (4.31)  [eV]

    const double schottky_coeff = 0.037947; // coeffecient for calculate schottky potenial from E field [eV/(MV^0.5)]

    if(ne == 0)
        return ;

    double Ez = 0;
    double obtain_erg = 0;
    double par_t = 0;
    for(int k = 0; k < ne; k++) {
        size_t n = itsBunch.getLocalNum() - k - 1;
        Vector_t externalE(0.0);
        Vector_t externalB(0.0);

        itsBunch.R[n] *= Vector_t(Physics::c * itsBunch.dt[n]);
        par_t = t + itsBunch.dt[n] / 2;
        itsOpalBeamline_m.getFieldAt(n, itsBunch.R[n], ls, par_t, externalE, externalB);
        Ez = externalE(2);

        // fabs(Ez): if the field of cathode surface is in the right direction, it will increase the
        // energy which electron obtain. If the field is in the wrong direction, this particle will
        // be back to the cathode surface and then be deleted automaticly by OPAL,  we don't add
        // another logical branch to handle this. So fabs is the simplest way to handle this
        obtain_erg = laser_erg - workFunction + schottky_coeff * sqrt(fabs(Ez) / 1E6);
        double schottkyScale = obtain_erg * obtain_erg * rescale_coeff;

        itsBunch.Q[n] *= schottkyScale;
        itsBunch.R[n] /= Vector_t(Physics::c * itsBunch.dt[n]);
    }
}

double ParallelTTracker::schottkyLoop(double rescale_coeff) {

    Inform msg("ParallelTTracker ");

    double recpgamma;
    double t = 0.0;
    double dt = itsBunch->getdT();
    Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    unsigned long long step = 0;
    unsigned int emissionSteps = 0;

    Vector_t um, a, s;
    Vector_t externalE, externalB;
    BorisPusher pusher(itsReference);
    Vector_t rmin, rmax;

    bool global_EOL;

    bool hasSwitchedToTEmission = false;
    bool hasSwitchedBackToTTrack = false;

    size_t totalParticles_i = itsBunch->getTotalNum();

    msg << "*****************************************************************" << endl;
    msg << " Estimate Schottky correction                                    " << endl;
    msg << "*****************************************************************" << endl;

    double margin = 0.0;
    if(!mpacflg_m) {
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            long &l = itsBunch->LastSection[i];
            l = -1;
            itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
            itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l), itsOpalBeamline_m.getSectionStart(l));
        }

        if(!(itsBunch->WeHaveEnergyBins())) {
            IpplTimings::startTimer(BinRepartTimer_m);
            itsBunch->do_binaryRepart();
            IpplTimings::stopTimer(BinRepartTimer_m);
            Ippl::Comm->barrier();
        }

        // Check if there are any particles in simulation. If there are,
        // as in a restart, use the usual function to calculate beam
        // parameters. If not, calculate beam parameters of the initial
        // beam distribution.
        if(totalParticles_i == 0) {// fixme: maybe cause nonsense output if initialized momenta=0; Q: by Chuan.
            itsBunch->calcBeamParametersInitial();
        } else {
            itsBunch->calcBeamParameters();
        }

        RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
        RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

        if(!OpalData::getInstance()->hasBunchAllocated()) {
            updateSpaceOrientation(false);  // vec{p} = (0,0,p_z), vec{r} = (0,0,z)
        }

        RefPartR_suv_m = itsBunch->get_rmean();
        RefPartP_suv_m = itsBunch->get_pmean();
        /* Activate all elements which influence the particles when the simulation starts;
         *  mark all elements which are already past.
         */

        /*
         increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
         before they are allocated when increasing the timestep in the gun.
         */
        itsBunch->get_bounds(rmin, rmax);
        margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + pSqr(RefPartP_suv_m, RefPartP_suv_m));
        margin = 0.01 > margin ? 0.01 : margin;
        itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);
    }

    double minBinEmitted  = 10.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar) {
        minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
        // This variable controls the number of particles to be emitted before we use
        // the space charge solver.
        msg << "MINBINEMITTED " << minBinEmitted << endl;
    }


    double minStepforReBin  = 10000.0;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br) {
        minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        msg << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

    int repartFreq = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep) {
        repartFreq = static_cast<int>(rep->getReal());  // this variable controls the minimal number of steps until we repartition the particles
        msg << "REPARTFREQ " << repartFreq << endl;
    }

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq = 1000000;

    size_t totalParticles_f = 0;

    for(; step < localTrackSteps_m; ++step) {
        global_EOL = true;  // check if any particle hasn't reached the end of the field from the last element

        itsOpalBeamline_m.resetStatus();

        IpplTimings::startTimer(timeIntegrationTimer1_m);

        // reset E and B to Vector_t(0.0) for every step

        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);

        Nimpact_m = 0; // Initial parallel plate benchmark variable.
        SeyNum_m = 0; // Initial parallel plate benchmark variable.

        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            //scale each particle with c*dt
            itsBunch->R[i] /= vscaleFactor;
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            // update local coordinate system of particleInform &PartBunc
            itsBunch->X[i] /= vscaleFactor;
            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
                        itsBunch->getdT());
            itsBunch->X[i] *= vscaleFactor;
        }

        if(totalParticles_i > minBinEmitted) {
            itsBunch->boundp();
        }

        IpplTimings::stopTimer(timeIntegrationTimer1_m);

        itsBunch->calcBeamParameters();


        /** \f[ Space Charge  \f]
         */
        if(itsBunch->hasFieldSolver() && totalParticles_i > minBinEmitted && fabs(itsBunch->getChargePerParticle()) > 0.0) {
            // Do repartition if we have enough particles.
            if(totalParticles_i > 1000 && (((step + 1) % repartFreq) == 0)) {
                INFOMSG("*****************************************************************" << endl);
                INFOMSG("do repartition because of repartFreq" << endl);
                INFOMSG("*****************************************************************" << endl);
                IpplTimings::startTimer(BinRepartTimer_m);
                itsBunch->do_binaryRepart();
                IpplTimings::stopTimer(BinRepartTimer_m);
                Ippl::Comm->barrier();
                INFOMSG("*****************************************************************" << endl);
                INFOMSG("do repartition done" << endl);
                INFOMSG("*****************************************************************" << endl);
            }

            // Calculate space charge.
            if(itsBunch->WeHaveEnergyBins()) {
                // When we have energy bins.
                itsBunch->calcGammas();
                ParticleAttrib<double> Q_back = itsBunch->Q;
                itsBunch->resetInterpolationCache();
                for(int binNumber = 0; binNumber <= itsBunch->getLastemittedBin() && binNumber < itsBunch->getNumBins(); ++binNumber) {
                    itsBunch->setBinCharge(binNumber);
                    itsBunch->computeSelfFields(binNumber);
                    itsBunch->Q = Q_back;
                }
            } else {
                // When we don't.
                itsBunch->computeSelfFields();
                /**
                 Need this maybe for the adaptive time integration scheme
                 pair<Vector_t,Vector_t> eExtrema = itsBunch->getEExtrema();
                 INFOMSG("maxE= " << eExtrema.first << " minE= " << eExtrema.second << endl);
                 */
            }
        }

        IpplTimings::startTimer(timeIntegrationTimer2_m);


        /*
         transport and emit particles
         that passed the cathode in the first
         half-step or that would pass it in the
         second half-step.

         to make IPPL and the field solver happy
         make sure that at least 10 particles are emitted

         also remember that node 0 has
         all the particles to be emitted

         this has to be done *after* the calculation of the
         space charges! thereby we neglect space charge effects
         in the very first step of a new-born particle.

         */

        if((itsBunch->WeHaveEnergyBins())) {

            // switch to TEmission
            if(!hasSwitchedToTEmission) {
                dt = itsBunch->GetEmissionDeltaT();
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                msg << "Changing emission time step to: " << dt << endl;
                hasSwitchedToTEmission = true;
            }

            int ne = 0;
            Vector_t externalE = Vector_t(0.0);
            Vector_t externalB = Vector_t(0.0);
            itsOpalBeamline_m.getFieldAt(Vector_t(0.0),
                                         Vector_t(0.0),
                                         itsBunch->getT() + 0.5 * itsBunch->getdT(),
                                         externalE,
                                         externalB);
            ne += itsBunch->EmitParticles(externalE[2]);

            if(Options::schottkyCorrection && !hasSwitchedBackToTTrack)
                applySchottkyCorrection(*itsBunch, ne, t, rescale_coeff);

            reduce(ne, ne, OpAddAssign());
            totalParticles_i += ne;

            //emission has finished, reset to TTrack
            if(itsBunch->getNumBins() == itsBunch->getLastemittedBin() &&
               !hasSwitchedBackToTTrack) {
                //dt = dtTrack;
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                msg << "Emission done. Switching back to track timestep: " << dt << endl;
                hasSwitchedBackToTTrack = true;
                break;
            }

        } else {
            //emission has finished, reset to TTrack
            if(!hasSwitchedBackToTTrack) {
                //dt = dtTrack;
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                msg << "Emission done. Switching back to track timestep: " << dt << endl;
                hasSwitchedBackToTTrack = true;
            }

        }

        // push the reference particle by a half step
        recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;

        //
        // get external fields for all particles
        //
        IpplTimings::startTimer(timeFieldEvaluation_m);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            //FIXME: rethink scaling!
            itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

            long ls = itsBunch->LastSection[i];
            itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
            if(ls != itsBunch->LastSection[i]) {
                if(!itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls)) {
                    itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
                }
                itsBunch->LastSection[i] = ls;
            }
            const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls, t + itsBunch->dt[i] / 2., externalE, externalB);

            global_EOL = global_EOL && (rtv & BEAMLINE_EOL);

            // skip rest of the particle push if the
            // particle is out of bounds i.e. does not see
            // a E or B field
            if(rtv & BEAMLINE_OOB)
                itsBunch->Bin[i] = -1;


            itsBunch->Ef[i] += externalE;
            itsBunch->Bf[i] += externalB;

            itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

            // in case a particle is pushed behind the emission surface, delete the particle

            if(itsBunch->R[i](2) < 0)
                itsBunch->Bin[i] = -1;

        }

        IpplTimings::stopTimer(timeFieldEvaluation_m);

        //        if(itsBunch->getLocalNum() == 0)
        //    global_EOL = false;

        /**
         Delete marked particles.
         */

        bool globPartOutOfBounds = (min(itsBunch->Bin) < 0);
        size_t ne = 0;
        if(globPartOutOfBounds) {
            ne = itsBunch->boundp_destroyT();
        }

        totalParticles_f = totalParticles_i - ne;
        if(ne > 0)
            msg << "* Deleted in Shotky " << ne << " particles, remaining " << totalParticles_f << " particles" << endl; //benchmark output

        kickParticles(pusher);

        if(totalParticles_f > 0) {
            // none of the particles is in a bending element
            updateReferenceParticle();
        }

        itsBunch->RefPart_R = RefPartR_zxy_m;
        itsBunch->RefPart_P = RefPartP_zxy_m;

        // calculate the dimensions of the bunch and add a small margin to them; then decide which elements have to be triggered
        // when an element is triggered memory is allocated and the field map is read in
        itsBunch->get_bounds(rmin, rmax);

        // trigger the elements
        margin = 3. * RefPartP_suv_m(2) * recpgamma;
        margin = 0.01 > margin ? 0.01 : margin;
        itsOpalBeamline_m.switchElements((rmin(2) - margin)*scaleFactor_m, (rmax(2) + margin)*scaleFactor_m);

        // start normal particle loop part 2 for simulation without boundary geometry.
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f]
             * \code
             * R[i] += 0.5 * P[i] * recpgamma;
             * \endcode
             */
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            //and scale back to dimensions
            itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
            // update local coordinate system
            itsBunch->X[i] /= vscaleFactor;
            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
            itsBunch->X[i] *= vscaleFactor;
            //reset time step if particle was emitted in the first half-step
            //the particle is now in sync with the simulation timestep
            itsBunch->dt[i] = itsBunch->getdT();
        }

        IpplTimings::stopTimer(timeIntegrationTimer2_m);

        if(totalParticles_f > minBinEmitted)
            itsBunch->boundp();

        totalParticles_i = itsBunch->getTotalNum();


        t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time

        itsBunch->setT(t);

        //IFF: cheap step dump regulation
        OPALTimer::Timer myt2;
        double sposRef = 0.0;
        if(totalParticles_f > 0) {
            sposRef = itsBunch->get_sPos();
            if(totalParticles_f <= minBinEmitted) {
                INFOMSG(myt2.time() << " Step " << step << "; only " << totalParticles_f << " particles emitted; t= " << t
                        << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl);
            } else if(std::isnan(sposRef) || std::isinf(sposRef)) {
                INFOMSG(myt2.time() << " Step " << step << "; there seems to be something wrong with the position of the bunch!" << endl);
            } else {
                INFOMSG(myt2.time() << " Step " << step << " at " << sposRef << " [m] t= " << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl);
                if(step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0) {
                    size_t nLoc = itsBunch->getLocalNum();
                    reduce(nLoc, nLoc, OpMultipplyAssign());
                    if((nLoc == 0) || ((step + 1) % repartFreq == 0)) {
                        INFOMSG("*****************************************************************" << endl);
                        INFOMSG("do repartition because of zero particles or repartition frequency" << endl);
                        IpplTimings::startTimer(BinRepartTimer_m);
                        itsBunch->do_binaryRepart();
                        IpplTimings::stopTimer(BinRepartTimer_m);
                        Ippl::Comm->barrier();
                        INFOMSG("done" << endl);
                        INFOMSG("*****************************************************************" << endl);
                    }
                }
            }
            /**
             Stop simulation if beyond zStop_m
             */
            if(sposRef > zStop_m) {
                localTrackSteps_m = step;
            }
        } else {
            INFOMSG("Step " << step << " no emission yet "  << " t= " << t << " [s]" << endl);
        }

        if(step > emissionSteps) {
            reduce(&global_EOL, &global_EOL + 1, &global_EOL, OpBitwiseAndAssign());
            if(global_EOL) {
                break;
            }
        }
        // this seams to fix Ticket #12
        //  Ippl::Comm->barrier();
        itsBunch->get_bounds(rmin, rmax);
        // trigger the elements
        RefPartP_suv_m = itsBunch->get_pmean();
        recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

        margin = 10. * RefPartP_suv_m(2) * recpgamma * scaleFactor_m;
        margin = 0.01 > margin ? 0.01 : margin;
    }
    OPALTimer::Timer myt3;
    OpalData::getInstance()->setLastStep(step);
    msg << "done executing Schottky loop " << myt3.time() << endl;
    return itsBunch->getCharge();
}


void ParallelTTracker::checkCavity(double s, Component *& comp, double & cavity_start_pos) {

    comp = NULL;
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if((fit != currently_ap_cavity_m)
           && ((*fit).getStart() <= s) && (s <= (*fit).getEnd())) {
            comp = (*fit).getElement();
            cavity_start_pos = (*fit).getStart();
            currently_ap_cavity_m = fit;
            return;
        }
    }
}


void ParallelTTracker::doOneStep(BorisPusher & pusher) {

    bends_m = 0;

    double t = itsBunch->getT();

    // increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
    // before they are allocated when increasing the timestep in the gun.

    switchElements(10.0);

    itsOpalBeamline_m.resetStatus();

    timeIntegration1(pusher);

    itsBunch->calcBeamParameters();

    // reset E and B to Vector_t(0.0) for every step
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    computeExternalFields();

    timeIntegration2(pusher);

    t += itsBunch->getdT();
    itsBunch->setT(t);
}


void ParallelTTracker::applyEntranceFringe(double angle, double curve,
                                           const BMultipoleField &field, double scale) {
}


void ParallelTTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale) {
}


void ParallelTTracker::showCavities(Inform &msg) {

    msg << "Found the following cavities:" << endl;

    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        msg << (*fit).getElement()->getName()
        << " from " << (*fit).getStart() << " to "
        << (*fit).getEnd() << " (m) phi=";
        if((*fit).getElement()->getType() == "TravelingWave")
            msg << static_cast<TravelingWave *>((*fit).getElement())->getPhasem() / Physics::pi * 180.0 << endl;
        else
            msg << static_cast<RFCavity *>((*fit).getElement())->getPhasem() / Physics::pi * 180.0 << endl;
    }
    msg << endl << endl;
}


void ParallelTTracker::updateRFElement(std::string elName, double maxPhi) {
    /**
     The maximum phase is added to the nominal phase of
     the element. This is done on all nodes except node 0 where
     the Autophase took place.
     */
    double phi  = 0.0;
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if((*fit).getElement()->getName() == elName) {
            if((*fit).getElement()->getType() == "TravelingWave") {
                phi  =  static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
                phi += maxPhi;
                static_cast<TravelingWave *>((*fit).getElement())->updatePhasem(phi);
            } else {
                phi  = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
                phi += maxPhi;
                static_cast<RFCavity *>((*fit).getElement())->updatePhasem(phi);
            }
        }
    }
}


void ParallelTTracker::updateAllRFElements(double phiShift) {
    /**
     All RF-Elements gets updated, where the phiShift is the
     global phase shift in units of seconds.
     */
    Inform msg("ParallelTTracker ");
    double phi = 0;
    double freq = 0.0;
    const double RADDEG = 1.0 / Physics::pi * 180.0;
    //   Inform m ("updateALLRFElements ",INFORM_ALL_NODES);
    msg << "\n-------------------------------------------------------------------------------------\n";
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if(fit != cavities_m.begin())
            msg << "\n";
        if((*fit).getElement()->getType() == "TravelingWave") {
            freq = static_cast<TravelingWave *>((*fit).getElement())->getFrequencym();
            phi = static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
            msg << (*fit).getElement()->getName()
            << ": phi= phi_nom + phi_maxE + global phase shift= " << (phi*RADDEG)-(phiShift*freq*RADDEG) << " degree, "
            << "(global phase shift= " << -phiShift *freq *RADDEG << " degree)\n";
            phi -= (phiShift * freq);
            static_cast<TravelingWave *>((*fit).getElement())->updatePhasem(phi);
        } else {
            freq = static_cast<RFCavity *>((*fit).getElement())->getFrequencym();
            phi = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
            msg << (*fit).getElement()->getName()
            << ": phi= phi_nom + phi_maxE + global phase shift= " << (phi*RADDEG)-(phiShift*freq*RADDEG) << " degree, "
            << "global phase shift= " << -phiShift *freq *RADDEG << " degree\n";
            phi -= (phiShift * freq);
            static_cast<RFCavity *>((*fit).getElement())->updatePhasem(phi);
        }
    }
    msg << "-------------------------------------------------------------------------------------\n"
	<< endl;
}


FieldList ParallelTTracker::executeAutoPhaseForSliceTracker() {
    Inform msg("executeAutoPhaseForSliceTracker ");

    double gPhaseSave;

    gPhaseSave = OpalData::getInstance()->getGlobalPhaseShift();
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    itsBeamline_m.accept(*this);
    // make sure that no monitor has overlap with two tracks
    FieldList monitors = itsOpalBeamline_m.getElementByType("Monitor");
    for(FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
        double zbegin, zend;
        it->getElement()->getDimensions(zbegin, zend);
        if(zbegin < zStop_m && zend >= zStop_m) {
            msg << "\033[0;31m"
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
            << "% Removing '" << it->getElement()->getName() << "' since it resides in two tracks.   %\n"
            << "% Please adjust zstop or place your monitor at a different position to prevent this. %\n "
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
            << "\033[0m"
            << endl;
            static_cast<Monitor *>(it->getElement())->moveBy(-zend - 0.001);
            itsOpalBeamline_m.removeElement(it->getElement()->getName());
        }
    }
    itsOpalBeamline_m.prepareSections();

    cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
    currently_ap_cavity_m = cavities_m.end();
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType("TravelingWave");
    cavities_m.merge(travelingwaves, OpalField::SortAsc);

    int tag = 101;
    int Parent = 0;
    Vector_t iniR(0.0);
    Vector_t iniP(0.0, 0.0, 1E-6);
    PID_t id;
    Ppos_t r, p, x;
    ParticleAttrib<double> q, dt;
    ParticleAttrib<int> bin;
    ParticleAttrib<long> ls;
    ParticleAttrib<short> ptype;

    double zStop = itsOpalBeamline_m.calcBeamlineLenght();

    msg << "Preparation done zstop= " << zStop << endl;

    if(Ippl::myNode() == 0)
      itsBunch->create(1);
    itsBunch->update();

    if(Ippl::myNode() == 0) {
      itsBunch->R[0] = iniR;
      itsBunch->P[0] = iniP;
      itsBunch->Bin[0] = 0;
      itsBunch->Q[0] = itsBunch->getChargePerParticle();
      itsBunch->PType[0] = 0;
      itsBunch->LastSection[0] = 0;

      RefPartP_suv_m = iniP;
      RefPartR_suv_m = iniR;
    }
    updateSpaceOrientation(false);
    executeAutoPhase(Options::autoPhase, zStop);
    if(Ippl::myNode() == 0) {
      RefPartP_suv_m = itsBunch->P[0];
      // need to rebuild for updateAllRFElements
      cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
      travelingwaves = itsOpalBeamline_m.getElementByType("TravelingWave");
      cavities_m.merge(travelingwaves, OpalField::SortAsc);

      // now send all max phases and names of the cavities to
      // all the other nodes for updating.
      Message *mess = new Message();
      putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

      for(std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
	putMessage(*mess, (*it).first);
	putMessage(*mess, (*it).second);
      }
      Ippl::Comm->broadcast_all(mess, tag);
    } else {
      // receive max phases and names and update the structures
      int nData = 0;
      Message *mess = Ippl::Comm->receive_block(Parent, tag);
      getMessage(*mess, nData);
      for(int i = 0; i < nData; i++) {
	std::string elName;
	double maxPhi;
	getMessage(*mess, elName);
	getMessage(*mess, maxPhi);
	updateRFElement(elName, maxPhi);
	OpalData::getInstance()->setMaxPhase(elName, maxPhi);
      }
    }

    if(Ippl::myNode() == 0)
      itsBunch->destroy(1, 0);
    itsBunch->update();

    OpalData::getInstance()->setGlobalPhaseShift(gPhaseSave);
    return cavities_m;
}


void ParallelTTracker::executeAutoPhase(int numRefs, double zStop) {
    Inform msg("Autophasing ");

    const double RADDEG = 180.0 / Physics::pi;

    Vector_t rmin, rmax;

    size_t maxStepsSave = localTrackSteps_m;
    size_t step = 0;

    double tSave = itsBunch->getT();
    double dTSave = itsBunch->getdT();
    double scaleFactorSave = scaleFactor_m;

    int dtfraction = 2;
    double newDt = itsBunch->getdT() / dtfraction;
    itsBunch->setdT(newDt);
    itsBunch->dt = newDt;

    bends_m = 0;  //fixme: AP and bends will not work at the moment

    scaleFactor_m = itsBunch->dt[0] * Physics::c;
    vscaleFactor_m = Vector_t(scaleFactor_m);

    BorisPusher pusher(itsReference);

    msg << "\n"
    << "start at t= " << itsBunch->getT() << " [s], zstop at: "
    << zStop << " [m], Nplocal= " << itsBunch->getLocalNum() << "\n"
    << "initial DT " << itsBunch->dt[0] << " [s], step= "
    << step << ", R =  " << itsBunch->R[0] << " [m]" << endl;

    //    showCavities(m);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
        itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l),
                                             itsOpalBeamline_m.getSectionStart(l));
    }

    RefPartR_suv_m = RefPartR_zxy_m = rmin = rmax = itsBunch->R[0];
    RefPartP_suv_m = RefPartP_zxy_m = itsBunch->P[0];

    /* Activate all elements which influence the particles when the simulation starts;
     * mark all elements which are already past.
     *
     * Increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
     * before they are allocated when increasing the timestep in the gun.
     */

    double margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

    margin = 0.01 > margin ? 0.01 : margin;

    itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin, true);

    double cavity_start = 0.0;
    Component *cavity = NULL;

    for(; step < localTrackSteps_m * dtfraction; ++step) {

        if (itsBunch->getLocalNum() != 0) {
            // let's do a drifting step to probe if the particle will reach element in next step
            Vector_t R_drift = itsBunch->R[0] + itsBunch->P[0] / sqrt(1.0 + dot(itsBunch->P[0],
                                                                                itsBunch->P[0])) * vscaleFactor_m;
            checkCavity(R_drift[2], cavity, cavity_start);
        }
        else
            cavity = NULL;

        if(cavity != NULL) {
            double orig_phi = 0.0;
            double Phimax = 0.0, Emax = 0.0;
            double PhiAstra;
            //////
            const double beta = sqrt(1. - 1 / (itsBunch->P[0](2) * itsBunch->P[0](2) + 1.));
            const double tErr  = (cavity_start - itsBunch->R[0](2)) / (Physics::c * beta);

            bool apVeto;

            INFOMSG("Found " << cavity->getName()
                    << " at " << itsBunch->R[0](2) << " [m], "
                    << "step  " << step << ", "
                    << "t= " << itsBunch->getT() << " [s],\n"
                    << "E= " << getEnergyMeV(itsBunch->P[0]) << " [MeV]\n"
                    << "start phase scan ... " << endl);


            INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
            if(cavity->getType() == "TravelingWave") {
                orig_phi = static_cast<TravelingWave *>(cavity)->getPhasem();
                apVeto = static_cast<TravelingWave *>(cavity)->getAutophaseVeto();
                if(apVeto) {
                    msg << " ----> APVETO -----> "
                    << static_cast<TravelingWave *>(cavity)->getName() <<  endl;
                    Phimax = orig_phi;
                }
                INFOMSG(cavity->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);

                if(!apVeto) {
                    TravelingWave *element = static_cast<TravelingWave *>(cavity);
                    Phimax = element->getAutoPhaseEstimate(getEnergyMeV(itsBunch->P[0]),
                                                           itsBunch->getT() + tErr,
                                                           itsReference.getQ(),
                                                           itsReference.getM() * 1e-6);
                }
            } else {
                orig_phi = static_cast<RFCavity *>(cavity)->getPhasem();
                apVeto = static_cast<RFCavity *>(cavity)->getAutophaseVeto();
                if(apVeto) {
                    msg << " ----> APVETO -----> "
                    << static_cast<RFCavity *>(cavity)->getName() << endl;
                    Phimax = orig_phi;
                }
                INFOMSG(cavity->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);

                if(!apVeto) {
                    RFCavity *element = static_cast<RFCavity *>(cavity);
                    Phimax = element->getAutoPhaseEstimate(getEnergyMeV(itsBunch->P[0]),
                                                           itsBunch->getT() + tErr,
                                                           itsReference.getQ(),
                                                           itsReference.getM() * 1e-6);
                }
            }

            double Phiini = Phimax;
            double phi = Phiini;
            double dphi = Physics::pi / 360.0;
            int j = -1;

            double E = APtrack(cavity, cavity_start, phi);
            if(!apVeto) {
                msg << "Did APtrack with phi= " << phi << " result E= " << E << endl;
                //                INFOMSG("do fine scan around effective max energy (" << E << " MeV)" << ", dphi= " << dphi << endl;);
                do {
                    j ++;
                    Emax = E;
                    Phiini = phi;
                    phi -= dphi;
                    //                    INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                    E = APtrack(cavity, cavity_start, phi);
                    //                    if(E > Emax) {
                        //                        INFOMSG(E - Emax << " MeV: accepted" << endl;);
                    //  } else {
                        //                        INFOMSG(E - Emax << " MeV: rejected" << " E= " << E << " Emax= " << Emax << endl;);
                    // }
                } while(E > Emax);

                if(j == 0) {
                    phi = Phiini;
                    E = Emax;
                    j = -1;
                    do {
                        j ++;
                        Emax = E;
                        Phiini = phi;
                        phi += dphi;
                        //                        INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                        E = APtrack(cavity, cavity_start, phi);
                        //                        if(E > Emax) {
                        //  INFOMSG(E - Emax << " MeV: accepted" << endl;);
                        //  } else {
                        //  INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        // }
                    } while(E > Emax);
                }
                for(int refinement_level = 0; refinement_level < numRefs; refinement_level ++) {
                    dphi /= 2.;
                    //                    INFOMSG("refinement level: " << refinement_level + 1 << ", dphi= " << dphi << endl;);
                    phi = Phiini - dphi;
                    // INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                    E = APtrack(cavity, cavity_start, phi);
                    if(E > Emax) {
                        //  INFOMSG(E - Emax << " MeV: accepted" << endl;);
                        Phiini = phi;
                        Emax = E;
                    } else {
                        // INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        phi = Phiini + dphi;
                        // INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                        E = APtrack(cavity, cavity_start, phi);
                        if(E > Emax) {
                            //  INFOMSG(E - Emax << " MeV: accepted" << endl;);
                            Phiini = phi;
                            Emax = E;
                        }
                        //else {
                        //  INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        //}
                    }
                }
                Phimax = Phiini;
                phi = Phimax + orig_phi;

                INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
            } else {
                msg << "Tracking with phi= " << orig_phi << " result E= " << E << endl;
                phi = orig_phi;
                Emax = E;
                E = APtrack(cavity, cavity_start, phi - Physics::pi);
                msg << "Tracking with phi= " << orig_phi - Physics::pi << " result E= " << E << endl;
                if (E > Emax) {
                    phi = orig_phi - Physics::pi;
                    Phimax = phi;
                    Emax = E;
                }
            }

            if(cavity->getType() == "TravelingWave") {
                static_cast<TravelingWave *>(cavity)->updatePhasem(phi);
            } else {
                static_cast<RFCavity *>(cavity)->updatePhasem(phi);
            }

            PhiAstra = (Phimax * RADDEG) + 90.0;
            PhiAstra -= floor(PhiAstra / 360.) * 360.;

            msg << cavity->getName() << "_phi= "  << Phimax << " rad / "
            << Phimax *RADDEG <<  " deg, AstraPhi= " << PhiAstra << " deg,\n"
            << "E= " << Emax << " (MeV), " << "phi_nom= " << orig_phi *RADDEG << endl;

            OpalData::getInstance()->setMaxPhase(cavity->getName(), Phimax);
        }

        // --- doOneStep(pusher);
        // --- timeIntegration1(pusher);
        IpplTimings::startTimer(timeIntegrationTimer1_m);
        double t = itsBunch->getT();

        // increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
        // before they are allocated when increasing the timestep in the gun.

        switchElements(10.0);
        itsOpalBeamline_m.resetStatus();


        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->getdT());
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->getdT());
        }
        IpplTimings::stopTimer(timeIntegrationTimer1_m);
        // --- timeIntegration1(pusher);

        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);

        // --- computeExternalFields();
        IpplTimings::startTimer(timeFieldEvaluation_m);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            long ls = itsBunch->LastSection[i];
            itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
            if(ls != itsBunch->LastSection[i]) {
                if(!itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls)) {
                    itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
                }
                itsBunch->LastSection[i] = ls;
            }

            Vector_t externalE, externalB;
            const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls, itsBunch->getT() + itsBunch->dt[i] / 2., externalE, externalB);
            if (rtv)
                ;
            itsBunch->Ef[i] += externalE;
            itsBunch->Bf[i] += externalB;
        }
        IpplTimings::stopTimer(timeFieldEvaluation_m);
        // --- computeExternalFields

        // --- timeIntegration2(pusher);
        IpplTimings::startTimer(timeIntegrationTimer2_m);
        double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;

        // kickParticlesAutophase(pusher);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            pusher.kick(itsBunch->R[i], itsBunch->P[i], itsBunch->Ef[i], itsBunch->Bf[i], itsBunch->dt[i]);
        }

        itsBunch->calcBeamParametersLight();

        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->getdT());
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->getdT());
        }
        t += itsBunch->getdT();
        itsBunch->setT(t);
        IpplTimings::stopTimer(timeIntegrationTimer2_m);
        // --- timeIntegration2(pusher);
        // --- doOneStep(pusher);

	double sposRef	= 0.0;

	if(!(step % 5000))	{

         if (itsBunch->getLocalNum() != 0)
           sposRef = itsBunch->R[0](2);

         reduce(sposRef,sposRef,OpAddAssign());

	 if(sposRef > zStop)
	   localTrackSteps_m = floor(step / dtfraction);


	  INFOMSG("step = " << step << ", spos = " << sposRef << " [m], t= " << itsBunch->getT() << " [s], "
		  << "E= " << itsBunch->get_meanEnergy() << " [MeV] " << endl);
	}
    }

    localTrackSteps_m = maxStepsSave;
    scaleFactor_m = scaleFactorSave;
    itsBunch->setT(tSave);
    itsBunch->setdT(dTSave);
}

double ParallelTTracker::APtrack(Component *cavity, double cavity_start_pos, const double &phi) const {
    double beta = std::max(sqrt(1. - 1 / (itsBunch->P[0](2) * itsBunch->P[0](2) + 1.)), 0.0001);
    double tErr  = (cavity_start_pos - itsBunch->R[0](2)) / (Physics::c * beta);

    //    INFOMSG("beta = " << beta << " tErr = " << tErr << endl;);

    double finalMomentum = 0.0;
    if(cavity->getType() == "TravelingWave") {
        TravelingWave *tws = static_cast<TravelingWave *>(cavity);
        tws->updatePhasem(phi);
        std::pair<double, double> pe = tws->trackOnAxisParticle(itsBunch->P[0](2),
                                                                itsBunch->getT() + tErr,
                                                                itsBunch->dt[0],
                                                                itsBunch->getQ(),
                                                                itsBunch->getM() * 1e-6);
        finalMomentum = pe.first;
    } else {
        RFCavity *rfc = static_cast<RFCavity *>(cavity);
        rfc->updatePhasem(phi);

        std::pair<double, double> pe = rfc->trackOnAxisParticle(itsBunch->P[0](2),
                                                                itsBunch->getT() + tErr,
                                                                itsBunch->dt[0],
                                                                itsBunch->getQ(),
                                                                itsBunch->getM() * 1e-6);
        finalMomentum = pe.first;
    }
    double finalGamma = sqrt(1.0 + finalMomentum * finalMomentum);
    double finalKineticEnergy = (finalGamma - 1.0) * itsBunch->getM() * 1e-6;
    return finalKineticEnergy;
}

#ifdef HAVE_AMR_SOLVER
void ParallelTTracker::Tracker_AMR()
{
    Inform msg("ParallelTTracker ");
    const Vector_t vscaleFactor_m = Vector_t(scaleFactor_m);
    BorisPusher pusher(itsReference);
    secondaryFlg_m = false;
    dtTrack_m = itsBunch->getdT();

    // upper limit of particle number when we do field emission and secondary emission
    // simulation. Could be reset to another value in input file with MAXPARTSNUM.
    maxNparts_m = 100000000;
    nEmissionMode_m = true;

    prepareSections();

    // do autophasing before tracking without a global phase shift!
    doAutoPhasing();

    numParticlesInSimulation_m = itsBunch->getTotalNum();

    OPALTimer::Timer myt1;

    setTime();

    double t = itsBunch->getT();

    unsigned long long step = itsBunch->getLocalTrackStep();

    msg << "Track start at: " << myt1.time() << ", t= " << t << "; zstop at: " << zStop_m << " [m]" << endl;

    gunSubTimeSteps_m = 10;
    prepareEmission();

    doSchottyRenormalization();

    msg << "Executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s];\n"
    << "max integration steps " << localTrackSteps_m << ", next step= " << step << endl;
    msg << "Using default (Boris-Buneman) integrator" << endl;

    // itsBeamline_m.accept(*this);
    // itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(msg);

    setupSUV();

    // increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
    // before they are allocated when increasing the timestep in the gun.
    switchElements(10.0);

    initializeBoundaryGeometry();

    setOptionalVariables();

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq_m = 1000000;

    wakeStatus_m = false;
    surfaceStatus_m = false;

    // reset E and B to Vector_t(0.0) for every step
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    for(; step < localTrackSteps_m; ++step)
    {
        bends_m = 0;
        numberOfFieldEmittedParticles_m = 0;

        itsOpalBeamline_m.resetStatus();

        // we dump later, after one step.
        // dumpStats(step, true, true);

        Real stop_time = -1.;

        std::cout << "                " << std::endl;
        std::cout << " ************** " << std::endl;
        std::cout << " DOING STEP ... " << step << std::endl;
        std::cout << " ************** " << std::endl;
        std::cout << "                " << std::endl;

        Real dt_from_amr = amrptr->coarseTimeStepDt(stop_time);

        std::cout << "                " << std::endl;
        std::cout << " ************** " << std::endl;
        std::cout << " COMPLETED STEP ... " << step << " WITH DT = " << dt_from_amr << std::endl;
        std::cout << " ************** " << std::endl;
        std::cout << "                " << std::endl;

        t += dt_from_amr;
        itsBunch->setT(t);

        bool const psDump = step % Options::psDumpFreq == 0;
        bool const statDump = step % Options::statDumpFreq == 0;
        dumpStats(step, psDump, statDump);

        if(hasEndOfLineReached()) break;

        double margin = 0.1;
        switchElements(margin);

        itsBunch->incTrackSteps();

        // These routines return the particle data for all of the particles and on all of the processes

        Array<int> particle_ids;
        amrptr->GetParticleIDs(particle_ids);

        Array<int> particle_cpu;
        amrptr->GetParticleCPU(particle_cpu);

        Array<Real> locs;
        amrptr->GetParticleLocations(locs);

        // Here we assume that we have stored, Q, V, ... in the particle data in TrackRun.cpp
        int start_comp = 1;
        int   num_comp = 3;
	Array<Real> Qs;
        Array<Real> vels;
        Array<Real> Evec;

        amrptr->GetParticleData(Qs,0,1);
        amrptr->GetParticleData(vels,start_comp,num_comp);
    	//amrptr->GetParticleData(vels,start_comp,6);
        amrptr->GetParticleData(Evec,4,num_comp);

        std::cout << "SIZE OF PARTICLE IDs "  << particle_ids.size() << std::endl;
        std::cout << "SIZE OF PARTICLE CPU "  << particle_cpu.size() << std::endl;
        std::cout << "SIZE OF PARTICLE LOCS " << locs.size() << std::endl;
        std::cout << "SIZE OF PARTICLE VELS " << vels.size() << std::endl;
        std::cout << "SIZE OF PARTICLE EFIELD " << Evec.size() << std::endl;


        int num_particles_total = particle_ids.size();

	double gamma=itsReference.getGamma();
	std:: cout << " GAMMA" << gamma << std::endl;

        Vector_t rmin;
        Vector_t rmax;
        itsBunch->get_bounds(rmin, rmax);

        FVector<double,6> six_vect;

	for (int i = 0; i < num_particles_total; i++)
        {
             if (i < 3 ) std::cout << "PARTICLE ID " << particle_ids[i] << "\n"
				    << Qs[i] << "\n"
                                    << "  " << locs[3*i  ] << " " << vels[3*i  ] << "\n"
		 		    << "  " << locs[3*i+1] << " " << vels[3*i+1] << "\n"
		 		    << "  " << locs[3*i+2] << " " << vels[3*i+2] << "\n"
#if 0
		 		    << "  " << locs[3*i]   << " " << vels[3*i+3] << "\n"
		 		    << "  " << locs[3*i+1] << " " << vels[3*i+4] << "\n"
		 		    << "  " << locs[3*i+2] << " " << vels[3*i+5] << "\n"
#endif
                                    << "  " << locs[3*i  ] << " " << Evec[3*i  ] << "\n"
		 		    << "  " << locs[3*i+1] << " " << Evec[3*i+1] << "\n"
		 		    << "  " << locs[3*i+2] << " " << Evec[3*i+2] << "\n"
		 		    << std::endl;
             if (particle_cpu[i] == Ippl::myNode())
             {
                 six_vect[0] = locs[3*i  ];
                 six_vect[1] = vels[3*i  ] * gamma / Physics::c;
                 six_vect[2] = locs[3*i+1];
                 six_vect[3] = vels[3*i+1] * gamma / Physics::c;
                 six_vect[4] = locs[3*i+2];
                 six_vect[5] = vels[3*i+2] * gamma / Physics::c;

                 // We subtract one from the particle_id because we added one to it when we
                 //    passed the particle into the AMR stuff.
                 // std::cout << "ON NODE ID " << Ippl::myNode() << " ADDING PARTICLE "
                 //          << particle_ids[i] << " WITH X,Y,Z "
                 //          << locs[3*i] << " " << locs[3*i+1] << " " << locs[3*i+2] << std::endl;
                 itsBunch->set_part(six_vect, particle_ids[i]-1);
		 for (int k=0; k<3; k++)
			itsBunch->Ef[particle_ids[i]-1](k) = Evec[3*i+k];
             }
	}
    }

    Vector_t rmin;
    Vector_t rmax;
    itsBunch->get_bounds(rmin, rmax);

    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
        numParticlesInSimulation_m = itsBunch->getTotalNum();
    }

    writePhaseSpace((step + 1), itsBunch->get_sPos(), true, true);
    msg << "Dump phase space of last step" << endl;
    OPALTimer::Timer myt3;
    itsOpalBeamline_m.switchElementsOff();
    msg << "done executing ParallelTTracker at " << myt3.time() << endl;
}
#endif

void ParallelTTracker::Tracker_Default() {
    Inform msg("ParallelTTracker ");
    const Vector_t vscaleFactor_m = Vector_t(scaleFactor_m);
    BorisPusher pusher(itsReference);
    secondaryFlg_m = false;
    dtTrack_m = itsBunch->getdT();

    double testPhi = 22.0 / 180.0 * Physics::pi;
    double testTheta = 112 / 180.0 * Physics::pi;
    Vector_t testA(1.0, 0,0);
    Vector_t testB(cos(testPhi)*sin(testTheta), sin(testPhi)*sin(testTheta), cos(testTheta));
    Quaternion testQ = getQuaternion(testA, testB);

    Vector_t testC = testQ.rotate(testA);
    msg << __FILE__ << ": " << __LINE__ << "\n"
        << testC << "\n"
        << testB << endl;

    // msg << testB << "\n"
    //     << testC << endl;
    // upper limit of particle number when we do field emission and secondary emission
    // simulation. Could be reset to another value in input file with MAXPARTSNUM.
    maxNparts_m = 100000000;
    nEmissionMode_m = true;

    prepareSections();

    // do autophasing before tracking without a global phase shift!
    doAutoPhasing();

    numParticlesInSimulation_m = itsBunch->getTotalNum();

    OPALTimer::Timer myt1;

    setTime();

    double t = itsBunch->getT();

    unsigned long long step = itsBunch->getLocalTrackStep();

    msg << "Track start at: " << myt1.time() << ", t= " << t << "; zstop at: " << zStop_m << " [m]" << endl;

    gunSubTimeSteps_m = 10;
    prepareEmission();

    doSchottyRenormalization();

    msg << "Executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s];\n"
    << "max integration steps " << localTrackSteps_m << ", next step= " << step << endl;
    msg << "Using default (Boris-Buneman) integrator" << endl;

    // itsBeamline_m.accept(*this);
    // itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(msg);

    setupSUV();

    // increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
    // before they are allocated when increasing the timestep in the gun.
    switchElements(10.0);

    initializeBoundaryGeometry();

    setOptionalVariables();

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq_m = 1000000;

    wakeStatus_m = false;
    surfaceStatus_m = false;

    for(; step < localTrackSteps_m; ++step) {
        bends_m = 0;
        numberOfFieldEmittedParticles_m = 0;

        itsOpalBeamline_m.resetStatus();

        // we dump later, after one step.
        // dumpStats(step, true, true);


        timeIntegration1(pusher);
        timeIntegration1_bgf(pusher);

        itsBunch->calcBeamParameters();

        // reset E and B to Vector_t(0.0) for every step
        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);

        if(step % repartFreq_m == 0 && step != 0) {
            doBinaryRepartition();
        }

        computeSpaceChargeFields();

        selectDT();
        emitParticles(step);
        selectDT();

        computeExternalFields();

        timeIntegration2(pusher);
        timeIntegration2_bgf(pusher);

        bgf_main_collision_test();

        //t after a full global timestep with dT "synchronization point" for simulation time
        t += itsBunch->getdT();
        itsBunch->setT(t);

        bool const psDump = step % Options::psDumpFreq == 0;
        bool const statDump = step % Options::statDumpFreq == 0;
        dumpStats(step, psDump, statDump);

        if(hasEndOfLineReached()) break;

        double margin = 0.1;
        switchElements(margin);

        itsBunch->incTrackSteps();

    }

    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
        numParticlesInSimulation_m = itsBunch->getTotalNum();
    }
    writePhaseSpace((step + 1), itsBunch->get_sPos(), true, true);
    msg << "Dump phase space of last step" << endl;
    OPALTimer::Timer myt3;
    itsOpalBeamline_m.switchElementsOff();
    msg << "done executing ParallelTTracker at " << myt3.time() << endl;
}

/**
 *  COMPONENTS
 */

double ParallelTTracker::getGlobalPhaseShift() {

    if(Options::autoPhase > 0) {
        double gPhaseSave = OpalData::getInstance()->getGlobalPhaseShift();
        OpalData::getInstance()->setGlobalPhaseShift(0.0);
        return gPhaseSave;
    } else
        return 0;

}

void ParallelTTracker::handleOverlappingMonitors() {
    // make sure that no monitor has overlap with two tracks
    Inform msg("ParallelTTracker ");
    FieldList monitors = itsOpalBeamline_m.getElementByType("Monitor");
    for(FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
        double zbegin, zend;
        it->getElement()->getDimensions(zbegin, zend);
        if(zbegin < zStop_m && zend >= zStop_m) {
            msg << "\033[0;31m"
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
            << "% Removing '" << it->getElement()->getName() << "' since it resides in two tracks.   %\n"
            << "% Please adjust zstop or place your monitor at a different position to prevent this. %\n "
            << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
            << "\033[0m"
            << endl;
            static_cast<Monitor *>(it->getElement())->moveBy(-zend - 0.001);
            itsOpalBeamline_m.removeElement(it->getElement()->getName());
        }
    }
}

void ParallelTTracker::prepareSections() {

    itsBeamline_m.accept(*this);
    handleOverlappingMonitors();
    itsOpalBeamline_m.prepareSections();

    cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType("TravelingWave");
    cavities_m.merge(travelingwaves, OpalField::SortAsc);
}


void ParallelTTracker::doAutoPhasing() {

    if(Options::autoPhase == 0) return;
    if(OpalData::getInstance()->inRestartRun()) {
        itsDataSink_m->retriveCavityInformation(OpalData::getInstance()->getInputFn());

        for(std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++)
            updateRFElement((*it).first, (*it).second);
    } else {
        if(OpalData::getInstance()->hasBunchAllocated()) {
            // we are in a followup track and the phase information is
            // already stored in the OPAL dictionary.
            for(std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
                updateRFElement((*it).first, (*it).second);
                INFOMSG("In follow-up track use saved phases for -> name: " << (*it).first << " phi= " << (*it).second << " (rad)" << endl);
            }
        } else {
            int tag = 101;
            int Parent = 0;

            itsBunch->stash();
            double zStop = itsOpalBeamline_m.calcBeamlineLenght();
            if(Ippl::myNode() == 0) {
                itsBunch->create(1);
                itsBunch->R[0] = Vector_t(0.0);
                Vector_t stashedP = itsBunch->getStashIniP();
                itsBunch->P[0] = Vector_t(0.0, 0.0, std::max(stashedP(2), 1e-6));
                itsBunch->Bin[0] = 0;
                itsBunch->Q[0] = itsBunch->getChargePerParticle();
                itsBunch->PType[0] = 0;
                itsBunch->LastSection[0] = 0;
                itsBunch->update();
                executeAutoPhase(Options::autoPhase, zStop);

                // now send all max phases and names of the cavities to
                // all the other nodes for updating.
                Message *mess = new Message();
                putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

                for(std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
                    putMessage(*mess, (*it).first);
                    putMessage(*mess, (*it).second);
                }
                Ippl::Comm->broadcast_all(mess, tag);

                itsBunch->destroy(1, 0);
            } else {
                // receive max phases and names and update the structure
                itsBunch->update();
                executeAutoPhase(Options::autoPhase, zStop);
                int nData = 0;
                Message *mess = Ippl::Comm->receive_block(Parent, tag);
                getMessage(*mess, nData);
                for(int i = 0; i < nData; i++) {
                    std::string elName;
                    double maxPhi;
                    getMessage(*mess, elName);
                    getMessage(*mess, maxPhi);
                    updateRFElement(elName, maxPhi);
                    OpalData::getInstance()->setMaxPhase(elName, maxPhi);
                }
            }
            itsBunch->update();
            itsBunch->pop();
            itsDataSink_m->storeCavityInformation();
        }
    }
    updateAllRFElements(OpalData::getInstance()->getGlobalPhaseShift());
    INFOMSG("finished autophasing" << endl);
}

void ParallelTTracker::bgf_main_collision_test() {

    Inform msg("ParallelTTracker ");

    if(!bgf_m) return;

    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();
    /**
     Here we check if a particles is
     outside the domain, flag it for
     deletion and create secondaries
     */
    if(secondaryFlg_m == 1) {
        /*
         entry for Furman-Pivi's secondary emission model
         itsBunch->getLocalNum() will change immediately, so we
         need Inc_num to record the local particle number before
         secondary emission, otherwise will be recursive generate
         secondaries and cause problem.
         */
        size_t Inc_num = itsBunch->getLocalNum();

        double dtime = 0.5 * itsBunch->getdT();

        double seyNum = 0;

        for(size_t i = 0; i < Inc_num; i++) {

            if(itsBunch->PType[i] == 3)
                // secondaries generated in last step will be set to be old
                // secondaries.

                itsBunch->PType[i] = 2;

            if(itsBunch->TriID[i] == 0) {
                /*
                 for primary bunch, primary dark current particles,
                 old secondaries in previous time steps and newly
                 generated secondaries which have no collision with
                 boundary in both first and second half step, do main
                 collision test and emit the secondaries.
                 */
                Vector_t intecoords = outr;
                int triId = 0;
                int res = bgf_m->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId);
                if(res == 0) {
                    res += bgf_m->doBGphysics(intecoords, triId, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum);
                }
                if(res >= 0) {
                    itsBunch->Bin[i] = -1;
                    Nimpact_m++;
                    SeyNum_m += seyNum;
                }
            } else {
                /*
                 Particles which collide the boundary in previous
                 two tests will not do main collision test and directly
                 call secondary emission module according to their
                 energy and momentum before collision. Attention, these
                 secondaries have not been kicked and are without new
                 momentum.
                 */

                int triId = itsBunch->TriID[i];

                int res = bgf_m->doBGphysics(itsBunch->R[i], triId, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum);
                if(res >= 0) {
                    itsBunch->Bin[i] = -1;
                    Nimpact_m++;
                    SeyNum_m += seyNum;
                }
            }
        }

        /*===========================
         Now we do fieldemission
         ============================== */
        if(itsBunch->getT() < surfaceEmissionStop_m) {
            numberOfFieldEmittedParticles_m += bgf_m->doFNemission(itsOpalBeamline_m, itsBunch, itsBunch->getT());
            itsBunch->boundp();
            numParticlesInSimulation_m = itsBunch->getTotalNum();
        } else
            msg << "* No field emission dT = " << itsBunch->getT() << endl;

    } else if(secondaryFlg_m != 0) {
        // entry for Vaughan's secondary emission model

        const int para_null = 0;// dummy parameter for overloading the Vaughan's version of BoundaryGeometry::doBGphysics();

        // itsBunch->getLocalNum() will change immediately, so we need Inc_num to record the
        // local particle number before secondary emission, otherwise will be recursive generate secondaries and cause problem.
        size_t Inc_num = itsBunch->getLocalNum();

        double dtime = 0.5 * itsBunch->getdT();

        double seyNum = 0;

        for(size_t i = 0; i < Inc_num; i++) {

            if(itsBunch->PType[i] == 3)
                itsBunch->PType[i] = 2;// secondaries generated in last step will be set to be old secondaries.
            // for primary bunch, primary dark current particles, old secondaries in previous time steps and newly generated
            // secondaries which have no collision with boundary in both first and second half step, do main collision test and emit the secondaries.
            if(itsBunch->TriID[i] == 0) {
                Vector_t intecoords = outr;
                int triId = 0;

                int res = bgf_m->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId);

                if(res == 0) {
                    res += bgf_m->doBGphysics(intecoords, triId, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum, para_null);
                }

                if(res >= 0) {
                    itsBunch->Bin[i] = -1;
                    Nimpact_m++;
                    SeyNum_m += seyNum;
                }
            } else {
                // Particles which collide the boundary in previous two tests will not do main collision test and directly call
                // secondary emission module according to their energy and momentum before collision. Attention, these secondaries have not
                // been kicked and are without new momentum.
                int triId = itsBunch->TriID[i];
                //assert(dot(itsBunch->P[i], bgf_m->TriNormals_m[triId]) < 0);
                int res = bgf_m->doBGphysics(itsBunch->R[i], triId, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum, para_null);

                if(res >= 0) {
                    itsBunch->Bin[i] = -1;
                    Nimpact_m++;
                    SeyNum_m += seyNum;
                }
            }
        }

        /*===========================
         Now we do fieldemission
         ============================== */
        if(itsBunch->getT() < surfaceEmissionStop_m) {
            numberOfFieldEmittedParticles_m += bgf_m->doFNemission(itsOpalBeamline_m, itsBunch, itsBunch->getT());
            itsBunch->boundp();
            numParticlesInSimulation_m = itsBunch->getTotalNum();
        } else
            msg << "* No field emission dT = " << itsBunch->getT() << endl;

    } else {// the case without secondary emission, i.e., secondaryFlg==0

        for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
            Vector_t intecoords = outr;
            if(itsBunch->TriID[i] == 0) { // Particles which do not collide the boundary in collision test after kick
                int triId = 0;
                int res = bgf_m->PartInside(itsBunch->R[i], itsBunch->P[i], itsBunch->getdT(), itsBunch->PType[i], itsBunch->Q[i], intecoords, triId);
                if(res == 0) {
                    res += bgf_m->doBGphysics(intecoords, triId);
                }
                if(res >= 0) {
                    itsBunch->Bin[i] = -1;
                    Nimpact_m++;
                }
            } else {// Particles which collide the boundary in collision test after kick will not do main collision test and directly call doBGphysics function.
                int triId = itsBunch->TriID[i];

                //assert(dot(itsBunch->P[i], bgf_m->TriNormals_m[triId]) < 0);
                int res = bgf_m->doBGphysics(intecoords, triId);

                if(res >= 0) {
		    itsBunch->Bin[i] = -1;
		    Nimpact_m++;
                }
            }
        }

        /*========================
	  Now we do fieldemission
	  =========================== */
        if(itsBunch->getT() < surfaceEmissionStop_m)
	  numberOfFieldEmittedParticles_m += bgf_m->doFNemission(itsOpalBeamline_m, itsBunch, itsBunch->getT());
        else
	  msg << "* No field emission dT = " << itsBunch->getT() << endl;
        /*  if (itsBunch->getTotalNum()!= 0) {
	    itsBunch->boundp();
	    msg<<"After boundp"<<endl;
	    }
	    numParticlesInSimulation_m = itsBunch->getTotalNum();
	*/
    }

    for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
      if (bgf_m->isOutsideApperture(itsBunch->R[i])) {
	itsBunch->Bin[i] = -1;
      }
    }

    itsBunch->boundp();

    numParticlesInSimulation_m = itsBunch->getTotalNum();
}

void ParallelTTracker::timeIntegration1(BorisPusher & pusher) {
    IpplTimings::startTimer(timeIntegrationTimer1_m);

    if(bgf_m != NULL && secondaryFlg_m > 0) return;

    itsBunch->switchToUnitlessPositions();
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        //scale each particle with c*dt
        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

        // update local coordinate system of particleInform &PartBunc
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
    }
    itsBunch->switchOffUnitlessPositions();

    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
    }

    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}

void ParallelTTracker::timeIntegration1_bgf(BorisPusher & pusher) {
    if(bgf_m == NULL || secondaryFlg_m == 0) return;

    IpplTimings::startTimer(timeIntegrationTimer1_m);

    /// We do collision test for newly generated secondaries before integration in the first half step of each time step.
    /// This is because only secondary emission model yield non zero inital momenta. The initial momenta of field emitted particles are zero.
    //  If hit, we set itsBunch->R[i] to intersection points, else we do normal integration.
    Nimpact_m = 0; // Initial parallel plate benchmark variable.
    SeyNum_m = 0; // Initial parallel plate benchmark variable.

    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();
    double dt = itsBunch->getdT();
    double bgf_scaleFactor = dt * Physics::c;
    Vector_t bgf_vscaleFactor = Vector_t(bgf_scaleFactor);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        bool particleHitBoundary = false;
        Vector_t intecoords = outr;
        int triId = 0;

        if(itsBunch->PType[i] == 3) { // only test newly generated secondaries
            particleHitBoundary = bgf_m->PartInside(itsBunch->R[i], itsBunch->P[i], 0.5 * itsBunch->dt[i], itsBunch->PType[i], itsBunch->Q[i], intecoords, triId) == 0;
        }


        if(particleHitBoundary) {// if hit, set particle position to intersection points coordinates and scale the position;
            // no scaling required
            itsBunch->R[i] = intecoords/bgf_vscaleFactor;
            itsBunch->TriID[i] = triId;
        } else {
            itsBunch->R[i] /= bgf_vscaleFactor;
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
        }
        // FIXME, is the local update necessary here?
        // update local coordinate system for particle
        itsBunch->X[i] /= bgf_vscaleFactor;
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());

        itsBunch->R[i] *= bgf_vscaleFactor;
        itsBunch->X[i] *= bgf_vscaleFactor;


    }


    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
    }
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}

void ParallelTTracker::timeIntegration2(BorisPusher & pusher) {
    if(bgf_m) return;
    IpplTimings::startTimer(timeIntegrationTimer2_m);

    /*
     transport and emit particles
     that passed the cathode in the first
     half-step or that would pass it in the
     second half-step.

     to make IPPL and the field solver happy
     make sure that at least 10 particles are emitted

     also remember that node 0 has
     all the particles to be emitted

     this has to be done *after* the calculation of the
     space charges! thereby we neglect space charge effects
     in the very first step of a new-born particle.

     */

    // push the reference particle by a half step
    double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;

    kickParticles(pusher);

    handleBends();

    //switchElements();

    itsBunch->switchToUnitlessPositions(true);
    // start normal particle loop part 2 for simulation without boundary geometry.
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} + \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f]
         * \code
         * R[i] += 0.5 * P[i] * recpgamma;
         * \endcode
         */
        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
        //and scale back to dimensions

        // update local coordinate system
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
        //reset time step if particle was emitted in the first half-step
        //the particle is now in sync with the simulation timestep

    }

    itsBunch->switchOffUnitlessPositions(true);

    fill(itsBunch->dt.begin(), itsBunch->dt.end(), itsBunch->getdT());

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTTracker::timeIntegration2_bgf(BorisPusher & pusher) {

    if(!bgf_m) return;

    /// After kick, we do collision test before integration in second half step with new momentum, if hit, then move collision particles to the position where collision occurs.

    IpplTimings::startTimer(timeIntegrationTimer2_m);

    // push the reference particle by a half step
    double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;


    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

	}
    kickParticles(pusher, 0);
    handleBends();
    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();

    double dtime = 0.5 * itsBunch->getdT();

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        bool particleHitBoundary = false;
        Vector_t intecoords = outr;
        int triId = 0;
        itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
        if(itsBunch->TriID[i] == 0) { // test all particles except those already have collided the boundary in the first half step.
            Vector_t scale_factor(0.0);

            particleHitBoundary =  bgf_m->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId) == 0;

            if(particleHitBoundary) {
                itsBunch->R[i] = intecoords / Vector_t(Physics::c * itsBunch->dt[i]);
                itsBunch->TriID[i] = triId;
                scale_factor = vscaleFactor_m;
            } else {//if no collision do normal push in the second half-step
                itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
                pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            }
            itsBunch->X[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
        }
        itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i]);
        itsBunch->X[i] *= Vector_t(Physics::c * itsBunch->dt[i]);
    }

    fill(itsBunch->dt.begin(), itsBunch->dt.end(), itsBunch->getdT());

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTTracker::selectDT() {

    double dt = dtTrack_m;
    itsBunch->setdT(dt);
    scaleFactor_m = dt * Physics::c;
    vscaleFactor_m = Vector_t(scaleFactor_m);

    if(itsBunch->GetIfBeamEmitting()) {
        dt = itsBunch->GetEmissionDeltaT();
        itsBunch->setdT(dt);
        scaleFactor_m = dt * Physics::c;
        vscaleFactor_m = Vector_t(scaleFactor_m);
    }

}


void ParallelTTracker::emitParticles(long long step) {

    if(!itsBunch->WeHaveEnergyBins()) return;

    if(itsBunch->GetIfBeamEmitting()) {
        int ne = 0;
        itsBunch->switchToUnitlessPositions(true);

        Vector_t externalE = Vector_t(0.0);
        Vector_t externalB = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(Vector_t(0.0),
                                     Vector_t(0.0),
                                     itsBunch->getT() + 0.5 * itsBunch->getdT(),
                                     externalE,
                                     externalB);
        ne += itsBunch->EmitParticles(externalE(2));

        if(Options::schottkyCorrection)
            applySchottkyCorrection(*itsBunch, ne, itsBunch->getT(), rescale_coeff_m);

        reduce(ne, ne, OpAddAssign());
        numParticlesInSimulation_m += ne;
        itsBunch->switchOffUnitlessPositions(true);

        if (numParticlesInSimulation_m > minBinEmitted_m)
            itsBunch->boundp();
    }

    if(step > minStepforReBin_m) {
        itsBunch->Rebin();
        itsBunch->resetInterpolationCache(true);
	//	if (itsBunch->getTotalNum() < specifiedNPart_m) {
	//  WARNMSG("Rebinning, but less particles emitted than specifyed. Either increase MINSTEPFORREBIN or you have particle losses!"<< endl);
	//}

    }

}


void ParallelTTracker::computeSpaceChargeFields() {
    if(numParticlesInSimulation_m <= minBinEmitted_m) return;

    // itsBunch->switchToUnitlessPositions(true);
    // // FIXME! why do we have to do compute self fields with unitless positions?
    // itsBunch->boundp();

    if(itsBunch->WeHaveEnergyBins()) {
        itsBunch->calcGammas();
        itsBunch->resetInterpolationCache();
        ParticleAttrib<double> Q_back = itsBunch->Q;
        for(int binNumber = 0; binNumber < itsBunch->GetLastEmittedEnergyBin() &&
            binNumber < itsBunch->GetNumberOfEnergyBins(); ++binNumber) {

            itsBunch->setBinCharge(binNumber);
            itsBunch->computeSelfFields(binNumber);
            itsBunch->Q = Q_back;

        }

    } else {
        itsBunch->computeSelfFields();
    }

    // itsBunch->switchOffUnitlessPositions(true);
    // // FIXME! (see above)
    // itsBunch->boundp();
}


void ParallelTTracker::computeExternalFields() {
    IpplTimings::startTimer(timeFieldEvaluation_m);
    Inform msg("ParallelTTracker ");

    unsigned long hasWake = 0;
    unsigned long hasSurfacePhysics = 0;
    long wfSection = 0;
    long sphysSection = 0;

    globalEOL_m = true;
    bool emission_in_progress = itsBunch->GetIfBeamEmitting();
    if(numParticlesInSimulation_m == 0 && emission_in_progress)
        globalEOL_m = false;

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

        long ls = itsBunch->LastSection[i];
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
        if(ls != itsBunch->LastSection[i]) {
            if(!(itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls) ||
                 itsOpalBeamline_m.section_is_glued_to(ls, itsBunch->LastSection[i]))) {
                itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
            }
            itsBunch->LastSection[i] = ls;
        }

        Vector_t externalE, externalB;
        const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls, itsBunch->getT() + itsBunch->dt[i] / 2., externalE, externalB);

        globalEOL_m = globalEOL_m && (rtv & BEAMLINE_EOL);
        if((rtv & BEAMLINE_WAKE) && hasWake == 0) {
            wfSection = ls;
            hasWake = 1;
        }

        if((rtv & BEAMLINE_SURFACEPHYSICS) && hasSurfacePhysics == 0) {
            sphysSection = ls;
            hasSurfacePhysics = 1;
        }

        bends_m = bends_m || (rtv & BEAMLINE_BEND);

        // skip rest of the particle push if the
        // particle is out of bounds i.e. does not see
        // a E or B field
        if(rtv & BEAMLINE_OOB)
            itsBunch->Bin[i] = -1;

        itsBunch->Ef[i] += externalE;
        itsBunch->Bf[i] += externalB;

        // in case a particle is pushed behind the emission surface, delete the particle
        if(itsBunch->R[i](2) < 0)
            itsBunch->Bin[i] = -1;

    }

    IpplTimings::stopTimer(timeFieldEvaluation_m);

    reduce(hasWake, hasWake, OpAddAssign());
    reduce(hasSurfacePhysics, hasSurfacePhysics, OpAddAssign());
    reduce(bends_m, bends_m, OpAddAssign());

    if(hasWake > 0) {
        IpplTimings::startTimer(WakeFieldTimer_m);
        reduce(wfSection, wfSection, OpMaxAssign());
        WakeFunction *wf = itsOpalBeamline_m.getWakeFunction(wfSection);
        /*--------- Added by Xiaoying Pang 04/22/2014 ---------------
         * If the CSR is turned on for a dipole, save its pointer to the CSRWakeFunction
         * and reuse it in the following drift.*/
        // xpang: start
        const ElementBase* element = itsOpalBeamline_m.getWakeFunctionOwner(wfSection);
        if(dynamic_cast<CSRWakeFunction*>(wf))
        {
          if(dynamic_cast<const RBend*>(element) || dynamic_cast<const SBend*>(element))
            wakeFunction_m = wf;
          if(dynamic_cast<const Drift*>(element))
            wf = wakeFunction_m;
        }
        // xpang: end
        if(!wakeStatus_m) {
            msg << "============== START WAKE CALCULATION =============" << endl;
            const ElementBase* element = itsOpalBeamline_m.getWakeFunctionOwner(wfSection);
            wf->initialize(element);
            wakeStatus_m = true;
        }
        if(wf == NULL) {
            INFOMSG("no wakefunction attached" << endl);
        } else {
            wf->apply(*itsBunch);
        }
        IpplTimings::stopTimer(WakeFieldTimer_m);

    } else if(wakeStatus_m) {
        msg << "=============== END WAKE CALCULATION ==============" << endl;
        wakeStatus_m = false;
    }

    if(hasSurfacePhysics > 0) {
        if(!surfaceStatus_m) {
            msg << "============== START SURFACE PHYSICS CALCULATION =============" << endl;
            surfaceStatus_m = true;
        }
        reduce(sphysSection, sphysSection, OpMaxAssign());
	if (sphys_m==NULL)
	  sphys_m = itsOpalBeamline_m.getSurfacePhysicsHandler(sphysSection);
	else {
	  /*
	   In case we have an other
	     handler, delete the first one and use the new one
	  */
	  SurfacePhysicsHandler *sphysNew = itsOpalBeamline_m.getSurfacePhysicsHandler(sphysSection);
	  if (sphysNew != sphys_m && sphysNew != NULL) {
	    //delete sphys_m; FixMe: should do this right!!!
	    sphys_m = sphysNew;
	  }
	}

        if(sphys_m == NULL) {
            INFOMSG("no surface physics attached" << endl);
        } else {
            sphys_m->apply(*itsBunch);
            if (itsBunch->getTotalNum()>10)
                sphys_m->print(msg);
        }
    } else if(surfaceStatus_m) {
        msg << "============== END SURFACE PHYSICS CALCULATION =============" << endl;
        surfaceStatus_m = false;
	/*
		 if (sphys_m) {
		 delete sphys_m;
		 sphys_m = NULL;
		 }
	*/
    }
    size_t ne = 0;
    bool globPartOutOfBounds = (min(itsBunch->Bin) < 0) && (itsBunch->getTotalNum() > 10);
    if(globPartOutOfBounds) {
      ne = itsBunch->boundp_destroyT();
      numParticlesInSimulation_m  = itsBunch->getTotalNum();
    }

    if (itsBunch->getTotalNum() > 10)
        itsBunch->update();

    if(ne > 0)
       msg << "* Deleted " << ne << " particles, "
	   << "remaining " << itsBunch->getTotalNum() << " particles" << endl;



}

void ParallelTTracker::handleBends() {
    if(bends_m == 0) {
        if(itsBunch->getTotalNum() > 0) {
            // none of the particles is in a bending element
            updateReferenceParticle();
        }
    } else {
        /*
         at least one of the elements bends the beam; until all
         particles have left the bending elements we track the
         reference particle as if it were a regular particle; from
         the moment when the reference particle has reached the
         bending field until it leaves it again we rotate the bunch
         about the position of the reference particle such that the
         momentum of the reference particle points in z direction;
         First update the momentum of the reference particle in zxy
         coordinate system, then update its position
         */

        updateSpaceOrientation(false);
        double recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

        RefPartP_zxy_m = dot(space_orientation_m, RefPartP_suv_m);
        RefPartR_zxy_m += RefPartP_zxy_m * recpgamma * scaleFactor_m / 2.;

        RefPartP_suv_m = Vector_t(0.0, 0.0, sqrt(dot(RefPartP_suv_m, RefPartP_suv_m)));
        RefPartR_suv_m += RefPartP_suv_m * recpgamma / 2. * Physics::c * itsBunch->getdT();
    }

    itsBunch->RefPart_R = RefPartR_zxy_m;
    itsBunch->RefPart_P = RefPartP_zxy_m;
}

void ParallelTTracker::switchElements(double scaleMargin) {
    // calculate the dimensions of the bunch and add a small margin to them; then decide which elements have to be triggered
    // when an element is triggered memory is allocated and the field map is read in
    Vector_t rmin(0.0), rmax(0.0);
    itsBunch->get_bounds(rmin, rmax);

    //FIXME: necessary
    RefPartP_suv_m = itsBunch->get_pmean();
    double recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

    // trigger the elements
    double margin = scaleMargin * RefPartP_suv_m(2) * recpgamma;
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements((rmin(2) - margin), (rmax(2) + margin));
}

void ParallelTTracker::doBinaryRepartition() {
    size_t particles_or_bins = std::max(minBinEmitted_m, size_t(1000));
    if(itsBunch->hasFieldSolver() && numParticlesInSimulation_m > particles_or_bins) {

        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition because of repartFreq_m" << endl);
        INFOMSG("*****************************************************************" << endl);
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition done" << endl);
        INFOMSG("*****************************************************************" << endl);
    }
}

void ParallelTTracker::dumpStats(long long step, bool psDump, bool statDump) {
    OPALTimer::Timer myt2;
    Inform msg("ParallelTTracker ");

    if(numParticlesInSimulation_m == 0) {
        msg << myt2.time() << " "
        << "Step " << setw(6) <<  itsBunch->getGlobalTrackStep() << "; "
        << "   -- no emission yet --     "
        << "t= "   << scientific << setprecision(3) << setw(10) << itsBunch->getT() << " [s]"
        << endl;
        return;
    }

    double sposRef = itsBunch->get_sPos();
    double sposPrint = sposRef;
    std::string sposUnit(" [m] ");
    double meanEnergy = itsBunch->get_meanEnergy();
    std::string meanEnergyUnit(" [MeV] ");

    if (sposRef < 1.0) {
        sposPrint = 1000.0*sposRef;
        sposUnit = std::string(" [mm] ");
    }

    if (meanEnergy < 1.0) {
        meanEnergy *= 1000.0;
        meanEnergyUnit = std::string(" [keV] ");
    }

    size_t totalParticles_f = numParticlesInSimulation_m;
    if(totalParticles_f <= minBinEmitted_m) {
        msg << myt2.time() << " "
        << "Step " << setw(6) << itsBunch->getGlobalTrackStep() << "; "
        << "only " << setw(4) << totalParticles_f << " particles emitted; "
        << "t= "   << scientific << setprecision(3) << setw(10) << itsBunch->getT() << " [s] "
        << "E="    << fixed      << setprecision(3) << setw(9) << meanEnergy << meanEnergyUnit
        << endl;
    } else if(std::isnan(sposRef) || std::isinf(sposRef)) {
        throw OpalException("ParallelTTracker::dumpStats()",
                            "there seems to be something wrong with the position of the bunch!");
    } else {
        msg << myt2.time() << " "
        << "Step " << setw(6) <<  itsBunch->getGlobalTrackStep() << " "
        << "at " << fixed      << setprecision(3) << setw(8) << sposPrint << sposUnit
        << "t= " << scientific << setprecision(3) << setw(10) << itsBunch->getT() << " [s] "
        << "E="  << fixed      << setprecision(3) << setw(9) << meanEnergy << meanEnergyUnit
        << endl;

        writePhaseSpace(step, sposRef, psDump, statDump);
    }

    if(bgf_m) {
        reduce(SeyNum_m, SeyNum_m, OpAddAssign());
        reduce(Nimpact_m, Nimpact_m, OpAddAssign());

        itsDataSink_m->writePartlossZASCII(*itsBunch, *bgf_m, std::string("data/Partloss-"));

        long long ustep = step;
        itsDataSink_m->writeImpactStatistics(*itsBunch, ustep, Nimpact_m, SeyNum_m, numberOfFieldEmittedParticles_m, nEmissionMode_m, std::string("data/PartStatistics"));

        if(((Options::surfDumpFreq) > 0) && ((step % Options::surfDumpFreq) == 0)) {
            itsDataSink_m->writeSurfaceInteraction(*itsBunch, ustep, *bgf_m, std::string("SurfaceInteraction"));
        }

        // If we are dealing with field emission and secondary emission, set upper
        // limit of particle number in simulation to prevent memory overflow.
        if(numParticlesInSimulation_m > maxNparts_m)
            localTrackSteps_m = step;

        // ada reset Nimpact_m, does not make sense to integrate this we obtain a rediculus large number !!
        Nimpact_m = 0;
    }

    if(sposRef > zStop_m)
        localTrackSteps_m = step;
}


void ParallelTTracker::setOptionalVariables() {
    Inform msg("ParallelTTracker ");

    minBinEmitted_m  = 10;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar)
        minBinEmitted_m = static_cast<size_t>(ar->getReal());
    msg << "MINBINEMITTED " << minBinEmitted_m << endl;

    minStepforReBin_m  = 200;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << "MINSTEPFORREBIN " << minStepforReBin_m << endl;

    surfaceEmissionStop_m  = 1000.0;
    RealVariable *cr = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("SURFACEEMISSIONSTOP"));
    if(cr)
        surfaceEmissionStop_m = static_cast<double>(cr->getReal());
    msg << "SURFACEEMISSIONSTOP after " << surfaceEmissionStop_m << " seconds" <<  endl;

    repartFreq_m = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep)
        repartFreq_m = static_cast<int>(rep->getReal());
    msg << "REPARTFREQ " << repartFreq_m << endl;

}


bool ParallelTTracker::hasEndOfLineReached() {
    reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
    return globalEOL_m;
}


void ParallelTTracker::doSchottyRenormalization() {
    Inform msg("ParallelTTracker ");
    double init_erg = itsBunch->getEkin();
    double tol_iter = 1e-5;
    rescale_coeff_m = 1 / init_erg / init_erg;

    if(Options::schottkyRennormalization > 0) {
        rescale_coeff_m = Options:: schottkyRennormalization;
        msg << "Set schottky scale coefficient to  " << rescale_coeff_m << endl;
    } else if(Options::schottkyCorrection) {
        while(true) {
            double real_charge = schottkyLoop(rescale_coeff_m);

            double total_charge = itsBunch->getTotalNum() * itsBunch->getChargePerParticle();
            msg << "Schottky scale coefficient " << rescale_coeff_m << ", actual emitted charge " << real_charge << " (Cb)" << endl;
            itsBunch->cleanUpParticles();
            itsBunch->setT(0);
            double scale_error = total_charge / real_charge - 1;
            // TODO : send scale_error to all nodes
            rescale_coeff_m *= (1.3 * scale_error + 1);
            if(fabs(scale_error) < tol_iter)
                break;
        }
        msg << "Schottky scan, final scale coefficient " << rescale_coeff_m << " ()" << endl;
    }
}

void ParallelTTracker::setupSUV() {

    if(mpacflg_m) return;

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
        itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l), itsOpalBeamline_m.getSectionStart(l));
    }

    if(!(itsBunch->WeHaveEnergyBins())) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
    }

    // Check if there are any particles in simulation. If there are,
    // as in a restart, use the usual function to calculate beam
    // parameters. If not, calculate beam parameters of the initial
    // beam distribution.
    if(numParticlesInSimulation_m == 0) {// fixme: maybe cause nonsense output if initialized momenta=0; Q: by Chuan.
        itsBunch->calcBeamParametersInitial();
    } else {
        itsBunch->calcBeamParameters();
    }

    RefPartP_suv_m = itsBunch->get_pmean();
    RefPartR_suv_m = itsBunch->get_rmean();
    RefPartP_zxy_m = RefPartP_suv_m;

    updateSpaceOrientation(false);
    RefPartP_suv_m = itsBunch->get_pmean();
}

void ParallelTTracker::setTime() {

    //if(Options::scan && OpalData::getInstance()->hasBunchAllocated())


    /*
     The track command set's T0 (or default 0.0)
     so we do not need this anymore

     if(!OpalData::getInstance()->hasBunchAllocated() &&
     !OpalData::getInstance()->inRestartRun())
     itsBunch->setT(0.0);
     */

    if(mpacflg_m) return;

    // set dt for all particles already in the simulation,
    // i.e. when doing a restarted simulation
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->dt[i] = itsBunch->getdT();
    }
}

void ParallelTTracker::prepareEmission() {
    Inform msg("ParallelTTracker ");

    if(mpacflg_m || !itsBunch->doEmission()) return;

    emissionSteps_m = itsBunch->GetNumberOfEmissionSteps();
    msg << "Do emission for " << itsBunch->getTEmission() << " [s] using "
          << itsBunch->GetNumberOfEnergyBins() << " energy bins " << endl
          << "Change dT from " <<  itsBunch->getdT() << " [s] to "
          <<  itsBunch->GetEmissionDeltaT() << " [s] during emission " << endl;;

}

void ParallelTTracker::initializeBoundaryGeometry() {

  Inform msg("ParallelTTracker ");

  /*
    For the moment, the Boundary geomentry must be attachedto the first element

  */

  bgf_m = itsOpalBeamline_m.getBoundaryGeometry(0);
  if (bgf_m) {
    Distribution *dist = NULL;
    Distribution *distrand = NULL;
    std::vector<std::string> distr_str = bgf_m->getDistributionArray();

    if(distr_str.size() == 0) {
      std::string distr = bgf_m->getDistribution();
      if(!distr.empty()) {
	msg << "* Find boundary geometry, start at: " << bgf_m->getS() << " (m) Distribution= " << bgf_m->getDistribution() << endl;
	dist = Distribution::find(bgf_m->getDistribution());
	msg << "* " << *dist << endl;
      } else {
	throw OpalException("ParallelTTracker::execute()",
			    "No distribution attached to BoundaryGeometry. Please check the input file... ...");
      }
    } else {
      msg << "************************************************************************************************* " << endl;
      msg <<  "* Find boundary geometry, start at: " << bgf_m->getS()  << " (m). " << endl;
      msg << "* Attached more than one distribution: " << endl;
      for(std::vector<std::string>::const_iterator dit = distr_str.begin(); dit != distr_str.end(); ++ dit) {
	Distribution *d = Distribution::find(*dit);
	msg << "* Distribution: " << *dit << " distribution type: " << d->GetTypeofDistribution() << endl;
	msg << "************************************************************************************************* " << endl;
	if(d->GetTypeofDistribution() == "SURFACEEMISSION") {
	  dist = d;
	  msg << *dist << endl;
	} else if(d->GetTypeofDistribution() == "SURFACERANDCREATE") {
	  distrand = d;
	  msg << *distrand << endl;
	  // here nbparts should be non zero as these particles will be the initialization of primary bunch.
	  size_t nbparts = distrand->GetNumberOfDarkCurrentParticles();
	  double darkinwardmargin = distrand->GetDarkCurrentParticlesInwardMargin();
	  double einitthreshold = distrand->GetEInitThreshold();
	  // make sure that the elements have already been switched on before initializing particles in position where the electric field > einitthreshold.
	  bgf_m->setEInitThreshold(einitthreshold);
	  if(!mpacflg_m) {
	    bgf_m->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
	    distrand->CreatePriPart(itsBunch, *bgf_m);
	    numParticlesInSimulation_m = itsBunch->getTotalNum();
	  } else {
	    /*
	      Multipacting flag set true. Generate primary particles.
	      Activate all elements (switch on the field map of elements in multipacting) in multipacting simulation
	    */

	    itsOpalBeamline_m.switchAllElements();
	    // it is possible to generate initial particles according to E field, since all elements switched on before we create particles.
	    bgf_m->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
	    // for Parallel Plate benchmark, Vw should be defined in input file and will be invoked by getVw method in createPriPart().
	    // for other multipacting simulation no need to define the Vw in SURFACERANDCREATE in input file.
	    distrand->CreatePriPart(itsBunch, *bgf_m);
	    numParticlesInSimulation_m = itsBunch->getTotalNum();
	    itsBunch->calcBeamParameters();
	    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
	      long &l = itsBunch->LastSection[i];
	      l = -1;
	      itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
	      itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l), itsOpalBeamline_m.getSectionStart(l));
	    }

	    // Check if there are any particles in simulation. If there are,
	    // as in a restart, use the usual function to calculate beam
	    // parameters. If not, calculate beam parameters of the initial
	    // beam distribution.

	    if(numParticlesInSimulation_m == 0) {
	      itsBunch->calcBeamParametersInitial();
	    } else {
	      itsBunch->calcBeamParameters();
	    }

	    //updateSpaceOrientation(false);
	    RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
	    RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();
	    msg << *itsBunch << endl;
	  }
	} else {
	  throw OpalException("ParallelTTracker::execute()",
			      "Unacceptable distribution type:\"" +
			      d->GetTypeofDistribution() + "\". Need to check the input file... ...");
	}
      }
    }

    /// this is still in BoundaryGeometry
    size_t nbparts = dist->GetNumberOfDarkCurrentParticles();
    double darkinwardmargin = dist->GetDarkCurrentParticlesInwardMargin();
    double workFunction = dist->GetWorkFunction();
    double fieldEnhancement = dist->GetFieldEnhancement();
    size_t maxfnemission = dist->GetMaxFNemissionPartPerTri();
    double fieldFNthreshold = dist->GetFieldFNThreshold();
    double parameterFNA = dist->GetFNParameterA();
    double parameterFNB = dist->GetFNParameterB();
    double parameterFNY = dist->GetFNParameterY();
    double parameterFNVYZe = dist->GetFNParameterVYZero();
    double parameterFNVYSe = dist->GetFNParameterVYSecond();

    secondaryFlg_m = dist->GetSecondaryEmissionFlag();
    nEmissionMode_m = dist->GetEmissionMode();
    bgf_m->setNEmissionMode(nEmissionMode_m);
    if(secondaryFlg_m) {
      if(secondaryFlg_m == 1) {
	int BoundaryMatType = dist->GetSurfMaterial();
	bgf_m->setBoundaryMatType(BoundaryMatType);
	if(Options::ppdebug) {
	  double vVThermal = dist->GetvVThermal();//return thermal velocity of Maxwellian distribution of secondaries for benchmark
	  bgf_m->setvVThermal(vVThermal);
	  double ppVw = dist->GetVw();
	  bgf_m->setVw(ppVw);
	} else {
	  bgf_m->setvVThermal(1.0);
	  bgf_m->setVw(1.0);
	}
      } else {
	/*
	  parameters for Vaughan's secondary model
	*/
	double vSeyZero = dist->GetvSeyZero();// return sey_0 in Vaughan's model
	double vEZero = dist->GetvEZero();// return the energy related to sey_0 in Vaughan's model
	double vSeyMax = dist->GetvSeyMax();// return sey max in Vaughan's model
	double vEmax = dist->GetvEmax();// return Emax in Vaughan's model
	double vKenergy = dist->GetvKenergy();// return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
	double vKtheta = dist->GetvKtheta();// return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
	double vVThermal = dist->GetvVThermal();// return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
	double ppVw = dist->GetVw();
	bgf_m->setVw(ppVw);
	bgf_m->setvSeyZero(vSeyZero);
	bgf_m->setvEZero(vEZero);
	bgf_m->setvSeyMax(vSeyMax);
	bgf_m->setvEmax(vEmax);
	bgf_m->setvKenergy(vKenergy);
	bgf_m->setvKtheta(vKtheta);
	bgf_m->setvVThermal(vVThermal);
      }
    }
    if(nbparts != 0) {
      //fixme: maybe need to be called in each time step for modeling creating darkcurrent in each time step
      bgf_m->createParticlesOnSurface(nbparts, darkinwardmargin, itsOpalBeamline_m, *itsBunch);
      dist->CreateBoundaryGeometry(*itsBunch, *bgf_m);
    }
    bgf_m->setWorkFunction(workFunction);
    bgf_m->setFieldEnhancement(fieldEnhancement);
    bgf_m->setMaxFN(maxfnemission);
    bgf_m->setFNTreshold(fieldFNthreshold);
    bgf_m->setFNParameterA(parameterFNA);
    bgf_m->setFNParameterB(parameterFNB);
    bgf_m->setFNParameterY(parameterFNY);
    bgf_m->setFNParameterVYZe(parameterFNVYZe);
    bgf_m->setFNParameterVYSe(parameterFNVYSe);
    numParticlesInSimulation_m = itsBunch->getTotalNum();
    if(numParticlesInSimulation_m > 0) {
      writePhaseSpace(0, 0, true, true); // dump the initial particles
    }
    itsDataSink_m->writeGeomToVtk(*bgf_m, std::string("data/testGeometry-00000.vtk"));
    //itsDataSink->writePartlossZASCII(*itsBunch, *bgf_m, std::string("vtk/PartlossZ-"));

    OpalData::getInstance()->setGlobalGeometry(bgf_m);

    RealVariable *maxnp = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MAXPARTSNUM"));
    if(maxnp) {
      maxNparts_m = static_cast<size_t>(maxnp->getReal());  // set upper limit of particle number in simulation
    }

    msg << "Boundary geometry initialized " << endl;
  }
}

void ParallelTTracker::execute() {
#ifdef HAVE_AMR_SOLVER
        Tracker_AMR();
#else
    if(timeIntegrator_m == 3) { // AMTS
        Tracker_AMTS();
    } else {
        Tracker_Default();
    }
#endif
}

void ParallelTTracker::Tracker_AMTS() {
    Inform msg("ParallelTTracker ");
    const Vector_t vscaleFactor_m = Vector_t(scaleFactor_m);
    dtTrack_m = itsBunch->getdT();

    // upper limit of particle number when we do field emission and secondary emission
    // simulation. Could be reset to another value in input file with MAXPARTSNUM.
    maxNparts_m = 100000000;

    prepareSections();

    // do autophasing before tracking without a global phase shift!
    doAutoPhasing();

    numParticlesInSimulation_m = itsBunch->getTotalNum();
    setTime();
    unsigned long long step = itsBunch->getLocalTrackStep();
    msg << "Track start at: " << OPALTimer::Timer().time() << ", t = " << itsBunch->getT() << "; zstop at: " << zStop_m << " [m]" << endl;
    msg << "Executing ParallelTTracker, next step = " << step << endl;
    msg << "Using AMTS (adaptive multiple-time-stepping) integrator" << endl;
    itsOpalBeamline_m.print(msg);
    setupSUV();

    itsOpalBeamline_m.switchAllElements();

    setOptionalVariables();

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq_m = 1000000;

    wakeStatus_m = false;
    surfaceStatus_m = false;

    // Count inner steps
    int totalInnerSteps = 0;

    itsBunch->boundp();
    itsBunch->calcBeamParameters();
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);
    computeSpaceChargeFields();
    if(itsBunch->WeHaveEnergyBins()) {
        itsBunch->Rebin();
        itsBunch->resetInterpolationCache(true);
    }

    // AMTS step size initialization
    double const dt_inner_target = itsBunch->getdT();
    msg << "AMTS initialization: dt_inner_target = " << dt_inner_target << endl;
    double dt_outer, deltaTau;
    if(itsBunch->deltaTau_m != -1.0) {
        // DTAU is set in the inputfile, calc initial outer time step from that
        deltaTau = itsBunch->deltaTau_m;
        dt_outer = calcG() * deltaTau;
    } else {
        // Otherwise use DTSCINIT
        dt_outer = itsBunch->dtScInit_m;
        deltaTau = dt_outer / calcG();
    }
    msg << "AMTS initialization: dt_outer = " << dt_outer << " deltaTau = " << deltaTau << endl;

    // AMTS calculation of stopping times
    double const tEnd = itsBunch->getT() + double(localTrackSteps_m - step) * dt_inner_target;
    double const psDumpInterval = double(Options::psDumpFreq) * dt_inner_target;
    double const statDumpInterval = double(Options::statDumpFreq) * dt_inner_target;
    double const repartInterval = double(repartFreq_m) * dt_inner_target;
    double const tTrackStart = itsBunch->getT() - double(step) * dt_inner_target; // we could be in a restarted simulation!
    double tNextPsDump = tTrackStart + psDumpInterval;
    while(tNextPsDump < itsBunch->getT()) tNextPsDump += psDumpInterval;
    double tNextStatDump = tTrackStart + statDumpInterval;
    while(tNextStatDump < itsBunch->getT()) tNextStatDump += statDumpInterval;
    double tDoNotRepartBefore = itsBunch->getT() + repartInterval;

    IpplTimings::startTimer(IpplTimings::getTimer("AMTS"));

    bool finished = false;
    for(; !finished; ++step) {
        itsOpalBeamline_m.resetStatus();

        // AMTS choose new timestep
        IpplTimings::startTimer(IpplTimings::getTimer("AMTS-TimestepSelection"));
        dt_outer = calcG() * deltaTau;
        double tAfterStep = itsBunch->getT() + dt_outer;
        double const tNextStop = std::min(std::min(tEnd, tNextPsDump), tNextStatDump);
        bool psDump = false, statDump = false;
        if(tAfterStep > tNextStop) {
            dt_outer = tNextStop - itsBunch->getT();
            tAfterStep = tNextStop;
        }
        double const eps = 1e-14; // To test approx. equality of times
        if(std::fabs(tAfterStep - tEnd) < eps) {
            finished = true;
        }
        if(std::fabs(tAfterStep - tNextPsDump) < eps) {
            psDump = true;
            tNextPsDump += psDumpInterval;
        }
        if(std::fabs(tAfterStep - tNextStatDump) < eps) {
            statDump = true;
            tNextStatDump += statDumpInterval;
        }
        msg << "AMTS: dt_outer = " << dt_outer;
        double numSubsteps = std::max(round(dt_outer / dt_inner_target), 1.0);
        msg << " numSubsteps = " << numSubsteps;
        double dt_inner = dt_outer / numSubsteps;
        msg << " dt_inner = " << dt_inner << endl;
        IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-TimestepSelection"));

        IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Kick"));
        if(itsBunch->hasFieldSolver()) {
            kick(0.5 * dt_outer);
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Kick"));

        for(int n = 0; n < numSubsteps; ++n) {
            bool const isFirstSubstep = (n == 0);
            bool const isLastSubstep = (n == numSubsteps - 1);
            borisExternalFields(dt_inner, isFirstSubstep, isLastSubstep);
            ++totalInnerSteps;
        }

        IpplTimings::startTimer(IpplTimings::getTimer("AMTS-SpaceCharge"));
        if(itsBunch->hasFieldSolver()) {
            itsBunch->boundp();
            itsBunch->Ef = Vector_t(0.0);
            itsBunch->Bf = Vector_t(0.0);
            if(itsBunch->getT() >= tDoNotRepartBefore) {
            	doBinaryRepartition();
            	tDoNotRepartBefore = itsBunch->getT() + repartInterval;
            }
            computeSpaceChargeFields();
            if(itsBunch->WeHaveEnergyBins()) {
                itsBunch->rebin();
                itsBunch->resetInterpolationCache(true);
            }
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-SpaceCharge"));

        IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Kick"));
        if(itsBunch->hasFieldSolver()) {
            kick(0.5 * dt_outer);
        }
        IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Kick"));

        IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Dump"));
        itsBunch->RefPart_R = RefPartR_zxy_m;
        itsBunch->RefPart_P = RefPartP_zxy_m;
        itsBunch->calcBeamParameters();
        dumpStats(step, psDump, statDump);
        IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Dump"));

        if(hasEndOfLineReached()) break;
        itsBunch->incTrackSteps();
    }

    IpplTimings::stopTimer(IpplTimings::getTimer("AMTS"));

    msg << "totalInnerSteps = " << totalInnerSteps << endl;

    itsBunch->boundp();
    numParticlesInSimulation_m = itsBunch->getTotalNum();
    writePhaseSpace((step + 1), itsBunch->get_sPos(), true, true);
    msg << "Dump phase space of last step" << endl;
    itsOpalBeamline_m.switchElementsOff();
    msg << "done executing ParallelTTracker at " << OPALTimer::Timer().time() << endl;
}

void ParallelTTracker::push(double h) {
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
        Vector_t const v = itsBunch->P[i] * Physics::c / gamma;
        itsBunch->R[i] += h * v;
        itsBunch->X[i] += h * TransformTo(v, itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i]));
    }

    // Push the reference particle
    double const gamma = sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    Vector_t v_zxy = RefPartP_zxy_m * Physics::c / gamma;
    RefPartR_zxy_m += h * v_zxy;

    itsBunch->setT(itsBunch->getT() + h);
}

void ParallelTTracker::kick(double h, bool avoidGammaCalc) {
    double const q = itsReference.getQ();
    double const M = itsReference.getM();
    double const h12Halfqc_M = 0.5 * h * q * Physics::c / M;
    double const h12Halfqcc_M = h12Halfqc_M * Physics::c;
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->P[i] += h12Halfqc_M * itsBunch->Ef[i];
    }
    if(avoidGammaCalc) {
        // Only calculate gamma if B-field nonzero...of course this check itself
        // gives some overhead, thus the caller can use the flag avoidGammaCalc
        Vector_t const zero = Vector_t(0.0);
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            if(itsBunch->Bf[i] != zero) {
                double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
                Vector_t const r = h12Halfqcc_M * itsBunch->Bf[i] / gamma;
                Vector_t const w = itsBunch->P[i] + cross(itsBunch->P[i], r);
                Vector_t const s = 2.0 / (1.0 + dot(r, r)) * r;
                itsBunch->P[i] += cross(w, s);
            }
        }
    } else {
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            double const gamma = sqrt(1.0 + dot(itsBunch->P[i], itsBunch->P[i]));
            Vector_t const r = h12Halfqcc_M * itsBunch->Bf[i] / gamma;
            Vector_t const w = itsBunch->P[i] + cross(itsBunch->P[i], r);
            Vector_t const s = 2.0 / (1.0 + dot(r, r)) * r;
            itsBunch->P[i] += cross(w, s);
        }
    }
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->P[i] += h12Halfqc_M * itsBunch->Ef[i];
    }
}

void ParallelTTracker::computeExternalFields_AMTS() {
    // We use a own function because the original computeExternalFields evaluates
    // the field at times itsBunch->getT() + itsBunch->dt[i] / 2, but we
    // need it evaluated at itsBunch->getT(). There is currently no need for individual time steps
    // in AMTS, because we don't support emission.
    // TODO: Unify these two functions
    Inform msg("ParallelTTracker ");
    bends_m = 0;
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        long ls = itsBunch->LastSection[i];
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
        if(ls != itsBunch->LastSection[i]) {
            if(!itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls)) {
                itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
            }
            itsBunch->LastSection[i] = ls;
        }
        Vector_t externalE, externalB;
        const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls, itsBunch->getT(), externalE, externalB);
        itsBunch->Ef[i] += externalE;
        itsBunch->Bf[i] += externalB;
        bends_m = bends_m || (rtv & BEAMLINE_BEND);
    }
    reduce(bends_m, bends_m, OpAddAssign());
}

void ParallelTTracker::borisExternalFields(double h, bool isFirstSubstep, bool isLastSubstep) {

    IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Push"));
    if(isFirstSubstep) {
        push(0.5 * h);
    }
    // Optimization: If this is not the first substep, our first half push was done in the last
    // substep.
    IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Push"));

    IpplTimings::startTimer(IpplTimings::getTimer("AMTS-EvalExternal"));
    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);
    computeExternalFields_AMTS();
    IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-EvalExternal"));

    IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Kick"));
    kick(h, true);

    // Update momentum of reference particle
    if(bends_m == 0) {
        RefPartP_suv_m = calcMeanP();
    } else {
        updateSpaceOrientation(false);
        RefPartP_suv_m = Vector_t(0.0, 0.0, sqrt(dot(RefPartP_suv_m, RefPartP_suv_m)));
    }
    RefPartP_zxy_m = dot(space_orientation_m, RefPartP_suv_m);
    IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Kick"));

    IpplTimings::startTimer(IpplTimings::getTimer("AMTS-Push"));
    if(isLastSubstep) {
        push(0.5 * h);
    } else {
        // Optimization: We include the first half push of the following substep here
        push(h);
    }
    IpplTimings::stopTimer(IpplTimings::getTimer("AMTS-Push"));

}

double ParallelTTracker::calcG() {
    if(!itsBunch->hasFieldSolver()) {
        return 1.0;
    }
    double maxAcc = 0.0;
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        // Assuming m = 1 and q = 1 (constants can be taken out of g)
        Vector_t const P = itsBunch->P[i];
        double const invGamma = 1.0 / sqrt(1.0 + dot(P, P));
        // dpdt is the force acting on the particle, derivative of the relativistic momentum
        Vector_t const dpdt = itsBunch->Ef[i] + cross(P * Physics::c * invGamma, itsBunch->Bf[i]);
        // Calculate acceleration as the second derivative of the position, i.e. derivative of
        // velocity
        Vector_t const acc = invGamma * (dpdt - P * dot(P, dpdt) * invGamma * invGamma);
        double const accSquared = dot(acc, acc);
        if(accSquared > maxAcc) {
            maxAcc = accSquared;
        }
    }
    reduce(maxAcc, maxAcc, OpMaxAssign());
    maxAcc = sqrt(maxAcc);
    double const exponent = 1.0;
    return std::pow(maxAcc, -0.5 * exponent);
}

Vector_t ParallelTTracker::calcMeanR() const {
    Vector_t meanR(0.0, 0.0, 0.0);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanR(d) += itsBunch->R[i](d);
        }
    }
    reduce(meanR, meanR, OpAddAssign());
    return meanR / Vector_t(itsBunch->getTotalNum());
}

Vector_t ParallelTTracker::calcMeanP() const {
    Vector_t meanP(0.0, 0.0, 0.0);
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        for(int d = 0; d < 3; ++d) {
            meanP(d) += itsBunch->P[i](d);
        }
    }
    reduce(meanP, meanP, OpAddAssign());
    return meanP / Vector_t(itsBunch->getTotalNum());
}
