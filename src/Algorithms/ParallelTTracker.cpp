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

#include "Algorithms/ParallelTTracker.h"

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
#include "Solvers/WakeFunction.hh"

#include "Structure/BoundaryGeometry.h"

#define PSdim 6

class PartData;

using namespace OPALTimer;
using Physics::c;

extern Inform *gmsg;
extern Inform *gmsg2all;

// Class ParallelTTracker
// ------------------------------------------------------------------------

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    itsOpalBeamline_m() {

}


ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunch &bunch,
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
    //    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    scaleFactor_m = itsBunch->getdT() * c;

    timeIntegrationTimer1_m  = IpplTimings::getTimer("TIntegration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("TIntegration2");
    timeFieldEvaluation_m  = IpplTimings::getTimer("Fieldeval");

    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");
    WakeFieldTimer_m   = IpplTimings::getTimer("WakeField");
}


ParallelTTracker::~ParallelTTracker() {
}

void ParallelTTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelTTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale) {
}

// 2007/04/19 CKR
void ParallelTTracker::execute() {
    double recpgamma, gamma;
    double t = 0.0;
    double dt = itsBunch->getdT();
    double dtTrack = dt;
    double tEmission = itsBunch->getTEmission();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    long long step;
    int gunSubTimeSteps = 10;
    int emissionSteps = 0;

    Vector_t um, a, s;
    Vector_t externalE, externalB;
    BorisPusher pusher(itsReference);
    Vector_t rmin, rmax, rtmp;

    bool global_EOL;
    unsigned long bends = 0;            // flag which indicates whether any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch

    unsigned long hasWake = 0;          // flag which indicates whether any particle is within the influence of a wake field
    int wfSection = -1;
    WakeFunction *wf;

    unsigned long hasBGeom = 0;          // flag which indicates whether any particle is within the influence of a boundary geometry
    int BGSection = -1;
    BoundaryGeometry *bgf = NULL;

    unsigned long hasSurfacePhysics = 0;
    int sphysSection = -1;
    SurfacePhysicsHandler *sphys;

    bool hasSwitchedToTEmission = false;
    bool hasSwitchedBackToTTrack = false;

    size_t totalParticles_i = itsBunch->getTotalNum();
    INFOMSG("totalParticle_i= " << totalParticles_i << endl);
    OPALTimer::Timer myt1;

    //    itsBunch->printBinHist();
    if(OPAL.inRestartRun()) {
        int prevDumpFreq = OPAL.getRestartDumpFreq();
        step = OPAL.getRestartStep() * prevDumpFreq + 1;
        // CKR if I say in my input file that my track should run for MAXSTEPS=400
        // then I want it to run for 400 steps no matter what the statement
        // above this comment yields for step
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    else if(OPAL.hasBunchAllocated() && Options::scan) {
        step = 1;
        if(!itsBunch->doEmission())
            writePhaseSpace(step - 1, 0.0); // write initial phase space
        itsBunch->setT(0.0);
    } else {
        // ADA not sure if  Options::statDumpFreq is right however
        // the starting step of the followup track off by a factor of 10
        // step = lround(t/dt)*Options::statDumpFreq;

        // CKR does not make sense to me; if I have a follow-up track I want the
        // step numbers to be continuous, not a jump by a factor of statDumpFreq
        step = OPAL.getLastStep() + 1;
        // CKR if I say in my input file that my track should run for MAXSTEPS=400
        // then I want it to run for 400 steps no matter what the statement
        // above this comment yields for step
        maxSteps_m += step;
        t = itsBunch->getT();
    }

    *gmsg << "Track start at: " << myt1.time() << " t= " << itsBunch->getT() << " zstop@ " << zstop_m << " [m]" << endl;
    if(itsBunch->doEmission()) {
        emissionSteps = static_cast<int>(itsBunch->pbin_m->getNBins()) * gunSubTimeSteps;
        *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl
              << "Change dT from " <<  itsBunch->getdT() << " [s] to "
              <<  itsBunch->getdT() << " [s] during emission " << endl;;

        if(Ippl::myNode() == 0) {
            // if we don't want opal to behave like impactt move the bunch to the cathode;
            // otherwise comment the next line out;
            RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("TSHIFT"));
            if(ar) {
                t = itsBunch->calcTimeDelay(ar->getReal());
                *gmsg << "The time delay is " << t << " [s]" << endl;
                itsBunch->setT(t);
            }
        }
    }

    // set dt for all particles already in the simulation,
    // i.e. when doing a restarted simulation
    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->dt[i] = itsBunch->getdT();
    }

    *gmsg << "Executing ParallelTTracker, initial DT " << itsBunch->getdT()
          << " [s]; max integration steps " << maxSteps_m << " step= " << step << endl;
    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(*gmsg);

    //    itsDataSink->storeFieldmaps();

    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
        itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l), itsOpalBeamline_m.getSectionStart(l));
    }


    if(!(itsBunch->weHaveBins())) {
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch->do_binaryRepart();
        IpplTimings::stopTimer(BinRepartTimer_m);
        Ippl::Comm->barrier();
    }

    // Check if there are any particles in simulation. If there are,
    // as in a restart, use the usual function to calculate beam
    // parameters. If not, calculate beam parameters of the initial
    // beam distribution.
    if(totalParticles_i == 0) {
        itsBunch->calcBeamParametersInitial();
    } else {
        itsBunch->calcBeamParameters();
    }

    RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
    RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

    if(!OPAL.hasBunchAllocated()) {
        updateSpaceOrientation(false);  // vec{p} = (0,0,p_z), vec{r} = (0,0,z)

    }

    RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
    RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();
    /* Activate all elements which influence the particles when the simulation starts;
     *  mark all elements which are already past.
     */
    /*
      increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
      before they are allocated when increasing the timestep in the gun.
    */
    itsBunch->get_bounds(rmin, rmax);
    double margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);

    bgf = itsOpalBeamline_m.getBoundaryGeometry(0);
    if(bgf) {
        *gmsg << "Found boundary geometry, end at: " << bgf->getS() << " (m) Distribution= " << bgf->getDistribution() << endl;

        Distribution *dist = Distribution::find(bgf->getDistribution());
        *gmsg << *dist << endl;

        /// this is still in BoundaryGeometry
        size_t nbparts = dist->getNumberOfDarkCurrentParticles();
        double darkinwardmargin = dist->getDarkCurrentParticlesInwardMargin();
        double workFunction = dist->getWorkFunction();
        double fieldEnhancement = dist->getFieldEnhancement();
        bgf->makeBoundaryIndexSet();
        bgf->createParticlesOnSurface(nbparts, darkinwardmargin, itsOpalBeamline_m, *itsBunch);
        bgf->setWorkFunction(workFunction);
        bgf->setFieldEnhancement(fieldEnhancement);
        dist->create(*itsBunch, *bgf);

        itsDataSink->writeGeomToVtk(*bgf, string("vtk/testGeometry-00000.vtk"));
        itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("vtk/PartlossZ-"));

        totalParticles_i = itsBunch->getTotalNum();

    }

    double minBinEmitted  = 10.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OPAL.find("MINBINEMITTED"));
    if(ar) {
        minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
        // This variable controls the number of particles to be emitted before we use
        // the space charge solver.
        *gmsg << "MINBINEMITTED " << minBinEmitted << endl;
    }

    double minStepforReBin  = 200.0;
    RealVariable *br = dynamic_cast<RealVariable *>(OPAL.find("MINSTEPFORREBIN"));
    if(br) {
        minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        *gmsg << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

    int repartFreq = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OPAL.find("REPARTFREQ"));
    if(rep) {
        repartFreq = static_cast<int>(rep->getReal());  // this variable controls the minimal number of steps until we repartition the particles
        *gmsg << "REPARTFREQ " << repartFreq << endl;
    }

    size_t totalParticles_f = 0;
    bool wakestatus = false;
    bool surfacestatus = false;
    for(step; step < maxSteps_m; ++step) {
        global_EOL = true;  //check if any particle hasn't reached the end of the field from the last element
        bends = 0;
        hasWake = 0;
        wfSection = -1;
        hasSurfacePhysics = 0;
        sphysSection = -1;

        itsOpalBeamline_m.resetStatus();

        IpplTimings::startTimer(timeIntegrationTimer1_m);

        //reset E and B to Vector_t(0.0) for every step

        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);

        for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
            //scale each particle with c*dt
            itsBunch->R[i] /= vscaleFactor;

            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

            // update local coordinate system of particle
            itsBunch->X[i] /= vscaleFactor;
            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
            itsBunch->X[i] *= vscaleFactor;
        }

        if(totalParticles_i > minBinEmitted)
            itsBunch->boundp();

        IpplTimings::stopTimer(timeIntegrationTimer1_m);

        itsBunch->calcBeamParameters();


        /** \f[ Space Charge  \f]
         */
        if(itsBunch->hasFieldSolver() && totalParticles_i > minBinEmitted && fabs(itsBunch->getChargePerParticle()) > 0.0) {
            // Do repartition if we have enough particles.
            if(totalParticles_i > 1000 && ((step % repartFreq) == 0)) {
                *gmsg << "**********************************************************" << endl;
                *gmsg << "do repartition because of repartFreq" << endl;
                *gmsg << "**********************************************************" << endl;
                IpplTimings::startTimer(BinRepartTimer_m);
                itsBunch->do_binaryRepart();
                IpplTimings::stopTimer(BinRepartTimer_m);
                Ippl::Comm->barrier();
                *gmsg << "**********************************************************" << endl;
                *gmsg << "do repartition done" << endl;
                *gmsg << "**********************************************************" << endl;
            }

            // Calculate space charge.
            if(itsBunch->weHaveBins()) {
                // When we have energy bins.
                itsBunch->calcGammas();
                for(int binNumber = 0; binNumber <= itsBunch->getLastemittedBin() && binNumber < itsBunch->getNumBins(); ++binNumber) {
                    itsBunch->setBinCharge(binNumber, itsBunch->getChargePerParticle());
                    itsBunch->computeSelfFields(binNumber);
                }
                // the next line is not needed anymore
                itsBunch->Q = itsBunch->getChargePerParticle();
            } else {
                // When we don't.
                itsBunch->computeSelfFields();
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

        if((itsBunch->weHaveBins())) {

            // switch to TEmission
            if(!hasSwitchedToTEmission) {
                dt = itsBunch->getTBin();
                itsBunch->setdT(dt);
                scaleFactor_m = dt * c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Changing emission time step to: " << dt << endl;
                hasSwitchedToTEmission = true;
            }

            int ne = 0;
            ne += itsBunch->emitParticles();

            reduce(ne, ne, OpAddAssign());
            totalParticles_i += ne;

            //emission has finished, reset to TTrack
            //if(ne == 0 && !hasSwitchedBackToTTrack) {
            //dt = dtTrack;
            //itsBunch->setdT(dt);
            //scaleFactor_m = dt * c;
            //vscaleFactor = Vector_t(scaleFactor_m);
            //*gmsg << "Emission done. Switching back to track timestep: " << dt << endl;
            //hasSwitchedBackToTTrack = true;
            //}

            if(step > minStepforReBin) {
                itsBunch->calcGammas();
                const double maxdE = abs(itsBunch->getMaxdEBins());
                if(maxdE < itsBunch->getRebinEnergy() || itsBunch->getNumBins()) {
                    if(itsBunch->getNumBins() > 1) {
                        *gmsg << "**********************************************************" << endl;
                        *gmsg << "maxdE < " << maxdE << " [keV] we rebin" << endl;
                        *gmsg << "**********************************************************" << endl;
                    } else {
                        *gmsg << "**********************************************************" << endl;
                        *gmsg << "Only one energy bin used, so we just rebin." << endl;
                        *gmsg << "**********************************************************" << endl;
                    }
                    itsBunch->rebin();
                    if(itsBunch->weHaveBins())
                        *gmsg << "--> still have bins " << endl;
                    else {
                        *gmsg << "--> ok have no bins " << endl;
                        itsBunch->setTEmission(0.0);
                    }
                    *gmsg << "**********************************************************" << endl;
                    *gmsg << "rebin  done " << endl;
                    *gmsg << "**********************************************************" << endl;
                }
            }
        } else {

            //emission has finished, reset to TTrack
            if(!hasSwitchedBackToTTrack) {
                dt = dtTrack;
                itsBunch->setdT(dt);
                scaleFactor_m = dt * c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Emission done. Switching back to track timestep: " << dt << endl;
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
            if((rtv & BEAMLINE_WAKE) && hasWake == 0) {
                wfSection = ls;
                hasWake = 1;
            }

            if((rtv & BEAMLINE_SURFACEPHYSICS) && hasSurfacePhysics == 0) {
                sphysSection = ls;
                hasSurfacePhysics = 1;
            }

            bends = bends || (rtv & BEAMLINE_BEND);

            // skip rest of the particle push if the
            // particle is out of bounds i.e. does not see
            // a E or B field
            if(rtv & BEAMLINE_OOB)
                itsBunch->Bin[i] = -1;

            itsBunch->Ef[i] += externalE;
            itsBunch->Bf[i] += externalB;

            itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
        }


        reduce(hasWake, hasWake, OpAddAssign());
        reduce(hasSurfacePhysics, hasSurfacePhysics, OpAddAssign());
        reduce(bends, bends, OpAddAssign());

        if(hasWake > 0) {
            if(!wakestatus) {
                *gmsg << "============== START WAKE CALCULATION =============" << endl;
                wakestatus = true;
            }
            reduce(wfSection, wfSection, OpMaxAssign());
            wf = itsOpalBeamline_m.getWakeFunction(wfSection);
        } else if(wakestatus) {
            *gmsg << "=============== END WAKE CALCULATION ==============" << endl;
            wakestatus = false;
        }

        if(hasSurfacePhysics > 0) {
            if(!surfacestatus) {
                *gmsg << "============== START SURFACE PHYSICS CALCULATION =============" << endl;
                surfacestatus = true;
            }
            reduce(sphysSection, sphysSection, OpMaxAssign());
            sphys = itsOpalBeamline_m.getSurfacePhysicsHandler(sphysSection);
        } else if(surfacestatus) {
            *gmsg << "============== END SURFACE PHYSICS CALCULATION =============" << endl;
            surfacestatus = false;
        }

        IpplTimings::stopTimer(timeFieldEvaluation_m);

        if(itsBunch->getLocalNum() == 0)
            global_EOL = false;

        /**
           Delete marked particles
        */

        bool globPartOutOfBounds = (min(itsBunch->Bin) < 0);
        size_t ne = 0;
        if(globPartOutOfBounds) {
            ne = itsBunch->boundp_destroyT();
        }
        totalParticles_f = totalParticles_i - ne;
        if(hasWake > 0) {
            IpplTimings::startTimer(WakeFieldTimer_m);
            if(!wf) {
                INFOMSG("no wakefunction attached" << endl);
            } else {
                wf->apply(*itsBunch);
            }
            IpplTimings::stopTimer(WakeFieldTimer_m);

        }

        if(hasSurfacePhysics > 0) {
            if(!sphys) {
                INFOMSG("no surface physics attached" << endl);
            } else {
                sphys->apply(*itsBunch);
            }
        }

        kickParticles(pusher);

        if(bends == 0) {
            if(totalParticles_f > 0) {
                // none of the particles is in a bending element
                updateReferenceParticle();
            }
        } else {
            /* at least one of the elements bends the beam; until all particles have left the bending elements we track the reference particle
                      * as if it were a regular particle; from the moment when the reference particle has reached the bending field until it leaves
                      * it again we rotate the bunch about the position of the reference particle such that the momentum of the reference particle points
                      * in z direction
                      */

            RefPartP_suv_m = itsBunch->get_pmean();
            RefPartR_suv_m = itsBunch->get_rmean();
            recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
            updateSpaceOrientation(false);

            // First update the momentum of the reference particle in zxy coordinate system, then update its position
            RefPartP_zxy_m(0) = SpaceOrientation_m[0] * RefPartP_suv_m(0) + SpaceOrientation_m[1] * RefPartP_suv_m(1) + SpaceOrientation_m[2] * RefPartP_suv_m(2);
            RefPartP_zxy_m(1) = SpaceOrientation_m[3] * RefPartP_suv_m(0) + SpaceOrientation_m[4] * RefPartP_suv_m(1) + SpaceOrientation_m[5] * RefPartP_suv_m(2);
            RefPartP_zxy_m(2) = SpaceOrientation_m[6] * RefPartP_suv_m(0) + SpaceOrientation_m[7] * RefPartP_suv_m(1) + SpaceOrientation_m[8] * RefPartP_suv_m(2);
            RefPartR_zxy_m += RefPartP_zxy_m * recpgamma * scaleFactor_m / 2.;

            RefPartP_suv_m = Vector_t(0.0, 0.0, sqrt(dot(RefPartP_suv_m, RefPartP_suv_m)));
            RefPartR_suv_m += RefPartP_suv_m * recpgamma / 2.;
            RefPartR_suv_m *= vscaleFactor;
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

        // start particle loop part 2
        for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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

        if(totalParticles_i != totalParticles_f) {
            *gmsg << "======================================================================" << endl;
            *gmsg << "Particle loss at t= " << t << " n= " << totalParticles_i - totalParticles_f << " Ntot= " << totalParticles_f << endl;
            *gmsg << "======================================================================" << endl;
            totalParticles_i = totalParticles_f;
        }

        t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time

        itsBunch->setT(t);


        if(bgf) {
            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();


            /**
               Here we check if a particles is
               outside the domain, maybe flag it for
               deletion later and create secondaries
            */
            for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
                Vector_t intecoords = outr;
                int triId = 0;
                int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], itsBunch->getdT(), itsBunch->PType[i], intecoords, triId);
                /*
                  check physics here
                  bgf->doBGphysics(intecoords,triId);
                  unsigned long OpalBeamline::getFieldAt(const unsigned int& index, const Vector_t& pos, const long& sindex, const double& t, Vector_t& E, Vector_t& B)

                */
                if(res == 0)
                    res += bgf->doBGphysics(intecoords, triId);

                if(res == 0)
                    itsBunch->Bin[i] = -1;
            }

            /*
              Now we do again fieldemission
            */
            Distribution *dist = Distribution::find(bgf->getDistribution());
            size_t nbparts = dist->getNumberOfDarkCurrentParticles();
            double darkinwardmargin = dist->getDarkCurrentParticlesInwardMargin();
            bgf->createParticlesOnSurface(nbparts, darkinwardmargin, itsOpalBeamline_m, *itsBunch);
            dist->create(*itsBunch, *bgf);

            //            bgf->callFNemission(itsOpalBeamline_m, itsBunch, t);
            totalParticles_i = itsBunch->getTotalNum();
            *gmsg << *itsBunch << endl;
        }



        //IFF: cheap step dump regulation
        OPALTimer::Timer myt2;
        double sposRef = 0.0;
        if(totalParticles_f > 0) {
            sposRef = itsBunch->get_sPos();
            if(totalParticles_f <= minBinEmitted) {
                *gmsg << myt2.time() << " Step " << step << "; only " << totalParticles_f << " particles emitted; t= " << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
            } else if(isnan(sposRef) || isinf(sposRef)) {
                *gmsg << myt2.time() << " Step " << step << "; there seems to be something wrong with the position of the bunch!" << endl;
            } else {
                *gmsg << myt2.time() << " Step " << step << " at " << sposRef << " [m] t= "
                      << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
                if(step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0) {
                    writePhaseSpace(step, sposRef);
                    if(bgf) {
                        itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("vtk/PartlossZ-"));
                    }
                }
            }
            /**
               Stop simulation if beyond zstop_m
            */
            if(sposRef > zstop_m) {
                maxSteps_m = step;
            }
        } else {
            *gmsg << "Step " << step << " no emission yet "  << " t= " << t << " [s]" << endl;
        }

        //         if (sposRef > zstop_m) break;

        if(step > emissionSteps) {
            //TEST the non-len reduce: reduce(&global_EOL, &global_EOL, OpBitwiseOrAssign());
            reduce(&global_EOL, &global_EOL + 1, &global_EOL, OpBitwiseAndAssign());
            if(global_EOL) {
                break;
            }
        }
        // this seams to fix Ticket #12
        //  Ippl::Comm->barrier();
    }

    OPALTimer::Timer myt3;
    OPAL.setLastStep(step);

    itsOpalBeamline_m.switchElementsOff();
    //itsOpalBeamline_m.switchElementsOff(numeric_limits<double>::min(), numeric_limits<double>::max());

    *gmsg << "done executing ParallelTTracker at " << myt3.time() << endl;
}

// 2007/04/19 CKR
