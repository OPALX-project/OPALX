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


#include "Algorithms/ParallelTTracker.h"

#include "AbstractObjects/OpalData.h"

#include "BasicActions/Option.h"

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
#include "Utilities/OpalException.h"
#include "Structure/BoundaryGeometry.h"

#define PSdim 6

class PartData;

using namespace std;
using namespace OPALTimer;

extern Inform *gmsg;
extern Inform *gmsg2all;

// Class ParallelTTracker
// ------------------------------------------------------------------------

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
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
    maxSteps_m(maxSTEPS),
    zstop_m(zstop) {
    //    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    scaleFactor_m = itsBunch->getdT() * Physics::c;
    timeIntegrationTimer1_m  = IpplTimings::getTimer("TIntegration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("TIntegration2");
    timeFieldEvaluation_m  = IpplTimings::getTimer("Fieldeval");

    BinRepartTimer_m   = IpplTimings::getTimer("Binaryrepart");
    WakeFieldTimer_m   = IpplTimings::getTimer("WakeField");
#ifdef DBG_SYM
    string SfileName = OpalData::getInstance()->getInputFn();
    int pdot = SfileName.find(string("."), 0);
    SfileName.erase(pdot, SfileName.size() - pdot);
    string fn = SfileName + string(".fields");
    of_m.open(fn.c_str(), ios::out);
    of_m.precision(9);
    of_m << "# spos Ex Ey Ez Bz By Bz at: (h,h),(h,-h),(-h,h)(-h,-h) h=0.001" << endl;
#endif
}


ParallelTTracker::~ParallelTTracker() {
#ifdef DBG_SYM
    of_m.close();
#endif
}

void ParallelTTracker::applySchottkyCorrection(PartBunch &itsBunch, int ne, double t, double rescale_coeff) {

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

    *gmsg << "*****************************************************************" << endl;
    *gmsg << " Estimate Schottky correction                                    " << endl;
    *gmsg << "*****************************************************************" << endl;

    double margin = 0.0;
    if(!mpacflg_m) {
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
        margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        margin = 0.01 > margin ? 0.01 : margin;
        itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);
    }

    double minBinEmitted  = 10.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar) {
        minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
        // This variable controls the number of particles to be emitted before we use
        // the space charge solver.
        *gmsg << "MINBINEMITTED " << minBinEmitted << endl;
    }


    double minStepforReBin  = 10000.0;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br) {
        minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        *gmsg << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

    int repartFreq = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep) {
        repartFreq = static_cast<int>(rep->getReal());  // this variable controls the minimal number of steps until we repartition the particles
        *gmsg << "REPARTFREQ " << repartFreq << endl;
    }

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq = 1000000;

    size_t totalParticles_f = 0;

    for(; step < maxSteps_m; ++step) {
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
            if(itsBunch->weHaveBins()) {
                // When we have energy bins.
                itsBunch->calcGammas();
                ParticleAttrib<double> Q_back = itsBunch->Q;
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

        if((itsBunch->weHaveBins())) {

            // switch to TEmission
            if(!hasSwitchedToTEmission) {
                dt = itsBunch->getTBin();
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Changing emission time step to: " << dt << endl;
                hasSwitchedToTEmission = true;
            }

            int ne = 0;
            ne += itsBunch->emitParticles();

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
                *gmsg << "Emission done. Switching back to track timestep: " << dt << endl;
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
            *gmsg << "* Deleted " << ne << " particles, remaining " << totalParticles_f << " particles" << endl; //benchmark output

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
               Stop simulation if beyond zstop_m
            */
            if(sposRef > zstop_m) {
                maxSteps_m = step;
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
    *gmsg << "done executing Schottky loop " << myt3.time() << endl;
    return itsBunch->getCharge();
}


pair<FieldList::iterator , bool> ParallelTTracker::checkCavity(double s) {

    pair<FieldList::iterator , bool> res(cavities_m.begin(), false);
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit)
        if(((*fit).getStart() <= s) && (s <= (*fit).getEnd())) {
            res.first = fit;
            res.second = true;
            break;
        }
    return res;
}

pair<FieldList::iterator, bool> ParallelTTracker::doOneStep(BorisPusher pusher) {
    bool global_EOL = true;  //check if any particle hasn't reached the end of the field from the last element
    unsigned long bends = 0;

    double recpgamma;
    double t = itsBunch->getT();

    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    Vector_t externalE, externalB;

    Vector_t rmin, rmax;

    itsOpalBeamline_m.resetStatus();

    //reset E and B to Vector_t(0.0) for every step

    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        //scale each particle with c*dt
        itsBunch->R[i] /= vscaleFactor;

        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

        // update local coordinate system of particle
        itsBunch->X[i] /= vscaleFactor;
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i],
                                                itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->dt[0]);
        itsBunch->X[i] *= vscaleFactor;
    }
    //    itsBunch->calcBeamParameters();

    // push the reference particle by a half step
    recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;

    //
    // get external fields for all particles
    //
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        //FIXME: rethink scaling!
        itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c
                                   * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

        long ls = itsBunch->LastSection[i];
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], ls);
        if(ls != itsBunch->LastSection[i]) {
            if(!itsOpalBeamline_m.section_is_glued_to(itsBunch->LastSection[i], ls)) {
                itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls),
                                                     itsOpalBeamline_m.getSectionStart(ls));
            }
            itsBunch->LastSection[i] = ls;
        }
        const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch->R[i], ls,
                                  t + itsBunch->dt[i] / 2., externalE, externalB);
        global_EOL = global_EOL && (rtv & BEAMLINE_EOL);

        bends = bends || (rtv & BEAMLINE_BEND);

        // skip rest of the particle push if the
        // particle is out of bounds i.e. does not see
        // a E or B field
        if(rtv & BEAMLINE_OOB)
            itsBunch->Bin[i] = -1;

        itsBunch->Ef[i] += externalE;
        itsBunch->Bf[i] += externalB;

        itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c
                                   * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
    }

    kickParticlesAutophase(pusher);

    if(bends == 0) {
        updateReferenceParticleAutophase();
    } else {
        /* at least one of the elements bends the beam; until all particles
         * have left the bending elements we track the reference particle
         * as if it were a regular particle; from the moment when the reference
         * particle has reached the bending field until it leaves
         * it again we rotate the bunch about the position of the reference
         * particle such that the momentum of the reference particle points
         * in z direction
         */
        RefPartP_suv_m = itsBunch->P[0];
        RefPartR_suv_m = itsBunch->R[0];
        recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        updateSpaceOrientation(false);

        // First update the momentum of the reference particle in zxy coordinate system, then update its position
        RefPartP_zxy_m = dot(space_orientation_m, RefPartP_suv_m);
        RefPartR_zxy_m += RefPartP_zxy_m * recpgamma * scaleFactor_m / 2.;

        RefPartP_suv_m = Vector_t(0.0, 0.0, sqrt(dot(RefPartP_suv_m, RefPartP_suv_m)));
        RefPartR_suv_m += RefPartP_suv_m * recpgamma / 2.;
        RefPartR_suv_m *= vscaleFactor;
    }

    itsBunch->RefPart_R = RefPartR_zxy_m;
    itsBunch->RefPart_P = RefPartP_zxy_m;

    // calculate the dimensions of the bunch and add a small margin to them;
    // then decide which elements have to be triggered
    // when an element is triggered memory is allocated and the field map is read in
    rmin = rmax = itsBunch->R[0];

    // trigger the elements
    double margin = 3. * RefPartP_suv_m(2) * recpgamma;
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements((rmin(2) - margin)*scaleFactor_m, (rmax(2) + margin)*scaleFactor_m, true);

    // start particle loop part 2
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        /** \f[ \vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{1}{2}\vec{v}_{n+1/2}\quad (= \vec{x}_{n+1/2} +
            \frac{\Delta t}{2} \frac{\vec{\beta}_{n+1/2}\gamma_{n+1/2}}{\gamma_{n+1/2}}) \f]
            * \code
            * R[i] += 0.5 * P[i] * recpgamma;
            * \endcode
            */
        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
        //and scale back to dimensions
        itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i],
                                   Physics::c * itsBunch->dt[i]);
        // update local coordinate system
        itsBunch->X[i] /= vscaleFactor;
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i],
                                                itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->dt[0]);
        itsBunch->X[i] *= vscaleFactor;
    }

    t += itsBunch->dt[0]; //t after a full global timestep with dT "synchronization point" for simulation time

    itsBunch->setT(t);

    return checkCavity(itsBunch->R[0](2));
}


void ParallelTTracker::applyEntranceFringe(double angle, double curve,
        const BMultipoleField &field, double scale) {
}


void ParallelTTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale) {
}


void ParallelTTracker::showCavities(Inform &m) {

    m << "Found the following cavities:" << endl;

    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        m << (*fit).getElement()->getName()
          << " from " << (*fit).getStart() << " to "
          << (*fit).getEnd() << " (m) phi=";
        if((*fit).getElement()->getType() == "TravelingWave")
            m << static_cast<TravelingWave *>((*fit).getElement())->getPhasem() / Physics::pi * 180.0 << endl;
        else
            m << static_cast<RFCavity *>((*fit).getElement())->getPhasem() / Physics::pi * 180.0 << endl;
    }
    m << endl << endl;
}


void ParallelTTracker::updateRFElement(string elName, double maxPhi) {
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
    double phi = 0;
    double freq = 0.0;
    const double RADDEG = 1.0 / Physics::pi * 180.0;
    //   Inform m ("updateALLRFElements ",INFORM_ALL_NODES);
    *gmsg << "\n-------------------------------------------------------------------------------------\n";
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if(fit != cavities_m.begin())
            *gmsg << "\n";
        if((*fit).getElement()->getType() == "TravelingWave") {
            freq = static_cast<TravelingWave *>((*fit).getElement())->getFrequencym();
            phi = static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
            *gmsg << (*fit).getElement()->getName()
                  << ": phi_orig= phi_nom + phi_maxE= " << phi *RADDEG << " degree, "
                  << "global phase shift= " << -phiShift *freq *RADDEG << " degree\n";
            phi -= (phiShift * freq);
            static_cast<TravelingWave *>((*fit).getElement())->updatePhasem(phi);
        } else {
            freq = static_cast<RFCavity *>((*fit).getElement())->getFrequencym();
            phi = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
            *gmsg << (*fit).getElement()->getName()
                  << ": phi_orig= phi_nom + phi_maxE= " << phi *RADDEG << " degree, "
                  << "global phase shift= " << -phiShift *freq *RADDEG << " degree\n";
            phi -= (phiShift * freq);
            static_cast<RFCavity *>((*fit).getElement())->updatePhasem(phi);
        }
    }
    *gmsg << "-------------------------------------------------------------------------------------\n"
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
        if(zbegin < zstop_m && zend >= zstop_m) {
            *gmsg << "\033[0;31m"
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
    travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
    cavities_m.merge(travelingwaves_m, OpalField::SortAsc);

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


    if(Ippl::myNode() == 0) {
        itsBunch->create(1);
        itsBunch->R[0] = iniR;
        itsBunch->P[0] = iniP;
        itsBunch->Bin[0] = 0;
        itsBunch->Q[0] = itsBunch->getChargePerParticle();
        itsBunch->PType[0] = 0;
        itsBunch->LastSection[0] = 0;

        executeAutoPhase(Options::autoPhase, zStop);
        itsBunch->destroy(1, 0);

        // need to rebuild for updateAllRFElements
        cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
        travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
        cavities_m.merge(travelingwaves_m, OpalField::SortAsc);


        // now send all max phases and names of the cavities to
        // all the other nodes for updating.
        Message *mess = new Message();
        putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

        for(vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
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
            string elName;
            double maxPhi;
            getMessage(*mess, elName);
            getMessage(*mess, maxPhi);
            updateRFElement(elName, maxPhi);
            OpalData::getInstance()->setMaxPhase(elName, maxPhi);
        }
    }

    OpalData::getInstance()->setGlobalPhaseShift(gPhaseSave);
    return cavities_m;
}


void ParallelTTracker::executeAutoPhase(int numRefs, double zStop) {
    Inform msg("Autophasing ");

    const double RADDEG = 180.0 / Physics::pi;

    Vector_t rmin, rmax;

    size_t step = 0;
    int dtfraction = 2;
    itsBunch->dt = itsBunch->getdT() / dtfraction;         // need to fix this and make the factor 2 selectable

    double scaleFactorSave = scaleFactor_m;
    scaleFactor_m = itsBunch->dt[0] * Physics::c;

    double tSave = itsBunch->getT();


    Vector_t vscaleFactor = Vector_t(scaleFactor_m);

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

    for(; step < maxSteps_m * dtfraction; ++step) {

        itsBunch->setTrackStep(step);

        // let's do a drifting step to probe if the particle will reach element in next step
        Vector_t R_drift = itsBunch->R[0] + itsBunch->P[0] / sqrt(1.0 + dot(itsBunch->P[0],
                           itsBunch->P[0])) * vscaleFactor;

        pair<FieldList::iterator, bool> res = checkCavity(R_drift(2));

        if(res.second) {
            double orig_phi = 0.0;
            double Phimax, Emax = 0.0;
            double PhiAstra;
            //////
            const double beta = sqrt(1. - 1 / (itsBunch->P[0](2) * itsBunch->P[0](2) + 1.));
            const double tErr  = ((*res.first).getStart() - itsBunch->R[0](2)) / (Physics::c * beta);

            bool apVeto;

            INFOMSG("Found " << (*res.first).getElement()->getName()
                    << " at " << itsBunch->R[0](2) << " [m], "
                    << "step  " << step << ", "
                    << "t= " << itsBunch->getT() << " [s],\n"
                    << "E= " << getEnergyMeV(itsBunch->P[0]) << " [MeV]\n"
                    << "start phase scan ... " << endl);

            INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
            if((*res.first).getElement()->getType() == "TravelingWave") {
                orig_phi = static_cast<TravelingWave *>((*res.first).getElement())->getPhasem();
                apVeto = static_cast<TravelingWave *>((*res.first).getElement())->getAutophaseVeto();
                if(apVeto)
                    msg << " ----> APVETO -----> "
                        << static_cast<TravelingWave *>((*res.first).getElement())->getName() <<  endl;

                INFOMSG((*res.first).getElement()->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);

                if(apVeto) {
                    Phimax = orig_phi;
                } else {
                    TravelingWave *element = static_cast<TravelingWave *>((*res.first).getElement());
                    Phimax = element->getAutoPhaseEstimate(getEnergyMeV(itsBunch->P[0]),
                                                           itsBunch->getT() + tErr,
                                                           itsReference.getQ(),
                                                           itsReference.getM() * 1e-6);
                }
            } else {
                orig_phi = static_cast<RFCavity *>((*res.first).getElement())->getPhasem();
                apVeto = static_cast<RFCavity *>((*res.first).getElement())->getAutophaseVeto();
                if(apVeto)
                    msg << " ----> APVETO -----> "
                        << static_cast<RFCavity *>((*res.first).getElement())->getName() << endl;

                INFOMSG((*res.first).getElement()->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);

                if(apVeto) {
                    Phimax = orig_phi;
                } else {
                    RFCavity *element = static_cast<RFCavity *>((*res.first).getElement());
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

            double E = APtrack(res, phi);

            msg << "Did APtrack with phi= " << phi << " result E= " << E << endl;

            if(!apVeto) {
                INFOMSG("do fine scan around effective max energy (" << E << " MeV)" << ", dphi= " << dphi << endl;);
                do {
                    j ++;
                    Emax = E;
                    Phiini = phi;
                    phi -= dphi;
                    INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                    E = APtrack(res, phi);
                    if(E > Emax) {
                        INFOMSG(E - Emax << " MeV: accepted" << endl;);
                    } else {
                        INFOMSG(E - Emax << " MeV: rejected" << " E= " << E << " Emax= " << Emax << endl;);
                    }
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
                        INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                        E = APtrack(res, phi);
                        if(E > Emax) {
                            INFOMSG(E - Emax << " MeV: accepted" << endl;);
                        } else {
                            INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        }
                    } while(E > Emax);
                }
                for(int refinement_level = 0; refinement_level < numRefs; refinement_level ++) {
                    dphi /= 2.;
                    INFOMSG("refinement level: " << refinement_level + 1 << ", dphi= " << dphi << endl;);
                    phi = Phiini - dphi;
                    INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                    E = APtrack(res, phi);
                    if(E > Emax) {
                        INFOMSG(E - Emax << " MeV: accepted" << endl;);
                        Phiini = phi;
                        Emax = E;
                    } else {
                        INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        phi = Phiini + dphi;
                        INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                        E = APtrack(res, phi);
                        if(E > Emax) {
                            INFOMSG(E - Emax << " MeV: accepted" << endl;);
                            Phiini = phi;
                            Emax = E;
                        } else {
                            INFOMSG(E - Emax << " MeV: rejected" << endl;);
                        }
                    }
                }
                Phimax = Phiini;
                INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
            } else {
                Phimax = 0.0;
            }


            if((*res.first).getElement()->getType() == "TravelingWave") {
                static_cast<TravelingWave *>((*res.first).getElement())->updatePhasem(Phimax + orig_phi);
            } else {
                static_cast<RFCavity *>((*res.first).getElement())->updatePhasem(Phimax + orig_phi);
            }

            PhiAstra = (Phimax * RADDEG) + 90.0;
            PhiAstra -= floor(PhiAstra / 360.) * 360.;

            msg << (*res.first).getElement()->getName() << "_phi= "  << Phimax << " rad / "
                << Phimax *RADDEG <<  " deg, AstraPhi= " << PhiAstra << " deg,\n"
                << "E= " << Emax << " (MeV), " << "phi_nom= " << orig_phi *RADDEG << endl;

            maxPhases_m.push_back(MaxPhasesT((*res.first).getElement()->getName(), Phimax));
            OpalData::getInstance()->setMaxPhase((*res.first).getElement()->getName(), Phimax);
            cavities_m.erase(res.first);
        }

        doOneStep(pusher);
        double sposRef = itsBunch->R[0](2);

        if(sposRef > zStop)
            maxSteps_m = floor(step / dtfraction);

        if(!(step % 1000)) {
            INFOMSG("step = " << step << ", spos = " << sposRef << " [m], t= " << itsBunch->getT() << " [s], "
                    << "E= " << getEnergyMeV(itsBunch->P[0]) << " [MeV] " << endl);
        }
    }

    scaleFactor_m = scaleFactorSave;
    itsBunch->setT(tSave);
}

// 2007/04/19 CKR
void ParallelTTracker::execute() {
    //    Inform m("ParallelTTrackDebug ", INFORM_ALL_NODES);
    double recpgamma;
    double t = 0.0;
    double dt = itsBunch->getdT();
    double dtTrack = dt;
    double tEmission = itsBunch->getTEmission();
    Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    unsigned long long step = 0;
    int gunSubTimeSteps = 10;
    unsigned int emissionSteps = 0;

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
    WakeFunction *wf = NULL;

    BoundaryGeometry *bgf = NULL;
    int secondaryFlg = 0;
    size_t maxNparts = 100000000;        // upper limit of particle number when we do field emission and secondary emission simulation. Could be reset to another value in input file with MAXPARTSNUM.
    bool nEmissionMode = true;

    unsigned long hasSurfacePhysics = 0;
    int sphysSection = -1;
    SurfacePhysicsHandler *sphys = NULL;

    bool hasSwitchedToTEmission = false;
    bool hasSwitchedBackToTTrack = false;

    double gPhaseSave = 0.0;


    size_t numberOfFieldEmittedParticles = 0;


    if(Options::autoPhase > 0) {
        gPhaseSave = OpalData::getInstance()->getGlobalPhaseShift();
        OpalData::getInstance()->setGlobalPhaseShift(0.0);
    }

    itsBeamline_m.accept(*this);

    // make sure that no monitor has overlap with two tracks
    FieldList monitors = itsOpalBeamline_m.getElementByType("Monitor");
    for(FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
        double zbegin, zend;
        it->getElement()->getDimensions(zbegin, zend);
        if(zbegin < zstop_m && zend >= zstop_m) {
            *gmsg << "\033[0;31m"
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
    travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
    cavities_m.merge(travelingwaves_m, OpalField::SortAsc);

    /**
       AUTOPHASE>0: do autophasing before tracking
       without a global phase shift!
    */

    if((!OpalData::getInstance()->inRestartRun()) && (Options::autoPhase > 0)) {
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

        size_t Nloc = itsBunch->getLocalNum();
        if(!OpalData::getInstance()->hasBunchAllocated() && Nloc > 0) {
            iniR = itsBunch->get_rmean();
            iniP = itsBunch->get_pmean();
            id.create(Nloc);
            id = itsBunch->ID;
            r.create(Nloc);
            r = itsBunch->R;
            p.create(Nloc);
            p = itsBunch->P;
            x.create(Nloc);
            x = itsBunch->X;
            q.create(Nloc);
            q = itsBunch->Q;
            bin.create(Nloc);
            bin = itsBunch->Bin;
            dt.create(Nloc);
            dt = itsBunch->dt;
            ls.create(Nloc);
            ls = itsBunch->LastSection;
            ptype.create(Nloc);
            ptype = itsBunch->PType;

            itsBunch->destroy(Nloc, 0);
            itsBunch->update();
        }
        if(Ippl::myNode() == 0) {
            double zStop = itsOpalBeamline_m.calcBeamlineLenght();
            if(!OpalData::getInstance()->hasBunchAllocated()) {
                itsBunch->create(1);
                itsBunch->R[0] = iniR;
                itsBunch->P[0] = iniP;
                itsBunch->Bin[0] = 0;
                itsBunch->Q[0] = itsBunch->getChargePerParticle();
                itsBunch->PType[0] = 0;
                itsBunch->LastSection[0] = 0;
                executeAutoPhase(Options::autoPhase, zStop);
                itsBunch->destroy(1, 0);
                // need to rebuild for updateAllRFElements
                cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
                travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
                cavities_m.merge(travelingwaves_m, OpalField::SortAsc);

            } else {
                // we are in a followup track and the phase information is
                // already stored in the OPAL dictionary.
                for(vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
                    updateRFElement((*it).first, (*it).second);
                    INFOMSG("In follow-up track use saved phases for -> name: " << (*it).first << " phi= " << (*it).second << " (rad)" << endl);
                }
            }

            // now send all max phases and names of the cavities to
            // all the other nodes for updating.
            Message *mess = new Message();
            putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

            for(vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
                putMessage(*mess, (*it).first);
                putMessage(*mess, (*it).second);
            }
            Ippl::Comm->broadcast_all(mess, tag);
        } else {
            // receive max phases and names and update the structure
            int nData = 0;
            Message *mess = Ippl::Comm->receive_block(Parent, tag);
            getMessage(*mess, nData);
            for(int i = 0; i < nData; i++) {
                string elName;
                double maxPhi;
                getMessage(*mess, elName);
                getMessage(*mess, maxPhi);
                updateRFElement(elName, maxPhi);
                OpalData::getInstance()->setMaxPhase(elName, maxPhi);
            }
        }

        if(!OpalData::getInstance()->hasBunchAllocated() && Nloc > 0) {
            itsBunch->update();

            itsBunch->create(Nloc);

            itsBunch->ID = id;
            itsBunch->R = r;
            itsBunch->P = p;
            itsBunch->X = x;
            itsBunch->Q = q;
            itsBunch->Bin = bin;
            itsBunch->dt = dt;
            itsBunch->LastSection = ls;
            itsBunch->PType = ptype;
        }

        itsBunch->update();

    } else if(OpalData::getInstance()->inRestartRun() && Options::autoPhase > 0) {
        itsDataSink->retriveCavityInformation(OpalData::getInstance()->getInputFn());

        for(vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++)
            updateRFElement((*it).first, (*it).second);
    }
    OpalData::getInstance()->setGlobalPhaseShift(gPhaseSave);
    //    if((OpalData::getInstance()->getGlobalPhaseShift() > 0.0)  && (Options::autoPhase > 0))
    if(Options::autoPhase > 0)
        updateAllRFElements(OpalData::getInstance()->getGlobalPhaseShift());

    if((!OpalData::getInstance()->inRestartRun()) && (Options::autoPhase > 0) && (!OpalData::getInstance()->hasBunchAllocated()))
        itsDataSink->storeCavityInformation();

    size_t totalParticles_i = itsBunch->getTotalNum();
    *gmsg << "totalParticle_i= " << totalParticles_i << endl;
    OPALTimer::Timer myt1;

    if(OpalData::getInstance()->inRestartRun()) {
        int prevDumpFreq = OpalData::getInstance()->getRestartDumpFreq();
        step = OpalData::getInstance()->getRestartStep() * prevDumpFreq + 1;
        t = itsBunch->getT();
    }

    else if(OpalData::getInstance()->hasBunchAllocated() && Options::scan) {
        step = 1;
        if(itsBunch->getLocalNum() != 0)
            writePhaseSpace(step - 1, 0.0); // write initial phase space
        itsBunch->setT(0.0);
    } else {
        step = OpalData::getInstance()->getLastStep();
        t = itsBunch->getT();
    }

    *gmsg << "Track start at: " << myt1.time() << ", t= " << itsBunch->getT() << "; zstop at: " << zstop_m << " [m]" << endl;

    if(!mpacflg_m) {
        if(itsBunch->doEmission()) {
            emissionSteps = static_cast<unsigned int>(itsBunch->pbin_m->getNBins()) * gunSubTimeSteps;
            *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl
                  << "Change dT from " <<  itsBunch->getdT() << " [s] to "
                  <<  itsBunch->getdT() << " [s] during emission " << endl;;
        }

        // set dt for all particles already in the simulation,
        // i.e. when doing a restarted simulation
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            itsBunch->dt[i] = itsBunch->getdT();
        }
    }

    double init_erg = itsBunch->getEkin();
    double tol_iter = 1e-5;
    double rescale_coeff = 1 / init_erg / init_erg;
    if(Options:: schottkyRennormalization > 0) {
        rescale_coeff = Options:: schottkyRennormalization;
        *gmsg << "Set schottky scale coefficient to  " << rescale_coeff << endl;
    } else if(Options::schottkyCorrection) {
        while(true) {
            double real_charge = schottkyLoop(rescale_coeff);

            double total_charge = itsBunch->getTotalNum() * itsBunch->getChargePerParticle();
            *gmsg << "Schottky scale coefficient " << rescale_coeff << ", actual emitted charge " << real_charge << " (Cb)" << endl;
            itsBunch->cleanUpParticles();
            itsBunch->setT(0);
            double scale_error = total_charge / real_charge - 1;
            // TODO : send scale_error to all nodes
            rescale_coeff *= (1.3 * scale_error + 1);
            if(fabs(scale_error) < tol_iter)
                break;
        }
        *gmsg << "Schottky scan, final scale coefficient " << rescale_coeff << " ()" << endl;
    }

    *gmsg << "Executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s];\n"
          << "max integration steps " << maxSteps_m << ", next step= " << step << endl;

    //    itsBeamline_m.accept(*this);
    //    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(*gmsg);
    double margin = 0.0;
    if(!mpacflg_m) {
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
        margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        margin = 0.01 > margin ? 0.01 : margin;
        itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);
    }

    //bgf = itsOpalBeamline_m.getBoundaryGeometry(0);
    for(unsigned int i = 0; i < itsOpalBeamline_m.sections_m.size(); i++) {

        bgf = itsOpalBeamline_m.getBoundaryGeometry(i);
        if(bgf) {
            Distribution *dist = NULL;
            Distribution *distrand = NULL;
            vector<string> distr_str = bgf->getDistributionArray();
            if(distr_str.size() == 0) {
                string distr = bgf->getDistribution();
                if(!distr.empty()) {
                    *gmsg << "* Find boundary geometry, start at: " << bgf->getS() << " (m) Distribution= " << bgf->getDistribution() << endl;
                    dist = Distribution::find(bgf->getDistribution());
                    *gmsg << "* " << *dist << endl;
                } else {
                    throw OpalException("ParallelTTracker::execute()",
                                        "No distribution attached to BoundaryGeometry. Please check the input file... ...");

                }
            } else {
                *gmsg << "************************************************************************************************* " << endl;
                *gmsg <<  "* Find boundary geometry, start at: " << bgf->getS()  << " (m). " << endl;
                *gmsg << "* Attached more than one distribution: " << endl;
                for(vector<string>::const_iterator dit = distr_str.begin(); dit != distr_str.end(); ++ dit) {
                    Distribution *d = Distribution::find(*dit);
                    *gmsg << "* Distribution: " << *dit << " distribution type: " << d->getTypeofDistribution() << endl;
                    *gmsg << "************************************************************************************************* " << endl;
                    if(d->getTypeofDistribution() == "SURFACEEMISSION") {
                        dist = d;
                        *gmsg << *dist << endl;

                    } else if(d->getTypeofDistribution() == "SURFACERANDCREATE") {
                        distrand = d;
                        *gmsg << *distrand << endl;
                        // here nbparts should be non zero as these particles will be the initialization of primary bunch.
                        size_t nbparts = distrand->getNumberOfDarkCurrentParticles();
                        double darkinwardmargin = distrand->getDarkCurrentParticlesInwardMargin();
                        double einitthreshold = distrand->getEInitThreshold();
                        // make sure that the elements have already been switched on before initializing particles in position where the electric field > einitthreshold.
                        bgf->setEInitThreshold(einitthreshold);
                        if(!mpacflg_m) {
                            bgf->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
                            distrand->createPriPart(itsBunch, *bgf);
                            totalParticles_i = itsBunch->getTotalNum();
                        } else {
                            /*
                            Multipacting flag set true. Generate primary particles.
                                                 Activate all elements (switch on the field map of elements in multipacting) in multipacting simulation
                                               */

                            itsOpalBeamline_m.switchAllElements();
                            // it is possible to generate initial particles according to E field, since all elements switched on before we create particles.
                            bgf->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
                            // for Parallel Plate benchmark, Vw should be defined in input file and will be invoked by getVw method in createPriPart().
                            // for other multipacting simulation no need to define the Vw in SURFACERANDCREATE in input file.
                            distrand->createPriPart(itsBunch, *bgf);
                            totalParticles_i = itsBunch->getTotalNum();
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

                            if(totalParticles_i == 0) {
                                itsBunch->calcBeamParametersInitial();
                            } else {
                                itsBunch->calcBeamParameters();
                            }

                            RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
                            RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

                            *gmsg << *itsBunch << endl;
                            RefPartR_suv_m = itsBunch->get_rmean();
                            RefPartP_suv_m = itsBunch->get_pmean();
                        }
                    } else {
                        throw OpalException("ParallelTTracker::execute()",
                                            "Unacceptable distribution type:\"" +
                                            d->getTypeofDistribution() + "\". Need to check the input file... ...");
                    }
                }
            }

            /// this is still in BoundaryGeometry
            size_t nbparts = dist->getNumberOfDarkCurrentParticles();
            double darkinwardmargin = dist->getDarkCurrentParticlesInwardMargin();
            double workFunction = dist->getWorkFunction();
            double fieldEnhancement = dist->getFieldEnhancement();
            size_t maxfnemission = dist->getMaxFNemissionPartPerTri();
            double fieldFNthreshold = dist->getFieldFNThreshold();
            double parameterFNA = dist->getFNParameterA();
            double parameterFNB = dist->getFNParameterB();
            double parameterFNY = dist->getFNParameterY();
            double parameterFNVYZe = dist->getFNParameterVYZero();
            double parameterFNVYSe = dist->getFNParameterVYSecond();


            secondaryFlg = dist->getSecondaryEmissionFlag();
            nEmissionMode = dist->getEmissionMode();
            bgf->setNEmissionMode(nEmissionMode);
            if(secondaryFlg) {
                if(secondaryFlg == 1) {
                    int BoundaryMatType = dist->getSurfMaterial();
                    bgf->setBoundaryMatType(BoundaryMatType);
                    if(Options::ppdebug) {
                        double vVThermal = dist->getvVThermal();//return thermal velocity of Maxwellian distribution of secondaries for benchmark
                        bgf->setvVThermal(vVThermal);
                        double ppVw = dist->getVw();
                        bgf->setVw(ppVw);

                    } else {
                        bgf->setvVThermal(1.0);
                        bgf->setVw(1.0);
                    }

                } else {
                    /*
                      parameters for Vaughan's secondary model
                    */
                    double vSeyZero = dist->getvSeyZero();// return sey_0 in Vaughan's model
                    double vEZero = dist->getvEZero();// return the energy related to sey_0 in Vaughan's model
                    double vSeyMax = dist->getvSeyMax();// return sey max in Vaughan's model
                    double vEmax = dist->getvEmax();// return Emax in Vaughan's model
                    double vKenergy = dist->getvKenergy();// return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
                    double vKtheta = dist->getvKtheta();// return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
                    double vVThermal = dist->getvVThermal();// return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
                    double ppVw = dist->getVw();
                    bgf->setVw(ppVw);
                    bgf->setvSeyZero(vSeyZero);
                    bgf->setvEZero(vEZero);
                    bgf->setvSeyMax(vSeyMax);
                    bgf->setvEmax(vEmax);
                    bgf->setvKenergy(vKenergy);
                    bgf->setvKtheta(vKtheta);
                    bgf->setvVThermal(vVThermal);
                }
            }
            if(nbparts != 0) {
                //fixme: maybe need to be called in each time step for modeling creating darkcurrent in each time step
                bgf->createParticlesOnSurface(nbparts, darkinwardmargin, itsOpalBeamline_m, *itsBunch);
                dist->create(*itsBunch, *bgf);
            }
            bgf->setWorkFunction(workFunction);
            bgf->setFieldEnhancement(fieldEnhancement);
            bgf->setMaxFN(maxfnemission);
            bgf->setFNTreshold(fieldFNthreshold);
            bgf->setFNParameterA(parameterFNA);
            bgf->setFNParameterB(parameterFNB);
            bgf->setFNParameterY(parameterFNY);
            bgf->setFNParameterVYZe(parameterFNVYZe);
            bgf->setFNParameterVYSe(parameterFNVYSe);
            totalParticles_i = itsBunch->getTotalNum();
            if(totalParticles_i > 0) {
                writePhaseSpace(0, 0);// dump the initial particles
            }
            itsDataSink->writeGeomToVtk(*bgf, string("data/testGeometry-00000.vtk"));
            //itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("vtk/PartlossZ-"));

            OpalData::getInstance()->setGlobalGeometry(bgf);

            RealVariable *maxnp = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MAXPARTSNUM"));
            if(maxnp) {
                maxNparts = static_cast<size_t>(maxnp->getReal());  // set upper limit of particle number in simulation
            }
            *gmsg << "Boundary geometry initialized " << endl;
            break;// only one boundary geometry allowed at present
        }
    }

    double minBinEmitted  = 10.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar) {
        minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
        // This variable controls the number of particles to be emitted before we use
        // the space charge solver.
        *gmsg << "MINBINEMITTED " << minBinEmitted << endl;
    }

    double minStepforReBin  = 200.0;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br) {
        minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        *gmsg << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

    double surfaceEmissionStop  = 1000.0;
    RealVariable *cr = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("SURFACEEMISSIONSTOP"));
    if(cr) {
        surfaceEmissionStop = cr->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        *gmsg << "SURFACEEMISSIONSTOP after " << surfaceEmissionStop << " seconds" <<  endl;
    }


    int repartFreq = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep) {
        repartFreq = static_cast<int>(rep->getReal());  // this variable controls the minimal number of steps until we repartition the particles
        *gmsg << "REPARTFREQ " << repartFreq << endl;
    }

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq = 1000000;

    size_t totalParticles_f = 0;
    bool wakestatus = false;
    bool surfacestatus = false;

    ofstream outfile("field.txt");
    for(; step < maxSteps_m; ++step) {
        global_EOL = true;  // check if any particle hasn't reached the end of the field from the last element
        bends = 0;
        hasWake = 0;
        wfSection = -1;
        hasSurfacePhysics = 0;
        sphysSection = -1;
        numberOfFieldEmittedParticles = 0;

        itsOpalBeamline_m.resetStatus();

        IpplTimings::startTimer(timeIntegrationTimer1_m);

        // reset E and B to Vector_t(0.0) for every step

        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);

        Nimpact_m = 0; // Initial parallel plate benchmark variable.
        SeyNum_m = 0; // Initial parallel plate benchmark variable.

        /// We do collision test for newly generated secondaries before integration in the first half step of each time step.
        /// This is because only secondary emission model yield non zero inital momenta. The initial momenta of field emitted particles are zero.
        //  If hit, we set itsBunch->R[i] to intersection points, else we do normal integration.
        /*==========================================Collision test for secondaries in 1st half time step==========================================*/
        if(bgf) {
            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();
            if(secondaryFlg) {
                for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
                    if(itsBunch->PType[i] == 3) { // only test newly generated secondaries
                        Vector_t intecoords = outr;
                        int triId = 0;
                        double Energy = 0.0;
                        int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], 0.5 * itsBunch->dt[i], itsBunch->PType[i], itsBunch->Q[i], intecoords, triId, Energy);

                        if(res == 0) {// if hit, set particle position to intersection points coordinates and scale the position;
                            itsBunch->R[i] = intecoords / vscaleFactor;

                            itsBunch->TriID[i] = triId;
                            // Fix me, is the local update necessary here?
                            itsBunch->X[i] /= vscaleFactor;
                            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                            itsBunch->X[i] *= vscaleFactor;

                        } else {// case: not hit, i.e.,res<0;
                            itsBunch->R[i] /= vscaleFactor;
                            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
                            // update local coordinate system for particle
                            itsBunch->X[i] /= vscaleFactor;
                            pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                            itsBunch->X[i] *= vscaleFactor;
                        }
                    } else {// The particles which are not the newly generated secondaries will do normal integration in the first half step.
                        itsBunch->R[i] /= vscaleFactor;
                        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
                        // update local coordinate system of particleInform &PartBunc
                        itsBunch->X[i] /= vscaleFactor;
                        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                        itsBunch->X[i] *= vscaleFactor;
                    }
                }
                /*===================================End of collision test for newly generated secondaries in first half step=============================*/

            } else {// Simulation without secondary emission module will do normal integration in the first half step.
                for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
                    // scale each particle with c*dt
                    itsBunch->R[i] /= vscaleFactor;
                    pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
                    // update local coordinate system of particleInform &PartBunc
                    itsBunch->X[i] /= vscaleFactor;
                    pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                    itsBunch->X[i] *= vscaleFactor;
                }
            }
        } else {// Simulation without BoundaryGeometry will do normal integration in the first half step.
            for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
                //scale each particle with c*dt
                itsBunch->R[i] /= vscaleFactor;
                pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
                // update local coordinate system of particleInform &PartBunc
                itsBunch->X[i] /= vscaleFactor;
                pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                itsBunch->X[i] *= vscaleFactor;
            }
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
            if(itsBunch->weHaveBins()) {
                // When we have energy bins.
                itsBunch->calcGammas();
                ParticleAttrib<double> Q_back = itsBunch->Q;
                for(int binNumber = 0; binNumber <= itsBunch->getLastemittedBin() &&
                    binNumber < itsBunch->getNumBins(); ++binNumber) {
                    itsBunch->setBinCharge(binNumber);
                    itsBunch->computeSelfFields(binNumber);
                    itsBunch->Q = Q_back;
                }
            } else {
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

        if((itsBunch->weHaveBins())) {
            // switch to TEmission
            if(!hasSwitchedToTEmission) {
                if(Options::fineEmission)
                    dt = itsBunch->getTSBin();
                else
                    dt = itsBunch->getTBin();
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Changing emission time step to: " << dt << endl;
                hasSwitchedToTEmission = true;
            }

            int ne = 0;
            if(Options::fineEmission)
                ne += itsBunch->emitParticlesNEW();
            else
                ne += itsBunch->emitParticles();

            if(Options::schottkyCorrection && !hasSwitchedBackToTTrack)
                applySchottkyCorrection(*itsBunch, ne, t, rescale_coeff);

            reduce(ne, ne, OpAddAssign());
            totalParticles_i += ne;

            //emission has finished, reset to TTrack
            if(itsBunch->getNumBins() == itsBunch->getLastemittedBin() &&
               !hasSwitchedBackToTTrack) {
                dt = dtTrack;
                itsBunch->setdT(dt);
                scaleFactor_m = dt * Physics::c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Emission done. Switching back to track timestep: " << dt << endl;
                hasSwitchedBackToTTrack = true;
            }

            if(step > minStepforReBin) {
                itsBunch->calcGammas();
                const double maxdE = abs(itsBunch->getMaxdEBins());
                if(maxdE < itsBunch->getRebinEnergy() && itsBunch->getNumBins()) {
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
                scaleFactor_m = dt * Physics::c;
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

            if(i == 1)
                outfile << itsBunch->Bf[i](0) << " " << itsBunch->Bf[i](1) << " " << itsBunch->Bf[i](2) << endl;


            itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

            // in case a particle is pushed behind the emission surface, delete the particle
            if(itsBunch->R[i](2) < 0)
                itsBunch->Bin[i] = -1;

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
                *gmsg << "============== START SURFACE PHYSICS CALCULATION ===Inform &PartBunc==========" << endl;
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
           Delete marked particles.
        */

        bool globPartOutOfBounds = (min(itsBunch->Bin) < 0);
        size_t ne = 0;
        if(globPartOutOfBounds) {
            ne = itsBunch->boundp_destroyT();
        }

        totalParticles_f = totalParticles_i - ne;
        if(ne > 0)
            *gmsg << "* Deleted " << ne << " particles, remaining " << totalParticles_f << " particles" << endl; //benchmark output

        if(hasWake > 0) {
            IpplTimings::startTimer(WakeFieldTimer_m);
            if(wf == NULL) {
                INFOMSG("no wakefunction attached" << endl);
            } else {
                wf->apply(*itsBunch);
            }
            IpplTimings::stopTimer(WakeFieldTimer_m);

        }

        if(hasSurfacePhysics > 0) {
            if(sphys == NULL) {
                INFOMSG("no surface physics attached" << endl);
            } else {
                sphys->apply(*itsBunch);
            }
        }
        if(bgf && secondaryFlg) {
            // Simulation with boundary geometry module which turns the secondary flag on  will
            // not kick those newly generated secondaries which have collided the boundary during the first half step integration.
            // These secondaries will be marked for deletion in the following main collision test part.
            int tmp = 0;
            kickParticles(pusher, tmp);
        } else {
            kickParticles(pusher);
        }

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
            RefPartP_zxy_m = dot(space_orientation_m, RefPartP_suv_m);
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

        /// After kick, we do collision test before integration in second half step with new momentum, if hit, then move collision particles to the position where collision occurs.
        if(bgf) {
            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();

            double dtime = 0.5 * itsBunch->getdT();
            for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {

                if(itsBunch->TriID[i] == 0) { // test all particles except those already have collided the boundary in the first half step.
                    Vector_t intecoords = outr;
                    int triId = 0;
                    double Energy = 0.0;
                    itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
                    int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId, Energy);

                    if(res == 0) {// if hit, set position to intersection points and do not need to scale here.
                        itsBunch->R[i] = intecoords;
                        itsBunch->TriID[i] = triId;
                        //  Fix me: is the updating of local coordinate system neccesary?
                        itsBunch->X[i] /= vscaleFactor;
                        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i], itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
                        itsBunch->X[i] *= vscaleFactor;
                        //reset time step if particle was emitted in the first half-step
                        //the particle is now in sync with the simulation timestep
                        itsBunch->dt[i] = itsBunch->getdT();

                    } else {//if no collision do normal push in the second half-step

                        itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);
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
                } else {
                    // the secondaries collide the boundary in the first half-step will not move and lie in the collision point and only scale the position.
                    // the collision particles in the first half-step will be scaled back.fix me the dt[i] is correct,
                    // do I need to update local coordinate for the collision particles in the first half-step?
                    itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);

                    //reset time step if particle was emitted in the first half-step
                    //the particle is now in sync with the simulation timestep
                    itsBunch->dt[i] = itsBunch->getdT();
                }
            }
        } else {
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
        }
        IpplTimings::stopTimer(timeIntegrationTimer2_m);

        ///Main collision test part.
        if(bgf) {
            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();
            /**
               Here we check if a particles is
               outside the domain, flag it for
               deletion and create secondaries
            */
            if(secondaryFlg == 1) { // entry for Furman-Pivi's secondary emission model
                // itsBunch->getLocalNum() will change immediately, so we need Inc_num to record the local particle number
                // before secondary emission, otherwise will be recursive generate secondaries and cause problem.
                size_t Inc_num = itsBunch->getLocalNum();

                double dtime = 0.5 * itsBunch->getdT();

                double seyNum = 0;

                for(size_t i = 0; i < Inc_num; i++) {

                    if(itsBunch->PType[i] == 3)
                        itsBunch->PType[i] = 2;// secondaries generated in last step will be set to be old secondaries.

                    if(itsBunch->TriID[i] == 0) {
                        // for primary bunch, primary dark current particles, old secondaries in previous time steps and newly
                        // generated secondaries which have no collision with boundary in both first and second half step, do main collision test and emit the secondaries.
                        Vector_t intecoords = outr;
                        int triId = 0;
                        double Energy = 0.0;
                        int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId, Energy);
                        if(res == 0) {
                            res += bgf->doBGphysics(intecoords, triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum);
                        }
                        if(res >= 0) {
                            itsBunch->Bin[i] = -1;
                            Nimpact_m++;
                            SeyNum_m += seyNum;
                        }
                    } else {
                        // Particles which collide the boundary in previous two tests will not do main collision test and directly call
                        // secondary emission module according to their energy and momentum before collision.
                        // Attention, these secondaries have not been kicked and are without new momentum.

                        double p_sq = dot(itsBunch->P[i], itsBunch->P[i]);
                        double Energy =  Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9;
                        int triId = itsBunch->TriID[i];

                        int res = bgf->doBGphysics(itsBunch->R[i], triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum);
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
                if(itsBunch->getT() < surfaceEmissionStop) {
                    numberOfFieldEmittedParticles += bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
                    itsBunch->boundp();
                    totalParticles_i = itsBunch->getTotalNum();
                } else
                    *gmsg << "* No field emission dT = " << t << endl;

            } else if(secondaryFlg != 0) {
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
                        double Energy = 0.0;
                        int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], dtime, itsBunch->PType[i], itsBunch->Q[i], intecoords, triId, Energy);

                        if(res == 0) {
                            res += bgf->doBGphysics(intecoords, triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum, para_null);
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
                        double p_sq = dot(itsBunch->P[i], itsBunch->P[i]);
                        double Energy =  Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9;
                        int triId = itsBunch->TriID[i];
                        //assert(dot(itsBunch->P[i], bgf->TriNormal_m[triId]) < 0);
                        int res = bgf->doBGphysics(itsBunch->R[i], triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum, para_null);

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
                if(itsBunch->getT() < surfaceEmissionStop) {
                    numberOfFieldEmittedParticles += bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
                    itsBunch->boundp();
                    totalParticles_i = itsBunch->getTotalNum();
                } else
                    *gmsg << "* No field emission dT = " << t << endl;

            } else {// the case without secondary emission, i.e., secondaryFlg==0
                for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
                    Vector_t intecoords = outr;
                    if(itsBunch->TriID[i] == 0) { // Particles which do not collide the boundary in collision test after kick
                        int triId = 0;
                        double Energy = 0.0;
                        int res = bgf->PartInside(itsBunch->R[i], itsBunch->P[i], itsBunch->getdT(), itsBunch->PType[i], itsBunch->Q[i], intecoords, triId, Energy);
                        if(res == 0) {
                            res += bgf->doBGphysics(intecoords, triId);
                        }
                        if(res >= 0) {
                            itsBunch->Bin[i] = -1;
                            Nimpact_m++;
                        }
                    } else {// Particles which collide the boundary in collision test after kick will not do main collision test and directly call doBGphysics function.
                        int triId = itsBunch->TriID[i];

                        //assert(dot(itsBunch->P[i], bgf->TriNormal_m[triId]) < 0);
                        int res = bgf->doBGphysics(intecoords, triId);

                        if(res >= 0) {
                            itsBunch->Bin[i] = -1;
                            Nimpact_m++;
                        }
                    }
                }

                /*========================
                  Now we do fieldemission
                  =========================== */
                if(itsBunch->getT() < surfaceEmissionStop)
                    numberOfFieldEmittedParticles += bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
                else
                    *gmsg << "* No field emission dT = " << t << endl;
                /*  if (itsBunch->getTotalNum()!= 0) {
                    itsBunch->boundp();
                    *gmsg<<"After boundp"<<endl;
                    }
                    totalParticles_i = itsBunch->getTotalNum();
                */
            }
        }

        if((totalParticles_f > minBinEmitted) || bgf) {
            itsBunch->boundp();

            // Until we do a boundp() itsBunch->getTotalNum() always returns zero particles.
            totalParticles_i = itsBunch->getTotalNum();
        }

        t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time

        itsBunch->setT(t);

        //IFF: cheap step dump regulation
        OPALTimer::Timer myt2;
        double sposRef = 0.0;
        if(totalParticles_f > 0) {
            sposRef = itsBunch->get_sPos();
            if(totalParticles_f <= minBinEmitted) {
                *gmsg << myt2.time() << " Step " << step << "; only " << totalParticles_f << " particles emitted; t= " << t
                      << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
            } else if(std::isnan(sposRef) || std::isinf(sposRef)) {
                *gmsg << myt2.time() << " Step " << step << "; there seems to be something wrong with the position of the bunch!" << endl;
            } else {
                *gmsg << myt2.time() << " Step " << step << " at " << sposRef << " [m] t= "
                      << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
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
                    writePhaseSpace(step, sposRef);
                    //          itsBunch->makHistograms();
                }

                //      itsBunch->makHistograms();


#ifdef DENSITY_W
                /*
                  Experimental write out density profile
                */

                stringstream filename_str;
                const int every = 10;
                bool print_criterion = (step + 1) % every == 0;

                if(print_criterion) {
                    *gmsg << "Experimental write out density profile ... ";
                    itsBunch->calcLineDensity();
                    itsBunch->getLineDensity(lineDensity_m);
                    const double h = itsBunch->get_hr()(2);
                    const unsigned int N = lineDensity_m.size();
                    if(Ippl::myNode() == 0) {
                        static unsigned int file_number = 0;
                        ++ file_number;
                        filename_str << "data/LineDensity-" << file_number << ".dat";
                        ofstream csr(filename_str.str().c_str());

                        csr << "# h= " << h << " (m)" << endl;
                        for(int i = 0; i < N; ++ i) {
                            csr << i *h << "\t"
                                << lineDensity_m[i] << "\t"
                                << lineDensity_m[i]*h << endl;
                        }
                        csr.close();
                    }
                    *gmsg << "** wrote " << filename_str.str() << endl;
                }
#endif
            }
            if(bgf) {
                reduce(SeyNum_m, SeyNum_m, OpAddAssign());
                reduce(Nimpact_m, Nimpact_m, OpAddAssign());
                itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("data/Partloss-"));



                long long ustep = step;
                itsDataSink->writeImpactStatistics(*itsBunch, ustep, Nimpact_m, SeyNum_m, numberOfFieldEmittedParticles, nEmissionMode, string("data/PartStatistics"));
                if(((Options::surfDumpFreq) > 0) && ((step % Options::surfDumpFreq) == 0)) {
                    itsDataSink->writeSurfaceInteraction(*itsBunch, ustep, *bgf, string("SurfaceInteraction"));
                }
            }
            /**
               Stop simulation if beyond zstop_m
            */
            if(sposRef > zstop_m) {
                maxSteps_m = step;
            }

            if(bgf) {   // If we are dealing with field emission and secondary emission, set upper limit of particle number in simulation to prevent memory overflow.
                if(totalParticles_i > maxNparts) {
                    maxSteps_m = step;
                }
            }
        } else {
            *gmsg << "Step " << step << " no emission yet "  << " t= " << t << " [s]" << endl;
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
        // free the memory allocated in monitors
        itsOpalBeamline_m.switchElementsOff(rmin(2) - margin, "Monitor");
    }
    bool doDump = true;
    writePhaseSpace((step + 1), itsBunch->get_sPos(), doDump);
    *gmsg << "Dump phase space of last step" << endl;
    OPALTimer::Timer myt3;
    OpalData::getInstance()->setLastStep(step);
    itsOpalBeamline_m.switchElementsOff();
    *gmsg << "done executing ParallelTTracker at " << myt3.time() << endl;
}
