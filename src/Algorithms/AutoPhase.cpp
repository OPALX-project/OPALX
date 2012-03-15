// ------------------------------------------------------------------------
// $RCSfile: AutoPhase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AutoPhase
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
#include <vector>

#include "Algorithms/AutoPhase.h"

#include "BeamlineGeometry/Euclid3D.h"
#include "BeamlineGeometry/PlanarArcGeometry.h"
#include "BeamlineGeometry/RBendGeometry.h"
#include "Beamlines/Beamline.h"
#include "Lines/Sequence.h"

#include "Utilities/NumToStr.h"
#include "Utilities/Timer.h"

#include "BeamlineCore/TravelingWaveRep.h"
#include "BeamlineCore/RFCavityRep.h"

#define PSdim 6

class PartData;

using namespace OPALTimer;
using Physics::c;

extern Inform *gmsg;
extern Inform *gmsg2all;

// Class AutoPhase
// ------------------------------------------------------------------------

AutoPhase::AutoPhase(const Beamline &beamline,
                     const PartData &reference,
                     bool revBeam,
                     bool revTrack):
    ParallelTTracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    itsOpalBeamline_m() {

}


AutoPhase::AutoPhase(const Beamline &beamline,
                     PartBunch &bunch,
                     DataSink &ds,
                     const PartData &reference,
                     bool revBeam,
                     bool revTrack,
                     int maxSTEPS,
                     double zstop,
                     long long actStep,
                     double actT,
                     int numRefinement):
    ParallelTTracker(beamline, reference, revBeam, revTrack),
    itsBeamline_m(beamline),
    maxSteps_m(maxSTEPS),
    zstop_m(zstop),
    step_m(actStep),
    actT_m(actT),
    numRefs_m(numRefinement) {
    //    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());
    itsBunch = &bunch;
    itsDataSink = &ds;
    scaleFactor_m = itsBunch->getdT() * c;

    timeIntegrationTimer1_m  = IpplTimings::getTimer("TIntegration1");
    timeIntegrationTimer2_m  = IpplTimings::getTimer("TIntegration2");
    timeFieldEvaluation_m  = IpplTimings::getTimer("Fieldeval");
}


AutoPhase::~AutoPhase() {
}


pair<FieldList::iterator , bool> AutoPhase::checkCavity(double s) {

    pair<FieldList::iterator , bool> res(cavities_m.begin(), false);
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit)
        if(((*fit).getStart() <= s) && (s <= (*fit).getEnd())) {
            res.first = fit;
            res.second = true;
            exit;
        }
    return res;
}


pair<FieldList::iterator, bool>  AutoPhase::doOneStep(BorisPusher pusher) {
    bool global_EOL = true;  //check if any particle hasn't reached the end of the field from the last element
    unsigned long bends = 0;

    double recpgamma, gamma;
    double t = itsBunch->getT();
    itsBunch->dt = itsBunch->getdT();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    Vector_t externalE, externalB;

    Vector_t rmin, rmax;

    itsOpalBeamline_m.resetStatus();

    //reset E and B to Vector_t(0.0) for every step

    itsBunch->Ef = Vector_t(0.0);
    itsBunch->Bf = Vector_t(0.0);

    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
        //scale each particle with c*dt
        itsBunch->R[i] /= vscaleFactor;

        pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

        // update local coordinate system of particle
        itsBunch->X[i] /= vscaleFactor;
        pusher.push(itsBunch->X[i], TransformTo(itsBunch->P[i],
                                                itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
        itsBunch->X[i] *= vscaleFactor;
    }

    itsBunch->calcBeamParameters();

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

    kickParticles(pusher);

    if(bends == 0) {
        updateReferenceParticle();
    } else {
        /* at least one of the elements bends the beam; until all particles
                  * have left the bending elements we track the reference particle
                  * as if it were a regular particle; from the moment when the reference
                  * particle has reached the bending field until it leaves
                  * it again we rotate the bunch about the position of the reference
                  * particle such that the momentum of the reference particle points
                  * in z direction
                  */
        RefPartP_suv_m = itsBunch->get_pmean();
        RefPartR_suv_m = itsBunch->get_rmean();
        recpgamma = 1. / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
        updateSpaceOrientation(false);

        // First update the momentum of the reference particle in zxy coordinate system, then update its position
        RefPartP_zxy_m(0) = SpaceOrientation_m[0] * RefPartP_suv_m(0) + SpaceOrientation_m[1] * RefPartP_suv_m(1) +
                            SpaceOrientation_m[2] * RefPartP_suv_m(2);
        RefPartP_zxy_m(1) = SpaceOrientation_m[3] * RefPartP_suv_m(0) + SpaceOrientation_m[4] * RefPartP_suv_m(1) +
                            SpaceOrientation_m[5] * RefPartP_suv_m(2);
        RefPartP_zxy_m(2) = SpaceOrientation_m[6] * RefPartP_suv_m(0) + SpaceOrientation_m[7] * RefPartP_suv_m(1) +
                            SpaceOrientation_m[8] * RefPartP_suv_m(2);
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
    itsBunch->get_bounds(rmin, rmax);

    // trigger the elements
    double margin = 3. * RefPartP_suv_m(2) * recpgamma;
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements((rmin(2) - margin)*scaleFactor_m, (rmax(2) + margin)*scaleFactor_m);

    // start particle loop part 2
    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
                                                itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])), itsBunch->getdT());
        itsBunch->X[i] *= vscaleFactor;
        //reset time step if particle was emitted in the first half-step
        //the particle is now in sync with the simulation timestep
        itsBunch->dt[i] = itsBunch->getdT();
    }

    t += itsBunch->getdT(); //t after a full global timestep with dT "synchronization point" for simulation time

    itsBunch->setT(t);

    return checkCavity(itsBunch->R[0](2));
}






void AutoPhase::execute() {
    double recpgamma, gamma;
    itsBunch->dt = itsBunch->getdT();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);



    BorisPusher pusher(itsReference);
    Vector_t rmin, rmax;

    bool global_EOL;
    unsigned long bends = 0; // flag which indicates whether any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch

    OPALTimer::Timer myt1;


    string fn = OPAL.getInputFn();

    int pos = fn.find(string("."), 0);
    fn.erase(pos, fn.size() - pos);

    string phasesFn = fn + string(".phases");

    ofstream maxphases;

    if(OPAL.hasBunchAllocated()) {
        maxphases.open(phasesFn.c_str(), ios::app);
    } else {
        itsBunch->setT(0.0);
        maxphases.open(phasesFn.c_str());
    }



    *gmsg << "Autophaseing start at: " << myt1.time() << " t= "
          << itsBunch->getT() << " zstop@ " << zstop_m << " [m]" << endl;

    *gmsg << "Executing AutoPhase, initial DT " << itsBunch->getdT()
          << " [s]; max integration steps " << maxSteps_m << ", step= " << step_m << ",\nR =  "
          << itsBunch->R[0] << " [m]" << endl;

    itsBeamline_m.accept(*this);
    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(*gmsg);

    // prepare lists
    cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
    travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
    cavities_m.merge(travelingwaves_m, OpalField::SortAsc);
    *gmsg << "Found the following cavities:" << endl;


    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit)
        *gmsg << (*fit).getElement()->getName()
              << " from " << (*fit).getStart() << " to "
              << (*fit).getEnd() << " (m)" << endl;

    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
        long &l = itsBunch->LastSection[i];
        l = -1;
        itsOpalBeamline_m.getSectionIndexAt(itsBunch->R[i], l);
        itsBunch->ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l),
                                             itsOpalBeamline_m.getSectionStart(l));
    }

    RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
    RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();

    /* Activate all elements which influence the particles when the simulation starts;
     * mark all elements which are already past.
     *
     * Increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
     * before they are allocated when increasing the timestep in the gun.
     */

    itsBunch->get_bounds(rmin, rmax);
    double margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    margin = 0.01 > margin ? 0.01 : margin;

    itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);
    for(; step_m < maxSteps_m; ++step_m) {
        itsBunch->setTrackStep(step_m);

        // let's do a drifting step to probe if the particle will reach element in next step
        Vector_t R_drift = itsBunch->R[0] + itsBunch->P[0] / sqrt(1.0 + dot(itsBunch->P[0], itsBunch->P[0])) * vscaleFactor;
        pair<FieldList::iterator, bool> res = checkCavity(R_drift(2));

        if(res.second) {
            double orig_phi = 0.0;
            Vector_t xi = itsBunch->R[0];
            Vector_t pi = itsBunch->P[0];
            double orig_t = itsBunch->getT();
            double Emax = 0.0;
            double Phimax = 0.0;
            double phi_ini = 0.0;
            double phi_end = 2. * Physics::pi;

            const double posErr  = (*res.first).getStart() - itsBunch->R[0](2);

            if((*res.first).getElement()->getType() == "TravelingWave") {
                orig_phi = static_cast<TravelingWave *>((*res.first).getElement())->getPhasem();
            } else {
                orig_phi = static_cast<RFCavity *>((*res.first).getElement())->getPhasem();
            }

            *gmsg << "Found " << (*res.first).getElement()->getName() << " at " << itsBunch->get_sPos() << " [m], step  " << step_m << " t= " << itsBunch->getT() << " [s],\n"
                  << "E= " << itsBunch->get_meanEnergy() << " [MeV], phi= " << orig_phi << endl;
            *gmsg << "Start phase scan ... " << endl;

            for(int j = 0; j < numRefs_m; ++ j) {  // numRefs_m == levels of refinement;
                double dphi = (phi_end - phi_ini) / 10.; //10 == number of sampling points
                // between phi_ini and phi_end. see in the opal user guide for an explanation
                // on the final accuracy.
                vector<PhaseEnT> e;

                for(double phi = phi_ini; phi < phi_end + dphi; phi += dphi) {

                    if((*res.first).getElement()->getType() == "TravelingWave") {
                        static_cast<TravelingWave *>((*res.first).getElement())->setPhasem(phi);
                    } else {
                        static_cast<RFCavity *>((*res.first).getElement())->setPhasem(phi);
                    }
                    doOneStep(pusher);  //this is the one step we replaced by a drifting step
                    while((itsBunch->R[0](2) < (*res.first).getEnd()) && (itsBunch->R[0](2) > (*res.first).getStart())) {
                        doOneStep(pusher);
                    }

                    // in case we made it until the end of the structure, save phase and momenta
                    if(itsBunch->R[0](2) >= (*res.first).getEnd())
                        e.push_back(PhaseEnT(phi, itsBunch->P[0]));

                    // reset initial conditions for next scan
                    itsBunch->R[0] = xi;
                    itsBunch->P[0] = pi;
                    itsBunch->setT(orig_t);
                }

                for(vector<PhaseEnT>::iterator eit = e.begin(); eit < e.end(); eit++) {
                    double E = ptoEMeV((*eit).second);
                    if(E > Emax) {
                        Emax = E;
                        Phimax = (*eit).first;
                    }
                }

                *gmsg << "found max energy at " << Phimax * 180. / Physics::pi << " degree (" << Emax << " MeV)" << endl;

                /*
                   at this stage we have done the scan and need to set the
                   correct phase for the final track before going on
                */
                phi_ini = Phimax - dphi;
                phi_end = Phimax + dphi;
            }

            if((*res.first).getElement()->getType() == "TravelingWave") {
                //static_cast<TravelingWave *>((*res.first).getElement())->setPhasem(Phimax - orig_phi);
                // do we want to track on-crest?? NO therefore we *HAVE* to subtract the phase that is
                // assigned to the element in the input file
                static_cast<TravelingWave *>((*res.first).getElement())->setPhasem(Phimax);
            } else {
                //static_cast<RFCavity *>((*res.first).getElement())->setPhasem(Phimax - orig_phi);
                // see above
                static_cast<RFCavity *>((*res.first).getElement())->setPhasem(Phimax);
            }

            itsBunch->R[0] = xi;
            itsBunch->P[0] = pi;
            maxphases << (*res.first).getElement()->getName() << "_phi= "  << Phimax << "; // E= " << Emax << " (MeV) Astra phases "
                      << (Phimax / Physics::pi * 180) + 90. << " (deg)" << endl;
            cavities_m.erase(res.first);
            itsBunch->calcBeamParameters();
        }
        doOneStep(pusher);
        double sposRef = itsBunch->get_sPos();

        if(sposRef > zstop_m)
            maxSteps_m = step_m;

        if(!(step_m % 1000)) {
            *gmsg << "step = " << step_m << ", spos = " << sposRef << " [m], t= " << itsBunch->getT() << " [s], "
                  << "E= " << itsBunch->get_meanEnergy() << " [MeV] " << endl;
        }

        if(!(step_m % 100)) {
            writePhaseSpace(step_m, sposRef);
            //            Vector_t FD(0.0);
            //itsDataSink->writeStatData(*itsBunch, &FD, sposRef, sposRef, sposRef);
        }

    }
    maxphases.close();
    OPALTimer::Timer myt3;
    OPAL.setLastStep(step_m);
    itsOpalBeamline_m.switchElementsOff();
    *gmsg << "done executing AutoPhase at " << myt3.time() << endl;
}



/**

Debug code


// lets save for debug purposes
stringstream pon;
pon << "phases_" << (*res.first).getElement()->getName() << "_stage-" << j << ".dat";
std::ofstream out(pon.str().c_str());

out << (*eit).first << " \t" << E * 1.0e-6 << endl;
out.close();

ofstream enlog("energy.dat");
enlog << itsBunch->R[0](2) << "\t" << ptoEMeV(itsBunch->P[0]) << endl;



*/
