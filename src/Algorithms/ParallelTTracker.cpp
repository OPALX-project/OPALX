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
#ifdef DBG_SYM
    string SfileName = OPAL.getInputFn();
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


pair<FieldList::iterator , bool> ParallelTTracker::checkCavity(double s) {

    pair<FieldList::iterator , bool> res(cavities_m.begin(), false);
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit)
        if(((*fit).getStart() <= s) && (s <= (*fit).getEnd())) {
            res.first = fit;
            res.second = true;
            exit;
        }
    return res;
}

pair<FieldList::iterator, bool> ParallelTTracker::doOneStep(BorisPusher pusher) {
    bool global_EOL = true;  //check if any particle hasn't reached the end of the field from the last element
    unsigned long bends = 0;

    double recpgamma, gamma;
    double t = itsBunch->getT();

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
    kickParticles(pusher);

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
    rmin = rmax = itsBunch->R[0];

    // trigger the elements
    double margin = 3. * RefPartP_suv_m(2) * recpgamma;
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements((rmin(2) - margin)*scaleFactor_m, (rmax(2) + margin)*scaleFactor_m, true);

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
	    m << static_cast<TravelingWave *>((*fit).getElement())->getPhasem()/Physics::pi*180.0 << endl;
	else 
	    m << static_cast<RFCavity *>((*fit).getElement())->getPhasem()/Physics::pi*180.0 << endl;
    }
    m << endl << endl;
}


void ParallelTTracker::updateRFElement(string elName,double maxPhi) {
    /**
       The maximum phase is added to the nominal phase of
       the element. This is done on all nodes except node 0 where
       the Autophase took place.
    */
    double phi  = 0.0;
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if ((*fit).getElement()->getName()==elName) {
	    if((*fit).getElement()->getType() == "TravelingWave") {
		phi  =  static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
		phi += maxPhi;
		static_cast<TravelingWave *>((*fit).getElement())->setPhasem(phi);
	    }
	    else {
		phi  = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
		phi += maxPhi;
		static_cast<RFCavity *>((*fit).getElement())->setPhasem(phi);
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
    const double RADDEG = 1.0/Physics::pi*180.0; 
    //   Inform m ("updateALLRFElements ",INFORM_ALL_NODES);
    *gmsg << "\n-------------------------------------------------------------------------------------\n";
    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if (fit != cavities_m.begin()) 
            *gmsg << "\n";
	if((*fit).getElement()->getType() == "TravelingWave") {
	    freq = static_cast<TravelingWave *>((*fit).getElement())->getFrequencym(); 
	    phi = static_cast<TravelingWave *>((*fit).getElement())->getPhasem();
	    *gmsg << (*fit).getElement()->getName() 
                  << ": phi_orig= phi_nom + phi_maxE= " << phi * RADDEG << " degree, "
                  << "global phase shift= " << -phiShift * freq * RADDEG << " degree\n";
	    phi -= (phiShift*freq);
	    static_cast<TravelingWave *>((*fit).getElement())->setPhasem(phi);
	}
	else {
	    freq = static_cast<RFCavity *>((*fit).getElement())->getFrequencym(); 
	    phi = static_cast<RFCavity *>((*fit).getElement())->getPhasem();
	    *gmsg << (*fit).getElement()->getName() 
                  << ": phi_orig= phi_nom + phi_maxE= " << phi*RADDEG << " degree, "
                  << "global phase shift= " << -phiShift * freq * RADDEG << " degree\n";
	    phi -= (phiShift*freq);
	    static_cast<RFCavity *>((*fit).getElement())->setPhasem(phi);
	}
    }
    *gmsg << "-------------------------------------------------------------------------------------\n"
          << endl;
}


void ParallelTTracker::executeAutoPhase(int numRefs, double zStop) {
    Inform msg("Autophasing ");

    const double RADDEG = 180.0/Physics::pi;

    double recpgamma, gamma;
    Vector_t rmin, rmax;
    bool global_EOL;

    size_t step = 0;
    int dtfraction = 2;
    itsBunch->dt = itsBunch->getdT() / dtfraction;         // need to fix this and make the factor 2 selectable 
    
    double scaleFactorSave = scaleFactor_m;
    scaleFactor_m = itsBunch->dt[0] * c;

    double tSave = itsBunch->getT();
    

    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    BorisPusher pusher(itsReference);

    unsigned long bends = 0; // flag which indicates whether any particle is within the influence of bending element.
    // if this is the case we track the reference particle as if it were a real particle,
    // otherwise the reference particle is defined as the centroid particle of the bunch

    msg << "\n"
        << "start at t= " << itsBunch->getT() << " [s], zstop at: " << zStop << " [m], Nplocal= " << itsBunch->getLocalNum() << "\n"
        << "initial DT " << itsBunch->dt[0] << " [s], step= " << step << ", R =  " << itsBunch->R[0] << " [m]" << endl;

    //    showCavities(m);

    for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
        Vector_t R_drift = itsBunch->R[0] + itsBunch->P[0] / sqrt(1.0 + dot(itsBunch->P[0], itsBunch->P[0])) * vscaleFactor;
     
	pair<FieldList::iterator, bool> res = checkCavity(R_drift(2));
	
        if(res.second) {
            double orig_phi = 0.0;
            double Phimax, Emax = 0.0;
            double PhiAstra;
            //////
            const double orig_t = itsBunch->getT();
            //////
            const double beta = sqrt(1. - 1 / (itsBunch->P[0](2)*itsBunch->P[0](2) + 1.));
            const double tErr  = ((*res.first).getStart() - itsBunch->R[0](2)) / (Physics::c * beta);

	    INFOMSG("Found " << (*res.first).getElement()->getName() 
                    << " at " << itsBunch->R[0](2) << " [m], " 
                    <<"step  " << step << ", "
                    << "t= " << itsBunch->getT() << " [s],\n" 
                    << "E= " << getEnergyMeV(itsBunch->P[0]) << " [MeV]\n"
                    << "start phase scan ... " << endl);
            INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
	    if((*res.first).getElement()->getType() == "TravelingWave") {
                orig_phi = static_cast<TravelingWave *>((*res.first).getElement())->getPhasem();
                INFOMSG((*res.first).getElement()->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);
                Phimax = static_cast<TravelingWave *>((*res.first).getElement())->getAutoPhaseEstimate(getEnergyMeV(itsBunch->P[0]), 
                                                                                                       itsBunch->getT() + tErr, 
                                                                                                       itsReference.getQ(), 
                                                                                                       itsReference.getM() * 1e-6);
           } else {
                orig_phi = static_cast<RFCavity *>((*res.first).getElement())->getPhasem();
                INFOMSG((*res.first).getElement()->getName() << ", "
                        << "start Ekin= " << getEnergyMeV(itsBunch->P[0]) << " MeV, "
                        << "t= " << itsBunch->getT() << " s, "
                        << "phi= " << orig_phi << ", " << endl;);
                Phimax = static_cast<RFCavity *>((*res.first).getElement())->getAutoPhaseEstimate(getEnergyMeV(itsBunch->P[0]), 
                                                                                                  itsBunch->getT() + tErr,
                                                                                                  itsReference.getQ(), 
                                                                                                  itsReference.getM() * 1e-6);
            }

            double Phiini = Phimax;
            double phi = Phiini;
            double dphi = Physics::pi / 360.0;
            int j = -1;

            double E = APtrack(res, phi);
            INFOMSG("do fine scan around effective max energy (" << E << " MeV)" << ", dphi= " << dphi << endl;);
            do {
                j ++;
                Emax = E;
                Phiini = phi;
                phi -= dphi;
                INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                E = APtrack(res, phi);
                if (E > Emax) {
		  INFOMSG(E - Emax << " MeV: accepted" << endl; );
                } else {
		  INFOMSG(E - Emax << " MeV: rejected" << endl;);
                }
            } while (E > Emax);
                
            if (j == 0) {
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
                    if (E > Emax) {
		      INFOMSG(E - Emax << " MeV: accepted" << endl;);
                    } else {
                      INFOMSG(E - Emax << " MeV: rejected" << endl;);
                    }
                } while (E > Emax);
            }
            for (int refinement_level = 0; refinement_level < numRefs; refinement_level ++) {
                dphi /= 2.;
                INFOMSG("refinement level: " << refinement_level + 1 << ", dphi= " << dphi << endl;);
                phi = Phiini - dphi;
                INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                E = APtrack(res, phi);
                if (E > Emax) {
		  INFOMSG(E - Emax << " MeV: accepted" << endl;);
                    Phiini = phi;
                    Emax = E;
                } else {
		  INFOMSG(E - Emax << " MeV: rejected" << endl;);
                    phi = Phiini + dphi;
                    INFOMSG("try phi= " << phi << " rad -> DEkin= ";);
                    E = APtrack(res, phi);
                    if (E > Emax) {
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


            if((*res.first).getElement()->getType() == "TravelingWave") {
                static_cast<TravelingWave *>((*res.first).getElement())->setPhasem(Phimax + orig_phi);
            } else {
                static_cast<RFCavity *>((*res.first).getElement())->setPhasem(Phimax + orig_phi);
            }

	    PhiAstra = (Phimax * RADDEG) + 90.0;
            PhiAstra -= floor(PhiAstra / 360.) * 360.;

            msg << (*res.first).getElement()->getName() << "_phi= "  << Phimax << " rad / " 
                << Phimax * RADDEG <<  " deg, AstraPhi= " << PhiAstra << " deg,\n"
                << "E= " << Emax << " (MeV), " << "phi_nom= " << orig_phi * RADDEG << endl;

	    maxPhases_m.push_back(MaxPhasesT((*res.first).getElement()->getName(),Phimax));
	    OPAL.setMaxPhase((*res.first).getElement()->getName(),Phimax);
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
    Inform m ("ParallelTTrackDebug ",INFORM_ALL_NODES);
    double recpgamma, gamma;
    double t = 0.0;
    double dt = itsBunch->getdT();
    double dtTrack = dt;
    double tEmission = itsBunch->getTEmission();
    const Vector_t vscaleFactor = Vector_t(scaleFactor_m);

    long long step = 0;
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

    unsigned long hasBGeom = 0;          // flag which indicates whether any particle is within the influence of a boundary geometry. Chuan: may be not in use?
    int BGSection = -1;                  // not in use
    BoundaryGeometry *bgf = NULL;
    int secondaryFlg = 0;
    size_t maxNparts = 100000000;        // upper limit of particle number when we do field emission and secondary emission simulation. Could be reset to another value in input file with MAXPARTSNUM.
    bool nEmissionMode = true;
    
    unsigned long hasSurfacePhysics = 0;
    int sphysSection = -1;
    SurfacePhysicsHandler *sphys;

    bool hasSwitchedToTEmission = false;
    bool hasSwitchedBackToTTrack = false;

    double gPhaseSave = 0.0;

    if(Options::autoPhase > 0) {
	gPhaseSave = OPAL.getGlobalPhaseShift();
	OPAL.setGlobalPhaseShift(0.0);
    }

    itsBeamline_m.accept(*this);
    // make sure that no monitor has overlap with two tracks
    FieldList monitors = itsOpalBeamline_m.getElementByType("Monitor");
    for (FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
        double zbegin, zend;
        it->getElement()->getDimensions(zbegin, zend);
        if (zbegin < zstop_m && zend >= zstop_m) {
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

    if ((!OPAL.inRestartRun()) && (Options::autoPhase > 0)) {
	int tag = 101;
	int Parent=0;
        Vector_t iniR(0.0);
        Vector_t iniP(0.0, 0.0, 1E-6);
        PID_t id;
        Ppos_t r, p, x;
        ParticleAttrib<double> q, dt;
        ParticleAttrib<int> bin;
        ParticleAttrib<long> ls;
        ParticleAttrib<short> ptype;

        size_t Nloc = itsBunch->getLocalNum();
        if (!OPAL.hasBunchAllocated() && Nloc > 0) {
            iniR = itsBunch->get_rmean();
            iniP = itsBunch->get_pmean();
            id.create(Nloc);    id = itsBunch->ID;
            r.create(Nloc);     r = itsBunch->R;
            p.create(Nloc);     p = itsBunch->P;
            x.create(Nloc);     x = itsBunch->X;
            q.create(Nloc);     q = itsBunch->Q;
            bin.create(Nloc);   bin = itsBunch->Bin;
            dt.create(Nloc);    dt = itsBunch->dt;
            ls.create(Nloc);    ls = itsBunch->LastSection;
            ptype.create(Nloc); ptype = itsBunch->PType;
            
            itsBunch->destroy(Nloc, 0);
            itsBunch->update();
        }
	if (Ippl::myNode() == 0) {
	    double zStop = itsOpalBeamline_m.calcBeamlineLenght();
	    if (!OPAL.hasBunchAllocated()) {
		itsBunch->create(1);
		itsBunch->R[0] = iniR;
		itsBunch->P[0] = iniP;
		itsBunch->Bin[0] = 0; 
		itsBunch->Q[0] = itsBunch->getChargePerParticle();
		itsBunch->PType[0] = 0;
		itsBunch->LastSection[0] = 0;
		executeAutoPhase(Options::autoPhase, zStop);
		itsBunch->destroy(1,0);	
		// need to rebuild for updateAllRFElements
		cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
		travelingwaves_m = itsOpalBeamline_m.getElementByType("TravelingWave");
		cavities_m.merge(travelingwaves_m, OpalField::SortAsc);
                
	    } else {
		// we are in a followup track and the phase information is 
		// already stored in the OPAL dictionary.
		for(vector<MaxPhasesT>::iterator it = OPAL.getFirstMaxPhases(); it < OPAL.getLastMaxPhases(); it++) {
		    updateRFElement((*it).first,(*it).second);
		    INFOMSG("In follow-up track use saved phases for -> name: " <<  (*it).first << " phi= " << (*it).second << " (rad)" << endl); 
		}
	    }
	    /* 
	       If needed we can write out 
	       the phase information here
 
	    string fn = OPAL.getInputFn();
	    int pos = fn.find(string("."), 0);
	    fn.erase(pos, fn.size() - pos);
	    string phasesFn = fn + string(".phases");

	    */

	    // now send all max phases and names of the cavities to 
	    // all the other nodes for updating.
	    Message *mess = new Message();
	    putMessage(*mess, OPAL.getNumberOfMaxPhases());

	    for(vector<MaxPhasesT>::iterator it = OPAL.getFirstMaxPhases(); it < OPAL.getLastMaxPhases(); it++) {
		putMessage(*mess, (*it).first);
		putMessage(*mess, (*it).second);
	    }
	    Ippl::Comm->broadcast_all(mess,tag);
	}
	else {
	    // receive max phases and names and update the structures 
	    int nData = 0;
	    Message *mess = Ippl::Comm->receive_block(Parent,tag);
	    getMessage(*mess,nData);
	    for(int i=0; i<nData; i++) {
		string elName;
		double maxPhi;
		getMessage(*mess,elName);
		getMessage(*mess,maxPhi);
		updateRFElement(elName,maxPhi);
		OPAL.setMaxPhase(elName,maxPhi);
	    }
	}

        if (!OPAL.hasBunchAllocated() && Nloc > 0) {
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
	OPAL.setGlobalPhaseShift(gPhaseSave);
    }
    else if (OPAL.inRestartRun() && Options::autoPhase > 0) { 
	/**
	   Restart and Autophase
	*/
	itsDataSink->retriveCavityInformation(OPAL.getRestartFileName());
	
	for(vector<MaxPhasesT>::iterator it = OPAL.getFirstMaxPhases(); it < OPAL.getLastMaxPhases(); it++) 
	    updateRFElement((*it).first,(*it).second);
    }
    
    if ((OPAL.getGlobalPhaseShift() > 0.0)  && (Options::autoPhase > 0))
	updateAllRFElements(OPAL.getGlobalPhaseShift());

    /** 
	save autophase information in order to skip 
	autophase in a restart run
	this is now going into teh file FN.phases, later should go the the h5
    */

    if ((!OPAL.inRestartRun()) && (Options::autoPhase > 0))
	itsDataSink->storeCavityInformation();
    

    size_t totalParticles_i = itsBunch->getTotalNum();
    *gmsg << "totalParticle_i= " << totalParticles_i << endl;
    OPALTimer::Timer myt1;

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
        if(itsBunch->getLocalNum() != 0)
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

    *gmsg << "Track start at: " << myt1.time() << ", t= " << itsBunch->getT() << "; zstop at: " << zstop_m << " [m]" << endl;

    if (!mpacflg_m) {
	if(itsBunch->doEmission()) {
	    emissionSteps = static_cast<int>(itsBunch->pbin_m->getNBins()) * gunSubTimeSteps;
	    *gmsg << "Do emission for " << tEmission << " [s] using " << itsBunch->pbin_m->getNBins() << " energy bins " << endl
		  << "Change dT from " <<  itsBunch->getdT() << " [s] to "
		  <<  itsBunch->getdT() << " [s] during emission " << endl;;
	}
	
	// set dt for all particles already in the simulation,
	// i.e. when doing a restarted simulation
	for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
	    itsBunch->dt[i] = itsBunch->getdT();
	}
    }
    
    *gmsg << "Executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s];\n"
          << "max integration steps " << maxSteps_m << ", step= " << step << ", Nplocal= " << itsBunch->getLocalNum() << endl;

    //    itsBeamline_m.accept(*this);
    //    itsOpalBeamline_m.prepareSections();
    itsOpalBeamline_m.print(*gmsg);
    double margin = 0.0;
    if (!mpacflg_m) {
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
	if(totalParticles_i == 0) {// fixme: maybe cause nonsense output if initialized momenta=0; Q: by Chuan.
	    itsBunch->calcBeamParametersInitial();
	} else {
	    itsBunch->calcBeamParameters();
	}
	
	RefPartR_suv_m = RefPartR_zxy_m = itsBunch->get_rmean();
	RefPartP_suv_m = RefPartP_zxy_m = itsBunch->get_pmean();
	if(!OPAL.hasBunchAllocated()) {
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
	m<<"rmax1= "<<rmax<<" rmin= "<<rmin<<endl;
	margin = 10. * RefPartP_suv_m(2) * scaleFactor_m / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
	margin = 0.01 > margin ? 0.01 : margin;
	itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin);
    }
    
    //bgf = itsOpalBeamline_m.getBoundaryGeometry(0);
    for (unsigned int i = 0; i<itsOpalBeamline_m.sections_m.size();i++) {
	
        bgf = itsOpalBeamline_m.getBoundaryGeometry(i);
	if(bgf) {
            Distribution *dist;
            Distribution *distrand;
            vector<string> distr_str = bgf->getDistributionArray();
            if(distr_str.size() == 0) {
                string distr = bgf->getDistribution();
                if(!distr.empty()) {
                    *gmsg << "* Find boundary geometry, start at: " << bgf->getS() << " (m) Distribution= " << bgf->getDistribution() << endl;
                    dist = Distribution::find(bgf->getDistribution());
                    *gmsg <<"* "<< *dist << endl;
                } else {
                    throw OpalException("ParallelTTracker::execute()",
					"No distribution attached to BoundaryGeometry. Please check the input file... ...");
		    
                }
            } else {
		*gmsg<< "************************************************************************************************* " << endl;
                *gmsg <<  "* Find boundary geometry, start at: " << bgf->getS()  <<" (m). "<< endl;
		*gmsg<< "* Attached more than one distribution: "<<endl;
		for(vector<string>::const_iterator dit = distr_str.begin(); dit != distr_str.end(); ++ dit) {
		  Distribution *d = Distribution::find(*dit);
		  *gmsg << "* Distribution: " << *dit << " distribution type: " << d->getTypeofDistribution() << endl;
		  *gmsg << "************************************************************************************************* "<< endl;
		  if(d->getTypeofDistribution() == "SURFACEEMISSION") {
		    dist = d;
		    *gmsg << *dist <<endl;

		  } else if(d->getTypeofDistribution() == "SURFACERANDCREATE") {
		    distrand = d;
		    *gmsg << *distrand <<endl;
		    size_t nbparts = distrand->getNumberOfDarkCurrentParticles();// here nbparts should be non zero as these particles will be the initialization of primary bunch.
		    double darkinwardmargin = distrand->getDarkCurrentParticlesInwardMargin();
		    double einitthreshold = distrand->getEInitThreshold();
		    bgf->setEInitThreshold(einitthreshold);// Caution: make sure that the elements have already been switched on before initializing particles in position where the electric field > einitthreshold. 
		    if (!mpacflg_m) {
			bgf->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
			distrand->createPriPart(itsBunch, *bgf );
			totalParticles_i = itsBunch->getTotalNum();
		    } else {// Multipacting flag set true. Generate primary particles.
			/* 
			   Activate all elements (switch on the field map of elements in multipacting) in multipacting simulation
			 */

			itsOpalBeamline_m.switchAllElements();

			bgf->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);// it is possible to generate initial particles according to E field, since all elements switched on before we create particles.
			distrand->createPriPart(itsBunch, *bgf );  //  For Parallel Plate benchmark, Vw should be defined in input file and will be invoked by getVw method in createPriPart(). For other multipacting simulation no need to define the Vw in SURFACERANDCREATE in input file.
			totalParticles_i = itsBunch->getTotalNum();
			itsBunch->calcBeamParameters();
			for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
					  d->getTypeofDistribution()+ "\". Need to check the input file... ...");
		    
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
                if(secondaryFlg==1) {
               
                    int BoundaryMatType = dist->getSurfMaterial();
                    bgf->setBoundaryMatType(BoundaryMatType);

                    if ( Options::ppdebug ) {
                
                        double vVThermal = dist->getvVThermal();//return thermal velocity of Maxwellian distribution of secondaries for benchmark
                        bgf->setvVThermal(vVThermal);
                        double ppVw = dist->getVw();
                        bgf->setVw(ppVw);
                   
                    } else {
                        bgf->setvVThermal(1.0);
                        bgf->setVw(1.0);
                    }
               
                }else {
                
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
           
                bgf->createParticlesOnSurface(nbparts, darkinwardmargin, itsOpalBeamline_m, *itsBunch);// fixme: maybe need to be called in each time step for modeling creating darkcurrent in each time step
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
	    if(totalParticles_i>0) {
		writePhaseSpace(0, 0);// dump the initial particles
	    }
            itsDataSink->writeGeomToVtk(*bgf, string("vtk/testGeometry-00000.vtk"));
            //itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("vtk/PartlossZ-"));

            
            OPAL.setGlobalGeometry(bgf);
	   
	    RealVariable *maxnp = dynamic_cast<RealVariable *>(OPAL.find("MAXPARTSNUM"));
	    if(maxnp) {
		maxNparts = static_cast<size_t>(maxnp->getReal());  // set upper limit of particle number in simulation
	    }
            *gmsg << "Boundary geometry initialized " << endl;
            break;// only one boundary geometry allowed at present
        }
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

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq = 1000000;

    size_t totalParticles_f = 0;
    bool wakestatus = false;
    bool surfacestatus = false;

  
    for(step; step < maxSteps_m; ++step) {
        global_EOL = true;  // check if any particle hasn't reached the end of the field from the last element
        bends = 0;
        hasWake = 0;
        wfSection = -1;
        hasSurfacePhysics = 0;
        sphysSection = -1;

        itsOpalBeamline_m.resetStatus();
       
        IpplTimings::startTimer(timeIntegrationTimer1_m);

        // reset E and B to Vector_t(0.0) for every step

        itsBunch->Ef = Vector_t(0.0);
        itsBunch->Bf = Vector_t(0.0);
       
        Nimpact_m=0;// Initial parallel plate benchmark variable.
        SeyNum_m=0;// Initial parallel plate benchmark variable.

        /// We do collision test for newly generated secondaries before integration in the first half step of each time step. 
	/// This is because only secondary emission model yield non zero inital momenta. The initial momenta of field emitted particles are zero.
	//  If hit, we set itsBunch->R[i] to intersection points, else we do normal integration.
        /*==========================================Collision test for secondaries in 1st half time step==========================================*/
        if(bgf) {

            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();

            if(secondaryFlg) {

                for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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

                for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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

            for(int i = 0; i < itsBunch->getLocalNum(); ++i) {
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
            if(totalParticles_i > 1000 && ((step % repartFreq) == 0)) {
                *gmsg << "*****************************************************************" << endl;
                *gmsg << "do repartition because of repartFreq" << endl;
                *gmsg << "*****************************************************************" << endl;
                IpplTimings::startTimer(BinRepartTimer_m);
                itsBunch->do_binaryRepart();
                IpplTimings::stopTimer(BinRepartTimer_m);
                Ippl::Comm->barrier();
                *gmsg << "*****************************************************************" << endl;
                *gmsg << "do repartition done" << endl;
                *gmsg << "*****************************************************************" << endl;
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
            if(itsBunch->getNumBins() == itsBunch->getLastemittedBin() && 
               !hasSwitchedBackToTTrack) {
                dt = dtTrack;
                itsBunch->setdT(dt);
                scaleFactor_m = dt * c;
                vscaleFactor = Vector_t(scaleFactor_m);
                *gmsg << "Emission done. Switching back to track timestep: " << dt << endl;
                hasSwitchedBackToTTrack = true;
                
		if(Options::correctOrbit == 1) {
		  // correct for centroid shift
		  itsBunch->boundp();
		  itsBunch->calcBeamParameters();
		  Vector_t cent = itsBunch->get_centroid();
		  Vector_t pmea = itsBunch->get_pmean();
		  
		  *gmsg << "Correct centroid after emission: Cold= " <<  cent;
		  cent(2)=0.0;
		  pmea(2)=0.0;
		  itsBunch->R -= cent;
		  itsBunch->P -= pmea;
		  itsBunch->boundp();
		  itsBunch->calcBeamParameters();
		  cent = itsBunch->get_centroid();
		  *gmsg << " ---->   centroid Cnew= " <<  cent << endl;
		}
	    }

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
            
            // in case a particle is pushed behind the emission surface, delete the particle
              
            if (itsBunch->R[i](2) < 0) 
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
	if (ne>0)
	  *gmsg<<"* Deleted " << ne << " particles, remaining " <<totalParticles_f<< " particles"<< endl; //benchmark output

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
        if(bgf && secondaryFlg) { // Simulation with boundary geometry module which turns the secondary flag on  will not kick those newly generated secondaries which have collided the boundary during the first half step integration. These secondaries will be marked for deletion in the following main collision test part. Chuan
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
      
        /// After kick, we do collision test before integration in second half step with new momentum, if hit, then move collision particles to the position where collision occurs.
        if(bgf) {

            const Vector_t outr = bgf->getmaxcoords() + bgf->gethr();


            double dtime = 0.5 * itsBunch->getdT();
            for(int i = 0; i < itsBunch->getLocalNum(); ++i) {

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
                } else {//the secondaries collide the boundary in the first half-step will not move and lie in the collision point and only scale the position.

                    itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i], Physics::c * itsBunch->dt[i]);// the collision particles in the first half-step will be scaled back.fix me the dt[i] is correct, do I need to update local coordinate for the collision particles in the first half-step?


                    //reset time step if particle was emitted in the first half-step
                    //the particle is now in sync with the simulation timestep
                    itsBunch->dt[i] = itsBunch->getdT();
                }

            }


            /*==================================================================================================================*/

        } else {// start normal particle loop part 2 for simulation without boundary geometry.
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
            if(secondaryFlg==1) {// entry for Furman-Pivi's secondary emission model

                size_t Inc_num = itsBunch->getLocalNum();// itsBunch->getLocalNum() will change immediately, so we need Inc_num to record the local particle number before secondary emission, otherwise will be recursive generate secondaries and cause problem.

                double dtime = 0.5 * itsBunch->getdT();
               
                double seyNum=0;
                
                for(size_t i = 0; i < Inc_num; i++) {

                    if(itsBunch->PType[i] == 3)
                        itsBunch->PType[i] = 2;// secondaries generated in last step will be set to be old secondaries.

                    if(itsBunch->TriID[i] == 0) { // for primary bunch, primary dark current particles, old secondaries in previous time steps and newly generated secondaries which have no collision with boundary in both first and second half step, do main collision test and emit the secondaries.
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
                            SeyNum_m+=seyNum;
			}
                    } else {// Particles which collide the boundary in previous two tests will not do main collision test and directly call secondary emission module according to their energy and momentum before collision. Attention, these secondaries have not been kicked and are without new momentum.
                        
                        double p_sq = dot(itsBunch->P[i], itsBunch->P[i]);
                        double Energy =  Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9;
                        int triId = itsBunch->TriID[i];
                                               
                        int res = bgf->doBGphysics(itsBunch->R[i], triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum);
                        if(res >= 0) {
                            itsBunch->Bin[i] = -1;
                            Nimpact_m++;
                            SeyNum_m+=seyNum;
			}
                    }
                }

                /*===========================
                  Now we do fieldemission
                ============================== */
                bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
                itsBunch->boundp();
                totalParticles_i = itsBunch->getTotalNum();
            } else if(secondaryFlg!=0) {// entry for Vaughan's secondary emission model

		const int para_null = 0;// dummy parameter for overloading the Vaughan's version of BoundaryGeometry::doBGphysics();

                size_t Inc_num = itsBunch->getLocalNum();// itsBunch->getLocalNum() will change immediately, so we need Inc_num to record the local particle number before secondary emission, otherwise will be recursive generate secondaries and cause problem.

                double dtime = 0.5 * itsBunch->getdT();
               
                double seyNum=0;
                
                for(size_t i = 0; i < Inc_num; i++) {

                    if(itsBunch->PType[i] == 3)
                        itsBunch->PType[i] = 2;// secondaries generated in last step will be set to be old secondaries.

                    if(itsBunch->TriID[i] == 0) { // for primary bunch, primary dark current particles, old secondaries in previous time steps and newly generated secondaries which have no collision with boundary in both first and second half step, do main collision test and emit the secondaries.
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
                            SeyNum_m+=seyNum;
                        }
                    } else {//Particles which collide the boundary in previous two tests will not do main collision test and directly call secondary emission module according to their energy and momentum before collision. Attention, these secondaries have not been kicked and are without new momentum.
                        double p_sq = dot(itsBunch->P[i], itsBunch->P[i]);
                        double Energy =  Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9;
                        int triId = itsBunch->TriID[i];
                        //assert(dot(itsBunch->P[i], bgf->TriNormal_m[triId]) < 0);
                        int res = bgf->doBGphysics(itsBunch->R[i], triId, Energy, itsBunch->Q[i], itsBunch->P[i], itsBunch, seyNum, para_null);

                        if(res >= 0) {
                            itsBunch->Bin[i] = -1;
                            Nimpact_m++;
                            SeyNum_m+=seyNum;
                        }
                    }
                }

                /*===========================
                  Now we do fieldemission
                ============================== */
                bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
                itsBunch->boundp();
                totalParticles_i = itsBunch->getTotalNum();

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
                        double p_sq = dot(itsBunch->P[i], itsBunch->P[i]);
                        double Energy =  Physics::m_e * (sqrt(1.0 + p_sq) - 1.0) * 1.0e9;
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
		bgf->doFNemission(itsOpalBeamline_m, itsBunch, t);
		/*	if (itsBunch->getTotalNum()!= 0) {
		    itsBunch->boundp();
		    *gmsg<<"After boundp"<<endl;
		}
                totalParticles_i = itsBunch->getTotalNum();*/
            }

        }
       	     
        if((totalParticles_f > minBinEmitted) || bgf)
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
                *gmsg << myt2.time() << " Step " << step << "; only " << totalParticles_f << " particles emitted; t= " << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
            } else if(isnan(sposRef) || isinf(sposRef)) {
                *gmsg << myt2.time() << " Step " << step << "; there seems to be something wrong with the position of the bunch!" << endl;
            } else {
                *gmsg << myt2.time() << " Step " << step << " at " << sposRef << " [m] t= "
                      << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl;
                if(step % Options::psDumpFreq == 0 || step % Options::statDumpFreq == 0) {
                    size_t nLoc = itsBunch->getLocalNum();
		    reduce(nLoc, nLoc, OpMultipplyAssign());
		    if( (nLoc == 0) || (step % repartFreq == 0) ) {
                        *gmsg << "*****************************************************************" << endl;
                        *gmsg << "do repartition because of zero particles or repartition frequency" << endl;
			IpplTimings::startTimer(BinRepartTimer_m);
                        itsBunch->do_binaryRepart();
                        IpplTimings::stopTimer(BinRepartTimer_m);
                        Ippl::Comm->barrier();
                        *gmsg << "done" << endl;
                        nLoc = itsBunch->getLocalNum();
                        reduce(nLoc, nLoc, OpMultipplyAssign());

                        if(nLoc == 0) {
                            *gmsg << "*****************************************************************" << endl;

                            *gmsg << "Zero Particles on a node, after repartitioning.           " << endl;
                            *gmsg << "Cannot go further with the simulation, please restart     " << endl;
                            *gmsg << "with less cores from the last h5 dump (-restart -1)       " << endl;
                            *gmsg << "*****************************************************************" << endl;
                            maxSteps_m = step;
                        }
                        *gmsg << "*****************************************************************" << endl;
                    }
		    
                    writePhaseSpace(step, sposRef);

                    
                }
            }
	    if(bgf) {
		reduce(SeyNum_m, SeyNum_m, OpAddAssign());
		reduce(Nimpact_m, Nimpact_m, OpAddAssign());
		//itsDataSink->writePartlossZASCII(*itsBunch, *bgf, string("vtk/Partloss-"));
		itsDataSink->writeImpactStatistics(*itsBunch, step, Nimpact_m, SeyNum_m, nEmissionMode, string("Part_statistics-"));
		if( ((Options::surfDumpFreq)>0) && ((step % Options::surfDumpFreq)==0) ) {
		    itsDataSink->writeSurfaceInteraction(*itsBunch, step, *bgf, string("SurfaceInteraction"));
		}
	    }
            /**
               Stop simulation if beyond zstop_m
            */
            if(sposRef > zstop_m) {
                maxSteps_m = step;
            }

	    if ( bgf ) {// If we are dealing with field emission and secondary emission, set upper limit of particle number in simulation to prevent memory overflow.
		

		if(totalParticles_i>maxNparts) {
		    maxSteps_m = step;
		}
	    }
        } else {
            *gmsg << "Step " << step << " no emission yet "  << " t= " << t << " [s]" << endl;
        }

        if(step > emissionSteps) {
            //TEST the non-len reduce: reduce(&global_EOL, &global_EOL, OpBitwiseOrAssign());
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


	if(hasSwitchedBackToTTrack && (Options::correctOrbit > 1) && (step % Options::correctOrbit == 0)) {
	  // correct for centroid shift
	  Vector_t cent = itsBunch->get_centroid();
	  Vector_t pmea = itsBunch->get_pmean();
	  *gmsg << "Correct centroid after " << Options::correctOrbit  << " timestepe: Cold= " <<  cent;
	  cent(2)=0.0;
	  pmea(2)=0.0;
	  itsBunch->R -= cent;
	  itsBunch->P -= pmea;
	  itsBunch->boundp();
	  itsBunch->calcBeamParameters();
	  cent = itsBunch->get_centroid();
	  *gmsg << " ---->   centroid Cnew= " <<  cent << endl;
	}
    }

    OPALTimer::Timer myt3;
    OPAL.setLastStep(step);

    itsOpalBeamline_m.switchElementsOff();

    *gmsg << "done executing ParallelTTracker at " << myt3.time() << endl;
   
}
