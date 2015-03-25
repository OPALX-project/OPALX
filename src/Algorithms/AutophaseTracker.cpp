#include "Algorithms/AutophaseTracker.h"
#include "Algorithms/PartPusher.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Monitor.h"
#include "AbsBeamline/RFCavity.h"
#include "AbsBeamline/TravelingWave.h"
#include "Algorithms/PartData.h"
#include "Utilities/OpalOptions.h"

AutophaseTracker::AutophaseTracker(const Beamline &beamline, const PartData& data, const double &T0):
    DefaultVisitor(beamline, false, false),
    itsBunch_m(&data),
    itsPusher_m(data)
{
    itsBunch_m.setT(T0);
    visitBeamline(beamline);
}

AutophaseTracker::~AutophaseTracker()
{ }

// FieldList AutophaseTracker::executeAutoPhaseForSliceTracker() {
//     Inform msg("executeAutoPhaseForSliceTracker ");

//     double gPhaseSave;

//     gPhaseSave = OpalData::getInstance()->getGlobalPhaseShift();
//     OpalData::getInstance()->setGlobalPhaseShift(0.0);

//     // make sure that no monitor has overlap with two tracks
//     FieldList monitors = itsOpalBeamline_m.getElementByType("Monitor");
//     for(FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
//         double zbegin, zend;
//         it->getElement()->getDimensions(zbegin, zend);
//         if(zbegin < zStop_m.front() && zend >= zStop_m.front()) {
//             msg << "\033[0;31m"
//                 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
//                 << "% Removing '" << it->getElement()->getName() << "' since it resides in two tracks.   %\n"
//                 << "% Please adjust zstop or place your monitor at a different position to prevent this. %\n "
//                 << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
//                 << "\033[0m"
//                 << endl;
//             static_cast<Monitor *>(it->getElement().get())->moveBy(-zend - 0.001);
//             itsOpalBeamline_m.removeElement(it->getElement()->getName());
//         }
//     }
//     itsOpalBeamline_m.prepareSections();

//     cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
//     FieldList::iterator currentAPCavity = cavities_m.end();
//     FieldList travelingwaves = itsOpalBeamline_m.getElementByType("TravelingWave");
//     cavities_m.merge(travelingwaves, ClassicField::SortAsc);

//     int tag = 101;
//     int Parent = 0;
//     Vector_t iniR(0.0);
//     Vector_t iniP(0.0, 0.0, 1E-6);
//     PID_t id;
//     Ppos_t r, p, x;
//     ParticleAttrib<double> q, dt;
//     ParticleAttrib<int> bin;
//     ParticleAttrib<long> ls;
//     ParticleAttrib<short> ptype;

//     if(Ippl::myNode() == 0)
//         itsBunch_m.create(1);
//     itsBunch_m.update();

//     if(Ippl::myNode() == 0) {
//         itsBunch_m.R[0] = iniR;
//         itsBunch_m.P[0] = iniP;
//         itsBunch_m.Bin[0] = 0;
//         itsBunch_m.Q[0] = itsBunch_m.getChargePerParticle();
//         itsBunch_m.PType[0] = 0;
//         itsBunch_m.LastSection[0] = 0;

//         RefPartP_suv_m = iniP;
//         RefPartR_suv_m = iniR;
//     }

//     updateSpaceOrientation(false);
//     execute(zStop);
//     if(Ippl::myNode() == 0) {
//         RefPartP_suv_m = itsBunch_m.P[0];
//         // need to rebuild for updateAllRFElements
//         cavities_m = itsOpalBeamline_m.getElementByType("RFCavity");
//         travelingwaves = itsOpalBeamline_m.getElementByType("TravelingWave");
//         cavities_m.merge(travelingwaves, ClassicField::SortAsc);

//         // now send all max phases and names of the cavities to
//         // all the other nodes for updating.
//         Message *mess = new Message();
//         putMessage(*mess, OpalData::getInstance()->getNumberOfMaxPhases());

//         for(std::vector<MaxPhasesT>::iterator it = OpalData::getInstance()->getFirstMaxPhases(); it < OpalData::getInstance()->getLastMaxPhases(); it++) {
//             putMessage(*mess, (*it).first);
//             putMessage(*mess, (*it).second);
//         }
//         Ippl::Comm->broadcast_all(mess, tag);

//         delete mess;
//     } else {
//         // receive max phases and names and update the structures
//         int nData = 0;
//         Message *mess = Ippl::Comm->receive_block(Parent, tag);
//         getMessage(*mess, nData);
//         for(int i = 0; i < nData; i++) {
//             std::string elName;
//             double maxPhi;
//             getMessage(*mess, elName);
//             getMessage(*mess, maxPhi);
//             updateRFElement(elName, maxPhi);
//             OpalData::getInstance()->setMaxPhase(elName, maxPhi);
//         }

//         delete mess;
//     }

//     if(Ippl::myNode() == 0)
//         itsBunch_m.destroy(1, 0);
//     itsBunch_m.update();

//     OpalData::getInstance()->setGlobalPhaseShift(gPhaseSave);
//     return cavities_m;
// }


void AutophaseTracker::execute(const std::queue<double> &dtAllTracks,
                               const std::queue<double> &maxZ,
                               const std::queue<unsigned long long> &maxTrackSteps) {
    Inform msg("Autophasing ");

    // Vector_t rmin, rmax;

    // std::queue<double> dtAutoPhasing(dtTracks);
    // std::queue<double> maxSPosAutoPhasing(maxZ);
    // std::queue<unsigned long long> maxStepsAutoPhasing(maxTrackSteps);

    // maxSPosAutoPhasing.back() = itsOpalBeamline_m.calcBeamlineLenght();

    // const int numRefinements = Options::autoPhase;

    // size_t step = 0;
    // const int dtfraction = 2;

    // bends_m = 0;  //fixme: AP and bends will not work at the moment

    // for(unsigned int i = 0; i < itsBunch_m.getLocalNum(); ++i) {
    //     long &l = itsBunch_m.LastSection[i];
    //     l = -1;
    //     itsOpalBeamline_m.getSectionIndexAt(itsBunch_m.R[i], l);
    //     itsBunch_m.ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(l),
    //                                           itsOpalBeamline_m.getSectionStart(l));
    // }

    // RefPartR_suv_m = RefPartR_zxy_m = rmin = rmax = itsBunch_m.R[0];
    // RefPartP_suv_m = RefPartP_zxy_m = itsBunch_m.P[0];

    // /* Activate all elements which influence the particles when the simulation starts;
    //  * mark all elements which are already past.
    //  *
    //  * Increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
    //  * before they are allocated when increasing the timestep in the gun.
    //  */

    // double margin = 10. * RefPartP_suv_m(2) * scaleFactor / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));

    // margin = 0.01 > margin ? 0.01 : margin;

    // itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin, true);

    // double cavity_start = 0.0;
    // Component *cavity = NULL;

    // while (maxStepsAutoPhasing.size() > 0) {
    //     maxStepsAutoPhasing.front() = maxStepsAutoPhasing.front() * dtfraction + step;
    //     double newDt = dtAutoPhasing.front() / dtfraction;
    //     double scaleFactor = newDt * Physics::c;
    //     double vscaleFactor = Vector_t(scaleFactor);
    //     itsBunch_m.setdT(newDt);
    //     itsBunch_m.dt = newDt;

    //     msg << "\n"
    //         << "**********************************************\n"
    //         << "new set of dt, zstop and max number of time steps\n"
    //         << "**********************************************\n\n"
    //         << "start at t = " << itsBunch_m.getT() << " [s], "
    //         << "zstop at z = " << maxSPosAutoPhasing.front() << " [m],\n "
    //         << "dt = " << itsBunch_m.dt[0] << " [s], "
    //         << "step = " << step << ", "
    //         << "R =  " << itsBunch_m.R[0] << " [m]" << endl;

    //     for(; step < maxStepsAutoPhasing.front(); ++step) {


    //         adjustCavityPhase(cavity);


    //         IpplTimings::startTimer(timeIntegrationTimer1_m);
    //         double t = itsBunch_m.getT();

    //         // increase margin from 3.*c*dt to 10.*c*dt to prevent that fieldmaps are accessed
    //         // before they are allocated when increasing the timestep in the gun.

    //         switchElements(10.0);
    //         itsOpalBeamline_m.resetStatus();


    //         for(unsigned int i = 0; i < itsBunch_m.getLocalNum(); ++i) {
    //             itsBunch_m.R[i] /= Vector_t(Physics::c * itsBunch_m.getdT());
    //             itsPusher_m.push(itsBunch_m.R[i], itsBunch_m.P[i], itsBunch_m.dt[i]);
    //             itsBunch_m.R[i] *= Vector_t(Physics::c * itsBunch_m.getdT());
    //         }
    //         IpplTimings::stopTimer(timeIntegrationTimer1_m);

    //         itsBunch_m.Ef = Vector_t(0.0);
    //         itsBunch_m.Bf = Vector_t(0.0);

    //         IpplTimings::startTimer(timeFieldEvaluation_m);
    //         for(unsigned int i = 0; i < itsBunch_m.getLocalNum(); ++i) {
    //             long ls = itsBunch_m.LastSection[i];
    //             itsOpalBeamline_m.getSectionIndexAt(itsBunch_m.R[i], ls);
    //             if(ls != itsBunch_m.LastSection[i]) {
    //                 if(!itsOpalBeamline_m.section_is_glued_to(itsBunch_m.LastSection[i], ls)) {
    //                     itsBunch_m.ResetLocalCoordinateSystem(i, itsOpalBeamline_m.getOrientation(ls), itsOpalBeamline_m.getSectionStart(ls));
    //                 }
    //                 itsBunch_m.LastSection[i] = ls;
    //             }

    //             Vector_t externalE, externalB;
    //             const unsigned long rtv = itsOpalBeamline_m.getFieldAt(i, itsBunch_m.R[i], ls, itsBunch_m.getT() + itsBunch_m.dt[i] / 2., externalE, externalB);
    //             if (rtv)
    //                 ;
    //             itsBunch_m.Ef[i] += externalE;
    //             itsBunch_m.Bf[i] += externalB;
    //         }
    //         IpplTimings::stopTimer(timeFieldEvaluation_m);

    //         IpplTimings::startTimer(timeIntegrationTimer2_m);
    //         double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    //         RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor;

    //         for(unsigned int i = 0; i < itsBunch_m.getLocalNum(); ++i) {
    //             itsPusher_m.kick(itsBunch_m.R[i], itsBunch_m.P[i], itsBunch_m.Ef[i], itsBunch_m.Bf[i], itsBunch_m.dt[i]);
    //         }

    //         itsBunch_m.calcBeamParametersLight();

    //         for(unsigned int i = 0; i < itsBunch_m.getLocalNum(); ++i) {
    //             itsBunch_m.R[i] /= Vector_t(Physics::c * itsBunch_m.getdT());
    //             itsPusher_m.push(itsBunch_m.R[i], itsBunch_m.P[i], itsBunch_m.dt[i]);
    //             itsBunch_m.R[i] *= Vector_t(Physics::c * itsBunch_m.getdT());
    //         }
    //         t += itsBunch_m.getdT();
    //         itsBunch_m.setT(t);
    //         IpplTimings::stopTimer(timeIntegrationTimer2_m);

    //         double sposRef	= 0.0;

    //         if(!(step % 5000))	{

    //             if (itsBunch_m.getLocalNum() != 0)
    //                 sposRef = itsBunch_m.R[0](2);

    //             reduce(sposRef,sposRef,OpAddAssign());

    //             if(sposRef > maxSPosAutoPhasing.front())
    //                 maxStepsAutoPhasing.front() = step;


    //             INFOMSG("step = " << step << ", spos = " << sposRef << " [m], t= " << itsBunch_m.getT() << " [s], "
    //                     << "E= " << itsBunch_m.get_meanEnergy() << " [MeV] " << endl);
    //         }
    //     }

    //     dtAutoPhasing.pop();
    //     maxSPosAutoPhasing.pop();
    //     maxStepsAutoPhasing.pop();
    // }
}

void AutophaseTracker::adjustCavityPhase() {
    if (itsBunch_m.getLocalNum() == 0) return;

    // // let's do a drifting step to probe if the particle will reach element in next step
    // Vector_t R_drift = itsBunch_m.R[0] + itsBunch_m.P[0] / sqrt(1.0 + dot(itsBunch_m.P[0],
    //                                                                       itsBunch_m.P[0])) * vscaleFactor;
    // ClassicField comp = checkCavity(R_drift[2]);

    // if (comp == NULL) return;

    // Component *cavity = comp->getElement().get();
    // double cavityStart = comp->getStart();
    // double orig_phi = 0.0;
    // double Phimax = 0.0, Emax = 0.0;
    // double PhiAstra;

    // const double beta = sqrt(1. - 1 / (itsBunch_m.P[0](2) * itsBunch_m.P[0](2) + 1.));
    // const double tErr  = (cavityStart - itsBunch_m.R[0](2)) / (Physics::c * beta);

    // bool apVeto;

    // INFOMSG("Found " << cavity->getName()
    //         << " at " << itsBunch_m.R[0](2) << " [m], "
    //         << "step  " << step << ", "
    //         << "t= " << itsBunch_m.getT() << " [s],\n"
    //         << "E= " << getEnergyMeV(itsBunch_m.P[0]) << " [MeV]\n"
    //         << "start phase scan ... " << endl);


    // INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
    // if(cavity->getType() == "TravelingWave") {
    //     TravelingWave *element = static_cast<TravelingWave *>(cavity);
    //     orig_phi = element->getPhasem();
    //     apVeto = element->getAutophaseVeto();
    //     if(apVeto) {
    //         msg << " ----> APVETO -----> "
    //             << element->getName() <<  endl;
    //         Phimax = orig_phi;
    //     }
    //     INFOMSG(cavity->getName() << ", "
    //             << "start Ekin= " << getEnergyMeV(itsBunch_m.P[0]) << " MeV, "
    //             << "t= " << itsBunch_m.getT() << " s, "
    //             << "phi= " << orig_phi << ", " << endl;);

    //     if(!apVeto) {
    //         Phimax = element->getAutoPhaseEstimate(getEnergyMeV(itsBunch_m.P[0]),
    //                                                itsBunch_m.getT() + tErr,
    //                                                itsBunch_m.getQ(),
    //                                                itsBunch_m.getM() * 1e-6);
    //     }
    // } else {
    //     RFCavity *element = static_cast<RFCavity *>(cavity);
    //     orig_phi = element->getPhasem();
    //     apVeto = element->getAutophaseVeto();
    //     if(apVeto) {
    //         msg << " ----> APVETO -----> "
    //             << element->getName() << endl;
    //         Phimax = orig_phi;
    //     }
    //     INFOMSG(cavity->getName() << ", "
    //             << "start Ekin= " << getEnergyMeV(itsBunch_m.P[0]) << " MeV, "
    //             << "t= " << itsBunch_m.getT() << " s, "
    //             << "phi= " << orig_phi << ", " << endl;);

    //     if(!apVeto) {
    //         Phimax = element->getAutoPhaseEstimate(getEnergyMeV(itsBunch_m.P[0]),
    //                                                itsBunch_m.getT() + tErr,
    //                                                itsBunch_m.getQ(),
    //                                                itsBunch_m.getM() * 1e-6);
    //     }
    // }

    // double Phiini = Phimax;
    // double phi = Phiini;
    // double dphi = Physics::pi / 360.0;
    // int j = -1;

    // double E = APtrack(cavity, cavityStart, phi);
    // if(!apVeto) {
    //     msg << "Did APtrack with phi= " << phi << " result E= " << E << endl;
    //     do {
    //         j ++;
    //         Emax = E;
    //         Phiini = phi;
    //         phi -= dphi;
    //         E = APtrack(cavity, cavityStart, phi);
    //     } while(E > Emax);

    //     if(j == 0) {
    //         phi = Phiini;
    //         E = Emax;
    //         j = -1;
    //         do {
    //             j ++;
    //             Emax = E;
    //             Phiini = phi;
    //             phi += dphi;
    //             E = APtrack(cavity, cavityStart, phi);
    //         } while(E > Emax);
    //     }
    //     for(int refinementLevel = 0; refinementLevel < numRefinements; refinementLevel ++) {
    //         dphi /= 2.;
    //         phi = Phiini - dphi;
    //         E = APtrack(cavity, cavityStart, phi);
    //         if(E > Emax) {
    //             Phiini = phi;
    //             Emax = E;
    //         } else {
    //             phi = Phiini + dphi;
    //             E = APtrack(cavity, cavityStart, phi);
    //             if(E > Emax) {
    //                 Phiini = phi;
    //                 Emax = E;
    //             }
    //         }
    //     }
    //     Phimax = Phiini;
    //     phi = Phimax + orig_phi;

    //     INFOMSG("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n");
    // } else {
    //     msg << "Tracking with phi= " << orig_phi << " result E= " << E << endl;
    //     phi = orig_phi;
    //     Emax = E;
    //     E = APtrack(cavity, cavityStart, phi - Physics::pi);
    //     msg << "Tracking with phi= " << orig_phi - Physics::pi << " result E= " << E << endl;
    //     if (E > Emax) {
    //         phi = orig_phi - Physics::pi;
    //         Phimax = phi;
    //         Emax = E;
    //     }
    // }

    // if(cavity->getType() == "TravelingWave") {
    //     static_cast<TravelingWave *>(cavity)->updatePhasem(phi);
    // } else {
    //     static_cast<RFCavity *>(cavity)->updatePhasem(phi);
    // }

    // PhiAstra = (Phimax * Physics::rad2deg) + 90.0;
    // PhiAstra -= floor(PhiAstra / 360.) * 360.;

    // msg << cavity->getName() << "_phi= "  << Phimax << " rad / "
    //     << Phimax *Physics::rad2deg <<  " deg, AstraPhi= " << PhiAstra << " deg,\n"
    //     << "E= " << Emax << " (MeV), " << "phi_nom= " << orig_phi *Physics::rad2deg << endl;

    // OpalData::getInstance()->setMaxPhase(cavity->getName(), Phimax);
}

double AutophaseTracker::APtrack(Component *cavity, double cavityStart_pos, const double &phi) const {
    double beta = std::max(sqrt(1. - 1 / (itsBunch_m.P[0](2) * itsBunch_m.P[0](2) + 1.)), 0.0001);
    double tErr  = (cavityStart_pos - itsBunch_m.R[0](2)) / (Physics::c * beta);

    //    INFOMSG("beta = " << beta << " tErr = " << tErr << endl;);

    double finalMomentum = 0.0;
    if(cavity->getType() == "TravelingWave") {
        TravelingWave *tws = static_cast<TravelingWave *>(cavity);
        tws->updatePhasem(phi);
        std::pair<double, double> pe = tws->trackOnAxisParticle(itsBunch_m.P[0](2),
                                                                itsBunch_m.getT() + tErr,
                                                                itsBunch_m.dt[0],
                                                                itsBunch_m.getQ(),
                                                                itsBunch_m.getM() * 1e-6);
        finalMomentum = pe.first;
    } else {
        RFCavity *rfc = static_cast<RFCavity *>(cavity);
        rfc->updatePhasem(phi);

        std::pair<double, double> pe = rfc->trackOnAxisParticle(itsBunch_m.P[0](2),
                                                                itsBunch_m.getT() + tErr,
                                                                itsBunch_m.dt[0],
                                                                itsBunch_m.getQ(),
                                                                itsBunch_m.getM() * 1e-6);
        finalMomentum = pe.first;
    }
    double finalGamma = sqrt(1.0 + finalMomentum * finalMomentum);
    double finalKineticEnergy = (finalGamma - 1.0) * itsBunch_m.getM() * 1e-6;
    return finalKineticEnergy;
}


double AutophaseTracker::getGlobalPhaseShift() {

    if(Options::autoPhase > 0) {
        double gPhaseSave = OpalData::getInstance()->getGlobalPhaseShift();
        OpalData::getInstance()->setGlobalPhaseShift(0.0);
        return gPhaseSave;
    } else
        return 0;

}

ClassicField* AutophaseTracker::checkCavity(double s) {

    for(FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++ fit) {
        if((fit != currentAPCavity_m)
           && ((*fit).getStart() <= s) && (s <= (*fit).getEnd())) {
            currentAPCavity_m = fit;

            return &(*fit);
        }
    }

    return NULL;
}