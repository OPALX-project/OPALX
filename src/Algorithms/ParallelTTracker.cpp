//
// Class ParallelTTracker
//   OPAL-T tracker.
//   The visitor class for tracking particles with time as independent
//   variable.
//
// Copyright (c) 200x - 2014, Christof Kraus, Paul Scherrer Institut, Villigen PSI, Switzerland
//               2015 - 2016, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020, Christof Metzger-Kraus
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Algorithms/ParallelTTracker.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "AbsBeamline/Monitor.h"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/OrbitThreader.h"
#include "Algorithms/CavityAutophaser.h"
#include "BasicActions/Option.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Distribution/Distribution.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Units.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/BoundingBox.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "ValueDefinitions/RealVariable.h"

extern Inform* gmsg;

class PartData;

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack):
    Tracker(beamline, reference, revBeam, revTrack),
    itsDataSink_m(nullptr),
    itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
    globalEOL_m(false),
    wakeStatus_m(false),
    wakeFunction_m(nullptr),
    pathLength_m(0.0),
    zstart_m(0.0),
    dtCurrentTrack_m(0.0),
    minStepforReBin_m(-1),
    repartFreq_m(-1),
    emissionSteps_m(std::numeric_limits<unsigned int>::max()),
    numParticlesInSimulation_m(0),
    timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
    timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
    fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
    BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart"))
{}

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunchBase<double, 3> *bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   const std::vector<unsigned long long> &maxSteps,
                                   double zstart,
                                   const std::vector<double> &zstop,
                                   const std::vector<double> &dt):
    Tracker(beamline, bunch, reference, revBeam, revTrack),
    itsDataSink_m(&ds),
    itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
    globalEOL_m(false),
    wakeStatus_m(false),
    wakeFunction_m(nullptr),
    pathLength_m(0.0),
    zstart_m(zstart),
    dtCurrentTrack_m(0.0),
    minStepforReBin_m(-1),
    repartFreq_m(-1),
    emissionSteps_m(std::numeric_limits<unsigned int>::max()),
    numParticlesInSimulation_m(0),
    timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
    timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
    fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
    BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart"))
{
    for (unsigned int i = 0; i < zstop.size(); ++ i) {
        stepSizes_m.push_back(dt[i], zstop[i], maxSteps[i]);
    }

    stepSizes_m.sortAscendingZStop();
    stepSizes_m.resetIterator();
}

ParallelTTracker::~ParallelTTracker() {}

void ParallelTTracker::visitBeamline(const Beamline &bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    if (fbl->getRelativeFlag()) {
        OpalBeamline stash(fbl->getOrigin3D(), fbl->getInitialDirection());
        stash.swap(itsOpalBeamline_m);
        fbl->iterate(*this, false);
        itsOpalBeamline_m.prepareSections();
        itsOpalBeamline_m.compute3DLattice();
        stash.merge(itsOpalBeamline_m);
        stash.swap(itsOpalBeamline_m);
    } else {
        fbl->iterate(*this, false);
    }
}

void ParallelTTracker::updateRFElement(std::string elName, double maxPhase) {
    FieldList cavities = itsOpalBeamline_m.getElementByType(ElementType::RFCAVITY);
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType(ElementType::TRAVELINGWAVE);
    cavities.insert(cavities.end(), travelingwaves.begin(), travelingwaves.end());

    for (FieldList::iterator fit = cavities.begin(); fit != cavities.end(); ++ fit) {
        if ((*fit).getElement()->getName() == elName) {

            RFCavity *element = static_cast<RFCavity *>((*fit).getElement().get());

            element->setPhasem(maxPhase);
            element->setAutophaseVeto();

            INFOMSG("Restored cavity phase from the h5 file. Name: " << element->getName() << ", phase: " << maxPhase << " rad" << endl);
            return;
        }
    }
}

void ParallelTTracker::saveCavityPhases() {
    itsDataSink_m->storeCavityInformation();
}

void ParallelTTracker::restoreCavityPhases() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() ||
        OpalData::getInstance()->inRestartRun()) {
        iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++ it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}

void ParallelTTracker::execute() {
    Inform msg("ParallelTTracker ", *gmsg);
    OpalData::getInstance()->setInPrepState(true);

    BorisPusher pusher(itsReference);
    const double globalTimeShift = itsBunch_m->weHaveEnergyBins()? OpalData::getInstance()->getGlobalPhaseShift(): 0.0;
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    // the time step needs to be positive in the setup
    itsBunch_m->setdT(std::abs(itsBunch_m->getdT()));
    dtCurrentTrack_m = itsBunch_m->getdT();

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::APPEND);
    }

    prepareSections();

    double minTimeStep = stepSizes_m.getMinTimeStep();

    itsOpalBeamline_m.activateElements();

    if ( OpalData::getInstance()->hasPriorTrack() ||
         OpalData::getInstance()->inRestartRun()) {

        pathLength_m = itsBunch_m->get_sPos();
        zstart_m = pathLength_m;

        restoreCavityPhases();
    } else {

        double momentum = euclidean_norm(itsBunch_m->get_pmean_Distribution());
        CoordinateSystemTrafo beamlineToLab = itsOpalBeamline_m.getCSTrafoLab2Local().inverted();
        itsBunch_m->toLabTrafo_m = beamlineToLab;

        itsBunch_m->RefPartR_m = beamlineToLab.transformTo(Vector_t(0.0));
        itsBunch_m->RefPartP_m = beamlineToLab.rotateTo(momentum * Vector_t(0, 0, 1));

        if (itsBunch_m->getTotalNum() > 0) {

            if (zstart_m > pathLength_m) {
                findStartPosition(pusher);
            }

            itsBunch_m->set_sPos(pathLength_m);
        }
    }

    stepSizes_m.advanceToPos(zstart_m);
    if (back_track) {
        itsBunch_m->setdT(-std::abs(itsBunch_m->getdT()));
        stepSizes_m.reverseDirection();
        if (pathLength_m < stepSizes_m.getZStop()) {
            ++ stepSizes_m;
        }
    }

    Vector_t rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNum() > 0) {
        itsBunch_m->get_bounds(rmin, rmax);
    }

    OrbitThreader oth(itsReference,
                      itsBunch_m->RefPartR_m,
                      itsBunch_m->RefPartP_m,
                      pathLength_m,
                      -rmin(2),
                      itsBunch_m->getT(),
                      (back_track? -minTimeStep: minTimeStep),
                      stepSizes_m,
                      itsOpalBeamline_m);

    oth.execute();
    BoundingBox globalBoundingBox = oth.getBoundingBox();

    numParticlesInSimulation_m = itsBunch_m->getTotalNum();

    setTime();

    double time = itsBunch_m->getT() - globalTimeShift;
    itsBunch_m->setT(time);

    *gmsg << level1 << *itsBunch_m << endl;

    unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << "* Track start at: " << myt1.time() << ", t= " << Util::getTimeString(time) << "; "
          << "zstart at: " << Util::getLengthString(pathLength_m)
          << endl;

    *gmsg << "* Executing ParallelTTracker\n"
          << "* Initial dt = " << Util::getTimeString(itsBunch_m->getdT()) << "\n"
          << "* Max integration steps = " << stepSizes_m.getMaxSteps()
          << ", next step = " << step << endl << endl;

    setOptionalVariables();

    globalEOL_m = false;
    wakeStatus_m = false;
    deletedParticles_m = false;
    OpalData::getInstance()->setInPrepState(false);
    while (!stepSizes_m.reachedEnd()) {
        unsigned long long trackSteps = stepSizes_m.getNumSteps() + step;
        dtCurrentTrack_m = stepSizes_m.getdT();
        changeDT(back_track);

        for (; step < trackSteps; ++ step) {
            Vector_t rmin(0.0), rmax(0.0);
            if (itsBunch_m->getTotalNum() > 0) {
                itsBunch_m->get_bounds(rmin, rmax);
            }

            timeIntegration1(pusher);

            itsBunch_m->Ef = Vector_t(0.0);
            itsBunch_m->Bf = Vector_t(0.0);

            computeSpaceChargeFields(step);

            selectDT(back_track);

            selectDT(back_track);

            computeExternalFields(oth);

            timeIntegration2(pusher);

            itsBunch_m->incrementT();

            if (itsBunch_m->getT() > 0.0 ||
                itsBunch_m->getdT() < 0.0) {
                updateReference(pusher);
            }

            if (deletedParticles_m) {
                deletedParticles_m = false;
            }
            itsBunch_m->set_sPos(pathLength_m);

            if (hasEndOfLineReached(globalBoundingBox)) break;

            bool const psDump = ((itsBunch_m->getGlobalTrackStep() % Options::psDumpFreq) + 1 == Options::psDumpFreq);
            bool const statDump = ((itsBunch_m->getGlobalTrackStep() % Options::statDumpFreq) + 1 == Options::statDumpFreq);
            dumpStats(step, psDump, statDump);

            itsBunch_m->incTrackSteps();

            double beta = euclidean_norm(itsBunch_m->RefPartP_m / Util::getGamma(itsBunch_m->RefPartP_m));
            double driftPerTimeStep = std::abs(itsBunch_m->getdT()) * Physics::c * beta;
            if (std::abs(stepSizes_m.getZStop() - pathLength_m) < 0.5 * driftPerTimeStep) {
                break;
            }
        }

        if (globalEOL_m)
            break;

        ++ stepSizes_m;
    }

    itsBunch_m->set_sPos(pathLength_m);

    numParticlesInSimulation_m = itsBunch_m->getTotalNum();

    bool const psDump = (((itsBunch_m->getGlobalTrackStep() - 1) % Options::psDumpFreq) + 1 != Options::psDumpFreq);
    bool const statDump = (((itsBunch_m->getGlobalTrackStep() - 1) % Options::statDumpFreq) + 1 != Options::statDumpFreq);
    writePhaseSpace((step + 1), psDump, statDump);

    msg << level2 << "Dump phase space of last step" << endl;

    itsOpalBeamline_m.switchElementsOff();

    OPALTimer::Timer myt3;
    *gmsg << endl << "* Done executing ParallelTTracker at "
          << myt3.time() << endl << endl;

    Monitor::writeStatistics();

    OpalData::getInstance()->setPriorTrack();
}

void ParallelTTracker::prepareSections() {

    itsBeamline_m.accept(*this);
    
    itsOpalBeamline_m.prepareSections();

    /*
    itsOpalBeamline_m.compute3DLattice();
    itsOpalBeamline_m.save3DLattice();
    itsOpalBeamline_m.save3DInput();
    */
    }

void ParallelTTracker::timeIntegration1(BorisPusher & pusher) {

    IpplTimings::startTimer(timeIntegrationTimer1_m);
    pushParticles(pusher);
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}


void ParallelTTracker::timeIntegration2(BorisPusher & pusher) {
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

    IpplTimings::startTimer(timeIntegrationTimer2_m);
    kickParticles(pusher);
    //switchElements();
    pushParticles(pusher);

    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTTracker::selectDT(bool backTrack) {

    double dt = dtCurrentTrack_m;
    itsBunch_m->setdT(dt);
    
    if (backTrack) {
        itsBunch_m->setdT(-std::abs(itsBunch_m->getdT()));
    }
}

void ParallelTTracker::changeDT(bool backTrack) {
    selectDT(backTrack);
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++ i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}

void ParallelTTracker::emitParticles(long long step) {

}


void ParallelTTracker::computeSpaceChargeFields(unsigned long long step) {
    if (!itsBunch_m->hasFieldSolver()) {
        return;
    }

    itsBunch_m->calcBeamParameters();
    Quaternion alignment = getQuaternion(itsBunch_m->get_pmean(), Vector_t(0, 0, 1));
    CoordinateSystemTrafo beamToReferenceCSTrafo(Vector_t(0, 0, pathLength_m), alignment.conjugate());
    CoordinateSystemTrafo referenceToBeamCSTrafo = beamToReferenceCSTrafo.inverted();
    const unsigned int localNum1 = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum1; ++ i) {
        itsBunch_m->R[i] = referenceToBeamCSTrafo.transformTo(itsBunch_m->R[i]);
    }

    itsBunch_m->boundp();

    if (step % repartFreq_m + 1 == repartFreq_m) {
        doBinaryRepartition();
    }

    itsBunch_m->setGlobalMeanR(itsBunch_m->get_centroid());
    itsBunch_m->computeSelfFields();
    
    const unsigned int localNum2 = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum2; ++ i) {
        itsBunch_m->R[i] = beamToReferenceCSTrafo.transformTo(itsBunch_m->R[i]);
        itsBunch_m->Ef[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Ef[i]);
        itsBunch_m->Bf[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Bf[i]);
    }
}


void ParallelTTracker::computeExternalFields(OrbitThreader &oth) {
    IpplTimings::startTimer(fieldEvaluationTimer_m);
    Inform msg("ParallelTTracker ", *gmsg);
    const unsigned int localNum = itsBunch_m->getLocalNum();
    bool locPartOutOfBounds = false, globPartOutOfBounds = false;
    Vector_t rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNum() > 0)
        itsBunch_m->get_bounds(rmin, rmax);
    IndexMap::value_t elements;

    try {
        elements = oth.query(pathLength_m + 0.5 * (rmax(2) + rmin(2)), rmax(2) - rmin(2));
    } catch(IndexMap::OutOfBounds &e) {
        globalEOL_m = true;
        return;
    }

    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        CoordinateSystemTrafo refToLocalCSTrafo = (itsOpalBeamline_m.getMisalignment((*it)) *
                                                   (itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * itsBunch_m->toLabTrafo_m));
        CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

        (*it)->setCurrentSCoordinate(pathLength_m + rmin(2));

        for (unsigned int i = 0; i < localNum; ++ i) {
            if (itsBunch_m->Bin[i] < 0) continue;

            itsBunch_m->R[i] = refToLocalCSTrafo.transformTo(itsBunch_m->R[i]);
            itsBunch_m->P[i] = refToLocalCSTrafo.rotateTo(itsBunch_m->P[i]);
            double &dt = itsBunch_m->dt[i];

            Vector_t localE(0.0), localB(0.0);

            if ((*it)->apply(i, itsBunch_m->getT() + 0.5 * dt, localE, localB)) {
                itsBunch_m->R[i] = localToRefCSTrafo.transformTo(itsBunch_m->R[i]);
                itsBunch_m->P[i] = localToRefCSTrafo.rotateTo(itsBunch_m->P[i]);
                itsBunch_m->Bin[i] = -1;
                locPartOutOfBounds = true;

                continue;
            }

            itsBunch_m->R[i] = localToRefCSTrafo.transformTo(itsBunch_m->R[i]);
            itsBunch_m->P[i] = localToRefCSTrafo.rotateTo(itsBunch_m->P[i]);
            itsBunch_m->Ef[i] += localToRefCSTrafo.rotateTo(localE);
            itsBunch_m->Bf[i] += localToRefCSTrafo.rotateTo(localB);
        }
    }

    IpplTimings::stopTimer(fieldEvaluationTimer_m);

    reduce(locPartOutOfBounds, globPartOutOfBounds, OpOrAssign());

    size_t ne = 0;
    if (globPartOutOfBounds) {
        if (itsBunch_m->hasFieldSolver()) {
            ne = itsBunch_m->boundp_destroyT();
        } else {
            ne = itsBunch_m->destroyT();
        }
        numParticlesInSimulation_m  = itsBunch_m->getTotalNum();
        deletedParticles_m = true;
    }

    size_t totalNum = itsBunch_m->getTotalNum();
    numParticlesInSimulation_m = totalNum;

    if (ne > 0) {
        msg << level1 << "* Deleted " << ne << " particles, "
            << "remaining " << numParticlesInSimulation_m << " particles" << endl;
    }
}

void ParallelTTracker::doBinaryRepartition() {
    if (itsBunch_m->hasFieldSolver()) {

        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition because of repartFreq_m" << endl);
        INFOMSG("*****************************************************************" << endl);
        IpplTimings::startTimer(BinRepartTimer_m);
        itsBunch_m->do_binaryRepart();
        Ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
        INFOMSG("*****************************************************************" << endl);
        INFOMSG("do repartition done" << endl);
        INFOMSG("*****************************************************************" << endl);
    }
}

void ParallelTTracker::dumpStats(long long step, bool psDump, bool statDump) {
    OPALTimer::Timer myt2;
    Inform msg("ParallelTTracker ", *gmsg);

    if (itsBunch_m->getGlobalTrackStep() % 1000 + 1 == 1000) {
        msg << level1;
    } else if (itsBunch_m->getGlobalTrackStep() % 100 + 1 == 100) {
        msg << level2;
    } else {
        msg << level3;
    }

    if (numParticlesInSimulation_m == 0) {
        msg << myt2.time() << " "
            << "Step " << std::setw(6) <<  itsBunch_m->getGlobalTrackStep() << "; "
            << "   -- no emission yet --     "
            << "t= "   << Util::getTimeString(itsBunch_m->getT())
            << endl;
        return;
    }

    itsBunch_m->calcEMean();
    size_t totalParticles_f = numParticlesInSimulation_m;
    if (std::isnan(pathLength_m) || std::isinf(pathLength_m)) {
        throw OpalException("ParallelTTracker::dumpStats()",
                            "there seems to be something wrong with the position of the bunch!");
    } else {

        msg << myt2.time() << " "
            << "Step " << std::setw(6)  <<  itsBunch_m->getGlobalTrackStep() << " "
            << "at "   << Util::getLengthString(pathLength_m) << ", "
            << "t= "   << Util::getTimeString(itsBunch_m->getT()) << ", "
            << "E="    << Util::getEnergyString(itsBunch_m->get_meanKineticEnergy())
            << endl;

        writePhaseSpace(step, psDump, statDump);
    }
}


void ParallelTTracker::setOptionalVariables() {
    Inform msg("ParallelTTracker ", *gmsg);

    minStepforReBin_m  = Options::minStepForRebin;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if (br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << level2 << "MINSTEPFORREBIN " << minStepforReBin_m << endl;

    // there is no point to do repartitioning with one node
    if (Ippl::getNodes() == 1) {
        repartFreq_m = std::numeric_limits<unsigned int>::max();
    } else {
        repartFreq_m = Options::repartFreq * 100;
        RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
        if (rep)
            repartFreq_m = static_cast<int>(rep->getReal());
        msg << level2 << "REPARTFREQ " << repartFreq_m << endl;
    }
}


bool ParallelTTracker::hasEndOfLineReached(const BoundingBox& globalBoundingBox) {
    reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
    globalEOL_m = globalEOL_m || globalBoundingBox.isOutside(itsBunch_m->RefPartR_m);
    return globalEOL_m;
}

void ParallelTTracker::setTime() {
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->dt[i] = itsBunch_m->getdT();
    }
}

void ParallelTTracker::writePhaseSpace(const long long /*step*/, bool psDump, bool statDump) {
    extern Inform *gmsg;
    Inform msg("OPAL ", *gmsg);
    Vector_t externalE, externalB;
    Vector_t FDext[2];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);

    if (psDump || statDump) {
        externalB = Vector_t(0.0);
        externalE = Vector_t(0.0);
        itsOpalBeamline_m.getFieldAt(itsBunch_m->RefPartR_m,
                                     itsBunch_m->RefPartP_m,
                                     itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(),
                                     externalE,
                                     externalB);
        FDext[0] = itsBunch_m->toLabTrafo_m.rotateFrom(externalB);
        FDext[1] = itsBunch_m->toLabTrafo_m.rotateFrom(externalE * Units::Vpm2MVpm);
    }

    if (statDump) {

        msg << level3 << "* Wrote beam statistics." << endl;
    }

    if (psDump && (itsBunch_m->getTotalNum() > 0)) {
        // Write fields to .h5 file.
        const size_t localNum = itsBunch_m->getLocalNum();
        double distToLastStop = stepSizes_m.getFinalZStop() - pathLength_m;
        Vector_t beta = itsBunch_m->RefPartP_m / Util::getGamma(itsBunch_m->RefPartP_m);
        Vector_t driftPerTimeStep = itsBunch_m->getdT() * Physics::c * itsBunch_m->toLabTrafo_m.rotateFrom(beta);
        bool driftToCorrectPosition = std::abs(distToLastStop) < 0.5 * euclidean_norm(driftPerTimeStep);
        Ppos_t stashedR;
        Vector_t stashedRefPartR;

        if (driftToCorrectPosition) {
            const double tau = distToLastStop / euclidean_norm(driftPerTimeStep) * itsBunch_m->getdT();
            if (localNum > 0) {
                stashedR.create(localNum);
                stashedR = itsBunch_m->R;
                stashedRefPartR = itsBunch_m->RefPartR_m;

                for (size_t i = 0; i < localNum; ++ i) {
                    itsBunch_m->R[i] += tau * (Physics::c * itsBunch_m->P[i] / Util::getGamma(itsBunch_m->P[i]) -
                                               driftPerTimeStep / itsBunch_m->getdT());
                }
            }

            driftPerTimeStep = itsBunch_m->toLabTrafo_m.rotateTo(driftPerTimeStep);
            itsBunch_m->RefPartR_m = itsBunch_m->RefPartR_m + tau * driftPerTimeStep / itsBunch_m->getdT();
            CoordinateSystemTrafo update(tau * driftPerTimeStep / itsBunch_m->getdT(),
                                         Quaternion(1.0, 0.0, 0.0, 0.0));
            itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();

            itsBunch_m->set_sPos(stepSizes_m.getFinalZStop());

            itsBunch_m->calcBeamParameters();
        }
        if (!statDump && !driftToCorrectPosition) itsBunch_m->calcBeamParameters();

        msg << *itsBunch_m << endl;
        itsDataSink_m->dumpH5(itsBunch_m, FDext);

        if (driftToCorrectPosition) {
            if (localNum > 0) {
                itsBunch_m->R = stashedR;
                stashedR.destroy(localNum, 0);
            }

            itsBunch_m->RefPartR_m = stashedRefPartR;
            itsBunch_m->set_sPos(pathLength_m);

            itsBunch_m->calcBeamParameters();
        }

        msg << level2 << "* Wrote beam phase space." << endl;
    }
}

void ParallelTTracker::updateReference(const BorisPusher &pusher) {
    updateReferenceParticle(pusher);
    updateRefToLabCSTrafo();
}

void ParallelTTracker::updateReferenceParticle(const BorisPusher &pusher) {
    const double direction = back_track? -1: 1;
    const double dt = direction * std::min(itsBunch_m->getT(),
                                           direction * itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;
    Vector_t Ef(0.0), Bf(0.0);

    itsBunch_m->RefPartR_m /= scaleFactor;
    pusher.push(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, dt);
    itsBunch_m->RefPartR_m *= scaleFactor;

    IndexMap::value_t elements = itsOpalBeamline_m.getElements(itsBunch_m->RefPartR_m);
    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        const CoordinateSystemTrafo &refToLocalCSTrafo = itsOpalBeamline_m.getCSTrafoLab2Local((*it));

        Vector_t localR = refToLocalCSTrafo.transformTo(itsBunch_m->RefPartR_m);
        Vector_t localP = refToLocalCSTrafo.rotateTo(itsBunch_m->RefPartP_m);
        Vector_t localE(0.0), localB(0.0);

        if ((*it)->applyToReferenceParticle(localR,
                                            localP,
                                            itsBunch_m->getT() - 0.5 * dt,
                                            localE,
                                            localB)) {
            *gmsg << level1 << "The reference particle hit an element" << endl;
            globalEOL_m = true;
        }

        Ef += refToLocalCSTrafo.rotateFrom(localE);
        Bf += refToLocalCSTrafo.rotateFrom(localB);
    }

    pusher.kick(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, Ef, Bf, dt);

    itsBunch_m->RefPartR_m /= scaleFactor;
    pusher.push(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, dt);
    itsBunch_m->RefPartR_m *= scaleFactor;
}

void ParallelTTracker::transformBunch(const CoordinateSystemTrafo &trafo) {
    const unsigned int localNum = itsBunch_m->getLocalNum();
    for (unsigned int i = 0; i < localNum; ++i) {
        itsBunch_m->R[i] = trafo.transformTo(itsBunch_m->R[i]);
        itsBunch_m->P[i] = trafo.rotateTo(itsBunch_m->P[i]);
        itsBunch_m->Ef[i] = trafo.rotateTo(itsBunch_m->Ef[i]);
        itsBunch_m->Bf[i] = trafo.rotateTo(itsBunch_m->Bf[i]);
    }
}

void ParallelTTracker::updateRefToLabCSTrafo() {
    Vector_t R = itsBunch_m->toLabTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
    Vector_t P = itsBunch_m->toLabTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);

    pathLength_m += std::copysign(1, itsBunch_m->getdT()) * euclidean_norm(R);

    CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t(0, 0, 1)));

    transformBunch(update);

    itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();
}

void ParallelTTracker::applyFractionalStep(const BorisPusher &pusher, double tau) {
    double t = itsBunch_m->getT();
    t += tau;
    itsBunch_m->setT(t);

    // the push method below pushes for half a time step. Hence the ref particle
    // should be pushed for 2 * tau
    itsBunch_m->RefPartR_m /= (Physics::c * 2 * tau);
    pusher.push(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, tau);
    itsBunch_m->RefPartR_m *= (Physics::c * 2 * tau);

    pathLength_m = zstart_m;

    Vector_t R = itsBunch_m->toLabTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
    Vector_t P = itsBunch_m->toLabTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);
    CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t(0, 0, 1)));
    itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();
}

void ParallelTTracker::findStartPosition(const BorisPusher &pusher) {

    StepSizeConfig stepSizesCopy(stepSizes_m);
    if (back_track) {
        stepSizesCopy.shiftZStopLeft(zstart_m);
    }

    double t = 0.0;
    itsBunch_m->setT(t);

    dtCurrentTrack_m = stepSizesCopy.getdT();
    selectDT();

    while (true) {
        autophaseCavities(pusher);

        t += itsBunch_m->getdT();
        itsBunch_m->setT(t);

        Vector_t oldR = itsBunch_m->RefPartR_m;
        updateReferenceParticle(pusher);
        pathLength_m += euclidean_norm(itsBunch_m->RefPartR_m - oldR);

        double speed = euclidean_norm(itsBunch_m->RefPartP_m * Physics::c / Util::getGamma(itsBunch_m->RefPartP_m));

        if (pathLength_m > stepSizesCopy.getZStop()) {
            ++ stepSizesCopy;

            if (stepSizesCopy.reachedEnd()) {
                -- stepSizesCopy;
                double tau = (stepSizesCopy.getZStop() - pathLength_m) / speed;
                applyFractionalStep(pusher, tau);

                break;
            }

            dtCurrentTrack_m = stepSizesCopy.getdT();
            selectDT();
        }

        if (std::abs(pathLength_m - zstart_m) <=  0.5 * itsBunch_m->getdT() * speed) {
            double tau = (zstart_m - pathLength_m) / speed;
            applyFractionalStep(pusher, tau);

            break;
        }
    }

    changeDT();
}

void ParallelTTracker::autophaseCavities(const BorisPusher &pusher) {

    double t = itsBunch_m->getT();
    Vector_t nextR = itsBunch_m->RefPartR_m / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, itsBunch_m->RefPartP_m, itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    for (auto element: elementSet) {
        if (element->getType() == ElementType::TRAVELINGWAVE) {
            const TravelingWave *TWelement = static_cast<const TravelingWave *>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, itsBunch_m->RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, itsBunch_m->RefPartP_m),
                                       t, itsBunch_m->getdT());
            }

        } else if (element->getType() == ElementType::RFCAVITY) {
            const RFCavity *RFelement = static_cast<const RFCavity *>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(itsOpalBeamline_m.transformToLocalCS(element, itsBunch_m->RefPartR_m),
                                       itsOpalBeamline_m.rotateToLocalCS(element, itsBunch_m->RefPartP_m),
                                       t, itsBunch_m->getdT());
            }
        }
    }
}

struct DistributionInfo {
    unsigned int who;
    unsigned int whom;
    unsigned int howMany;
};
