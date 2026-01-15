//
// Class ParallelTracker
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
#include "Algorithms/ParallelTracker.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <boost/numeric/ublas/io.hpp>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/CavityAutophaser.h"
#include "Algorithms/OrbitThreader.h"
#include "BasicActions/Option.h"

#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedBeamline.h"
#include "Distribution/Distribution.h"
#include "Physics/Units.h"

#include "Structure/BoundaryGeometry.h"
#include "Structure/BoundingBox.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "ValueDefinitions/RealVariable.h"

#include "AbsBeamline/Offset.h"
#include "AbsBeamline/PluginElement.h"
#include "AbsBeamline/VerticalFFAMagnet.h"

extern Inform* gmsg;

class PartData;
/* ============================== Constructors ============================== */
ParallelTracker::ParallelTracker(
    const Beamline& beamline, const PartData& reference, bool revBeam, bool revTrack)
    : Tracker(beamline, reference, revBeam, revTrack),
      itsDataSink_m(nullptr),
      itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
      opalRing_m(nullptr),
      globalEOL_m(false),
      wakeStatus_m(false),
      wakeFunction_m(nullptr),
      pathLength_m(0.0),
      zstart_m(0.0),
      dtCurrentTrack_m(0.0),
      minStepforReBin_m(-1),
      repartFreq_m(0),
      emissionSteps_m(std::numeric_limits<unsigned int>::max()),
      numParticlesInSimulation_m(0),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
      PluginElemTimer_m(IpplTimings::getTimer("PluginElements")),
      BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
      OrbThreader_m(IpplTimings::getTimer("OrbThreader")) {
}

ParallelTracker::ParallelTracker(
    const Beamline& beamline, PartBunch_t* bunch, DataSink& ds, const PartData& reference,
    bool revBeam, bool revTrack, const std::vector<unsigned long long>& maxSteps, double zstart,
    const std::vector<double>& zstop, const std::vector<double>& dt)
    : Tracker(beamline, bunch, reference, revBeam, revTrack),
      itsDataSink_m(&ds),
      itsOpalBeamline_m(beamline.getOrigin3D(), beamline.getInitialDirection()),
      opalRing_m(nullptr),
      globalEOL_m(false),
      wakeStatus_m(false),
      wakeFunction_m(nullptr),
      pathLength_m(0.0),
      zstart_m(zstart),
      dtCurrentTrack_m(0.0),
      minStepforReBin_m(-1),
      repartFreq_m(0),
      emissionSteps_m(std::numeric_limits<unsigned int>::max()),
      numParticlesInSimulation_m(0),
      timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
      timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
      fieldEvaluationTimer_m(IpplTimings::getTimer("External field eval")),
      BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
      OrbThreader_m(IpplTimings::getTimer("OrbThreader")) {
    
      for (unsigned int i = 0; i < zstop.size(); ++i) {
          stepSizes_m.push_back(dt[i], zstop[i], maxSteps[i]);
      }

      stepSizes_m.sortAscendingZStop();
      stepSizes_m.resetIterator();
}

ParallelTracker::~ParallelTracker() {
}
/* ========================================================================== */
/* =========================== Visit Functions ============================== */
/**
 * @brief Iterates over the list of elements in TBeamline& bl and calls
 * the overloaded accept() function for each element.
 *
 * @param bl A reference to a TBeamline object, which holds the list of elements
 */
void ParallelTracker::visitBeamline(const Beamline& bl) {
    const FlaggedBeamline* fbl = static_cast<const FlaggedBeamline*>(&bl);
    /*
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
    */
    fbl->iterate(*this, false);
}

void ParallelTracker::visitScalingFFAMagnet(const ScalingFFAMagnet& bend) {
    *gmsg << "Adding ScalingFFAMagnet" << endl;
    /*
    if (opalRing_m != nullptr) {
        ScalingFFAMagnet* newBend = bend.clone(); // setup the end field, if required
        newBend->setupEndField();
        opalRing_m->appendElement(*newBend);
    } else {
        throw OpalException("ParallelCyclotronTracker::visitScalingFFAMagnet",
                            "Need to define a RINGDEFINITION to use ScalingFFAMagnet element");
    }
    */
}

void ParallelTracker::visitRing(const Ring& ring) {
    *gmsg << "* ----------------------------- Ring ------------------------------------- *" << endl;

    delete opalRing_m;

    opalRing_m = dynamic_cast<Ring*>(ring.clone());

    myElements.push_back(opalRing_m);

    opalRing_m->initialise(itsBunch_m);

    referenceR     = opalRing_m->getBeamRInit();
    referencePr    = opalRing_m->getBeamPRInit();
    referenceTheta = opalRing_m->getBeamPhiInit();

    if (referenceTheta <= -180.0 || referenceTheta > 180.0) {
        throw OpalException(
            "Error in ParallelTracker::visitRing", "PHIINIT is out of [-180, 180)!");
    }

    referenceZ  = 0.0;
    referencePz = 0.0;

    referencePtot = itsReference.getGamma() * itsReference.getBeta();
    referencePt   = std::sqrt(referencePtot * referencePtot - referencePr * referencePr);

    if (referencePtot < 0.0)
        referencePt *= -1.0;

    sinRefTheta_m = std::sin(referenceTheta * Units::deg2rad);
    cosRefTheta_m = std::cos(referenceTheta * Units::deg2rad);

    double BcParameter[8] = {};  // zero initialise array

    buildupFieldList(BcParameter, ElementType::RING, opalRing_m);

    // Finally print some diagnostic
    *gmsg << "* Initial beam radius = " << referenceR << " [mm] " << endl;
    *gmsg << "* Initial gamma = " << itsReference.getGamma() << endl;
    *gmsg << "* Initial beta  = " << itsReference.getBeta() << endl;
    *gmsg << "* Total reference momentum      = " << referencePtot << " [beta gamma]" << endl;
    *gmsg << "* Reference azimuthal momentum  = " << referencePt << " [beta gamma]" << endl;
    *gmsg << "* Reference radial momentum     = " << referencePr << " [beta gamma]" << endl;
    *gmsg << "* " << opalRing_m->getSymmetry() << " fold field symmetry " << endl;
    *gmsg << "* Harmonic number h = " << opalRing_m->getHarmonicNumber() << " " << endl;
}

void ParallelTracker::visitVerticalFFAMagnet(const VerticalFFAMagnet& mag) {
    *gmsg << "Adding Vertical FFA Magnet" << endl;
    if (opalRing_m != nullptr)
        opalRing_m->appendElement(mag);
    else
        throw OpalException(
            "ParallelCyclotronTracker::visitVerticalFFAMagnet",
            "Need to define a RINGDEFINITION to use VerticalFFAMagnet element");
}

/** * * * @param off */
void ParallelTracker::visitOffset(const Offset& off) {
    if (opalRing_m == nullptr)
        throw OpalException(
            "ParallelCylcotronTracker::visitOffset",
            "Attempt to place an offset when Ring not defined");
    Offset* offNonConst = const_cast<Offset*>(&off);
    offNonConst->updateGeometry(opalRing_m->getNextPosition(), opalRing_m->getNextNormal());
    opalRing_m->appendElement(off);
}


/* ========================================================================== */
/* =========================== Start Simulation ============================= */
void ParallelTracker::execute() {
    // Inform object
    Inform msg("ParallelTracker ", *gmsg);

    // Set preparation state in OpalData
    OpalData::getInstance()->setInPrepState(true);

    // Flag to indicate forward tracking
    bool back_track = false;

    // Initialize the Boris particle pusher
    BorisPusher pusher(itsReference);

    // Reset the global phase shift
    OpalData::getInstance()->setGlobalPhaseShift(0.0);

    // Ensure the time step is positive for the setup phase
    itsBunch_m->setdT(std::abs(itsBunch_m->getdT()));

    // Set the time step for the current track
    dtCurrentTrack_m = itsBunch_m->getdT();
    
    // If there is prior tracking data or if in a restart run, set the data
    // sink to append mode
    if (OpalData::getInstance()->hasPriorTrack() || 
        OpalData::getInstance()->inRestartRun()) 
    {
        OpalData::getInstance()->setOpenMode(OpalData::OpenMode::APPEND);
    }

    // Populate the OpalBeamline and calculate coordinate transformations 
    // for non-straight sections
    prepareSections();

    // Select the minimal time step from the configuration
    double minTimeStep = stepSizes_m.getMinTimeStep();

    // Activate all beamline elements (sets Component::online_m = true)
    itsOpalBeamline_m.activateElements();

    // Calculate the coordinate transformation from the beamline origin to the lab frame
    CoordinateSystemTrafo beamlineToLab = 
        itsOpalBeamline_m.getCSTrafoLab2Local().inverted();

    itsBunch_m->toLabTrafo_m = beamlineToLab;

    // Transform reference particle coordinates and momentum to the lab frame
    itsBunch_m->RefPartR_m = 
        beamlineToLab.transformTo(itsBunch_m->getParticleContainer()->getMeanR());
    itsBunch_m->RefPartP_m = 
        beamlineToLab.rotateTo(itsBunch_m->getParticleContainer()->getMeanP());

    // If the bunch contains particles and the desired starting position (zstart_m) 
    // is ahead of the current path length, integrate the bunch forward to the start position
    if (itsBunch_m->getTotalNum() > 0) {
        if (zstart_m > pathLength_m) {
            findStartPosition(pusher);
        }
        itsBunch_m->set_sPos(pathLength_m);
    }

    // Advance through the time step configurations, skipping any that end 
    // before the start position (zstart_m)
    stepSizes_m.advanceToPos(zstart_m);

    // Retrieve the bunch spatial bounds 
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNum() > 0) {
        itsBunch_m->get_bounds(rmin, rmax);
    }

    *gmsg << "ParallelTrack: momentum=  " 
        << itsBunch_m->getParticleContainer()->getMeanP()(2) << endl;
    *gmsg << "itsBunch_m->RefPartR_m= " 
        << itsBunch_m->RefPartR_m << endl;
    *gmsg << "itsBunch_m->RefPartP_m= " 
        << itsBunch_m->RefPartP_m << endl;

    
    // Start timing for the OrbitThreader section
    IpplTimings::startTimer(OrbThreader_m);

    // Flags to control phase space and statistics dumping for the initial step
    bool const psDump0 = 0;
    bool const statDump0 = 0;

    // Write initial phase space and statistics
    writePhaseSpace(0, psDump0, statDump0);
    msg << level2 << "Dump initial phase space" << endl;

    // Create an OrbitThreader object to handle orbit threading and element queries
    OrbitThreader oth(
        itsReference,                               // Reference particle data
        itsBunch_m->RefPartR_m,                     // Reference particle position
        itsBunch_m->RefPartP_m,                     // Reference particle momentum
        pathLength_m,                               // Current path length
        -rmin(2),                                   // Negative minimum z bound
        itsBunch_m->getT(),                         // Current bunch time
        (back_track ? -minTimeStep : minTimeStep),  // Time step direction
        stepSizes_m,                                // Step size configuration
        itsOpalBeamline_m);                         // OpalBeamline object

    // Execute the orbit threading (queries elements and updates state)
    oth.execute();

    // Stop timing for the OrbitThreader section
    IpplTimings::stopTimer(OrbThreader_m);
    
    // Get bounding box
    BoundingBox globalBoundingBox = oth.getBoundingBox();

    // Total particle number
    numParticlesInSimulation_m = itsBunch_m->getTotalNum();

    // Set the time view of the particle bunch
    setTime();

    // Reset the bunch time?
    double time = itsBunch_m->getT();
    itsBunch_m->setT(time);

    // Get the current global tracking step
    unsigned long long step = itsBunch_m->getGlobalTrackStep();
    OPALTimer::Timer myt1;
    *gmsg << "* Track start at: " << myt1.time() << ", t= " << Util::getTimeString(time) << "; "
          << "zstart at: " << Util::getLengthString(pathLength_m) << endl;
    *gmsg << "* Executing ParallelTracker\n"
          << "* Initial dt = " << Util::getTimeString(itsBunch_m->getdT()) << "\n"
          << "* Max integration steps = " << stepSizes_m.getMaxSteps() << ", next step = " << step
          << endl
          << endl;

    setOptionalVariables();

    globalEOL_m        = false;
    wakeStatus_m       = false;
    deletedParticles_m = false;
    OpalData::getInstance()->setInPrepState(false);

    stepSizes_m.printDirect(*gmsg);

    /// Directly before the tracker loop, perform bunch sanity checks
    this->itsBunch_m->performBunchSanityChecks();
    
    // Main tracking loop over step size configurations
    while (!stepSizes_m.reachedEnd()) {

        // Set the number of steps for the current track
        unsigned long long trackSteps = stepSizes_m.getNumSteps() + step;
        dtCurrentTrack_m              = stepSizes_m.getdT();
        changeDT(back_track);

        // Inner loop over the number of steps for the current configuration
        for (; step < trackSteps; ++step) {

            // Get the bunch spatial bounds
            Vector_t<double, 3> rmin(0.0), rmax(0.0);
            if (itsBunch_m->getTotalNum() > 0) {
                itsBunch_m->calcBeamParameters();
                itsBunch_m->get_bounds(rmin, rmax);
            }
 
            // First half of the time integration
            timeIntegration1(pusher);

            // Reset E and B fields
            resetFields();
            
            // Space charge field computation
            //computeSpaceChargeFields(step);
           
            // External field computation
            computeExternalFields(oth);

            // Second half of the time integration
            timeIntegration2(pusher);
            
            // \todo emitParticles(step);

            // Backtracking?
            selectDT(back_track);
            
            // Update the bunch time
            itsBunch_m->incrementT();

            // Update reference particle
            if (itsBunch_m->getT() > 0.0 || itsBunch_m->getdT() < 0.0) {
                updateReference(pusher);
            }

            // Delete particles?
            if (deletedParticles_m) {
                // \todo doDelete
                deletedParticles_m = false;
            }

            // Update the path length
            itsBunch_m->set_sPos(pathLength_m);

            // if (hasEndOfLineReached(globalBoundingBox)) break;
  
            // Dump phase space and statistics at configured intervals
            bool const psDump =
                ((itsBunch_m->getGlobalTrackStep() % Options::psDumpFreq) + 1
                 == Options::psDumpFreq);
            bool const statDump =
                ((itsBunch_m->getGlobalTrackStep() % Options::statDumpFreq) + 1
                 == Options::statDumpFreq);
            dumpStats(step, psDump, statDump);

            itsBunch_m->incTrackSteps();

            ippl::Vector<double,3> pdivg = itsBunch_m->RefPartP_m / Util::getGamma(itsBunch_m->RefPartP_m);
            double beta = euclidean_norm(pdivg);
            double driftPerTimeStep = std::abs(itsBunch_m->getdT()) * Physics::c * beta;

            if (std::abs(stepSizes_m.getZStop() - pathLength_m) < 0.5 * driftPerTimeStep) {
                break;
            }
        }

        if (globalEOL_m)
            break;
        ++stepSizes_m;    
     }
    itsBunch_m->set_sPos(pathLength_m);

    numParticlesInSimulation_m = itsBunch_m->getTotalNum();

    bool const psDump =
        (((itsBunch_m->getGlobalTrackStep() - 1) % Options::psDumpFreq) + 1 != Options::psDumpFreq);
    bool const statDump =
        (((itsBunch_m->getGlobalTrackStep() - 1) % Options::statDumpFreq) + 1
         != Options::statDumpFreq);

    writePhaseSpace((step + 1), psDump, statDump);

    *gmsg << "* Dump phase space of last step" << endl;

    itsOpalBeamline_m.switchElementsOff();

    // Ensure all Kokkos operations are complete 
    Kokkos::fence();

    OPALTimer::Timer myt3;
    *gmsg << endl << "* Done executing ParallelTracker at " << myt3.time() << endl << endl;
}
/* ========================================================================== */
/* =========================== PIC Functions ================================ */
void ParallelTracker::timeIntegration1(BorisPusher& pusher) {
    IpplTimings::startTimer(timeIntegrationTimer1_m);
    pushParticles(pusher);
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}

void ParallelTracker::timeIntegration2(BorisPusher& pusher) {
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
    // switchElements();
    pushParticles(pusher);

    
    auto dtview  = itsBunch_m->getParticleContainer()->dt.getView();
    double newdT = itsBunch_m->getdT();

    Kokkos::parallel_for("changeDT", ippl::getRangePolicy(dtview),
    KOKKOS_LAMBDA(const int i) {
        dtview(i) = newdT;
    });                     
    
    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}
void ParallelTracker::computeSpaceChargeFields(unsigned long long step) {

    if (!itsBunch_m->hasFieldSolver()) {
        *gmsg << "no solver avaidable " << endl;
        return;
    }
        
    itsBunch_m->calcBeamParameters();

    Quaternion alignment = getQuaternion(itsBunch_m->get_pmean(), Vector_t<double, 3>(0, 0, 1));

    CoordinateSystemTrafo beamToReferenceCSTrafo(
        Vector_t<double, 3>(0, 0, pathLength_m), alignment.conjugate());

    CoordinateSystemTrafo referenceToBeamCSTrafo = beamToReferenceCSTrafo.inverted();

    /**
       \brief referernce to beam coordinate system transformation 
       old: referenceToBeamCSTrafo.transformTo

       x = Rot * (x-o)

     */

    
    const OpalMatrix_t                rot = referenceToBeamCSTrafo.getRotationMatrix();
    const ippl::Vector<double, 3> org = referenceToBeamCSTrafo.getOrigin();


    typedef Kokkos::View<double**>  ViewMatrixType;
    ViewMatrixType Rot("Rot", 3, 3);
    ViewMatrixType::HostMirror h_Rot = Kokkos::create_mirror_view(Rot);

    // Initialize M matrix on host.
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 0; j < 3; ++j ) {
            h_Rot( i, j ) = rot(i, j);
        }
    }
    
    Kokkos::deep_copy( Rot, h_Rot );

    auto Rview  = itsBunch_m->getParticleContainer()->R.getView();
    auto Eview  = itsBunch_m->getParticleContainer()->E.getView();
    auto Bview  = itsBunch_m->getParticleContainer()->B.getView();

    Kokkos::parallel_for(
        "referenceToBeamCSTrafo", ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const size_t k) {
            ippl::Vector<double, 3> x = Rview(k); // x({Rview(i)[0],Rview(i)[1],Rview(i)[2]});
                x = x - org;
                for (size_t i = 0; i < 3; ++i ) {
                    Rview(k)[i] = 0.0;
                    for (size_t j = 0; j < 3; ++j ) {
                        Rview(k)[i] += Rot(i, j) * x(j);
                    }
                }
                // Eview(k) = 0; // ippl::Vector<double, 3>(0.0);   // was done outside of the routine in the past
                // Bview(k) = 0; // ippl::Vector<double, 3>(0.0); 
        });         

    itsBunch_m->boundp();

    if (step % repartFreq_m + 1 == repartFreq_m) {
        doBinaryRepartition();
    }

    itsBunch_m->setGlobalMeanR(itsBunch_m->get_centroid());

    itsBunch_m->computeSelfFields();

    /**
        \brief beam to referernce coordinate system transformation and field rotation  
        old: itsBunch_m->R[i]  = beamToReferenceCSTrafo.transformTo(itsBunch_m->R[i]);
            itsBunch_m->Ef[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Ef[i]);
            itsBunch_m->Bf[i] = beamToReferenceCSTrafo.rotateTo(itsBunch_m->Bf[i]);

            x = M^T(x+o)
            prod_boost_vector(boost::numeric::ublas::trans(rotationMatrix_m)

            Update 17.12.2025: It looks like, there were a few errors with this transformation after trying to get
            rid of boost. This should now be fixed.
     */
    
    Kokkos::parallel_for(
        "beamToReferenceCSTrafo", ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const size_t k) {                           
            ippl::Vector<double, 3> x = Rview(k); // ({Rview(k)[0],Rview(k)[1],Rview(k)[2]});
            ippl::Vector<double, 3> e = Eview(k); // ({Eview(k)[0],Eview(k)[1],Eview(k)[2]});
            ippl::Vector<double, 3> b = Bview(k); // ({Bview(k)[0],Bview(k)[1],Bview(k)[2]});
            
            // beamToReferenceCSTrafo.rotateTo
            for (size_t i = 0; i < 3; ++i ) {
                Eview(k)[i] = 0.0;
                Bview(k)[i] = 0.0;
                for (size_t j = 0; j < 3; ++j ) {
                    Eview(k)[i] += Rot(j,i) * e(j);
                    Bview(k)[i] += Rot(j,i) * b(j);
                }
            }
            // beamToReferenceCSTrafo.transformTo
            x = x + org;
            for (size_t i = 0; i < 3; ++i ) {
                Rview(k)[i] = 0.0;
                for (size_t j = 0; j < 3; ++j ) {
                    Rview(k)[i] += Rot(j,i) * x(j);
                }
            }
        });         

}

void ParallelTracker::computeExternalFields(OrbitThreader& oth) {
    IpplTimings::startTimer(fieldEvaluationTimer_m);
    Inform msg("ParallelTracker ", *gmsg);

    // local # particles
    const unsigned int localNum = itsBunch_m->getLocalNum();

    // Flag for out-of-bounds particles, locally and globally
    bool locPartOutOfBounds = false, globPartOutOfBounds = false;
    
    // Bunch bounds
    Vector_t<double, 3> rmin(0.0), rmax(0.0);
    if (itsBunch_m->getTotalNum() > 0)
        itsBunch_m->get_bounds(rmin, rmax);

    // Get elements at bunch position
    IndexMap::value_t elements;
    try {
        elements = oth.query(pathLength_m + 0.5 * 
            (rmax(2) + rmin(2)), rmax(2) - rmin(2));
    } catch (IndexMap::OutOfBounds& e) {
        globalEOL_m = true;
        return;
    }

    // Iterator all elements at the current position
    IndexMap::value_t::const_iterator it        = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    // Iterate over all elements
    for (; it != end; ++it) {

        // Determine transformation from bunch to element 
        CoordinateSystemTrafo refToLocalCSTrafo = 
            (itsOpalBeamline_m.getMisalignment((*it)) * 
            (itsOpalBeamline_m.getCSTrafoLab2Local((*it)) * 
            itsBunch_m->toLabTrafo_m));

        // Determine transformation from element to bunch
        CoordinateSystemTrafo localToRefCSTrafo = refToLocalCSTrafo.inverted();

        (*it)->setCurrentSCoordinate(pathLength_m + rmin(2));   

        // Transform from reference particle to element frame
        refToLocalCSTrafo.transformBunchTo(
            itsBunch_m->getParticleContainer()->R.getView());
        
        refToLocalCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->P.getView());

        refToLocalCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->E.getView());

        refToLocalCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->B.getView());

        // Apply element
        (*it)->apply(); 

        // Transform from element to reference particle frame
        localToRefCSTrafo.transformBunchTo(
            itsBunch_m->getParticleContainer()->R.getView());
        localToRefCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->P.getView());
        localToRefCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->E.getView());
        localToRefCSTrafo.rotateBunchTo(
            itsBunch_m->getParticleContainer()->B.getView());

    }


    IpplTimings::stopTimer(fieldEvaluationTimer_m);

    ippl::Comm->reduce(locPartOutOfBounds, globPartOutOfBounds, 1, std::logical_or<bool>());

    size_t ne = 0;
    if (globPartOutOfBounds) {
        if (itsBunch_m->hasFieldSolver()) {
            ne = itsBunch_m->boundp_destroyT();
        }

        // \todo
        // else {
        //    ne = itsBunch_m->destroyT();
        // }
        numParticlesInSimulation_m = itsBunch_m->getTotalNum();
        deletedParticles_m = true;
    }

    size_t totalNum            = itsBunch_m->getTotalNum();
    numParticlesInSimulation_m = totalNum;

    if (ne > 0) {
        msg << level1 << "* Deleted " << ne << " particles, "
            << "remaining " << numParticlesInSimulation_m << " particles" << endl;
    }
}
/**
 * @brief Resets the E and B field views to 0
 */
void ParallelTracker::resetFields() {    
    itsBunch_m->getParticleContainer()->E = 0;
    itsBunch_m->getParticleContainer()->B = 0;
}

/**
 * @brief Pushes particles 
 * 
 * @param pusher The BorisPusher
 */
void ParallelTracker::pushParticles(const BorisPusher& pusher) {

    /// \todo use false for now, since I am not sure how well integrated "dt_per_particle" is (needs to be consistent with particle emission later!).
    itsBunch_m->switchToUnitlessPositions(true);

    auto Rview  = itsBunch_m->getParticleContainer()->R.getView();
    auto Pview  = itsBunch_m->getParticleContainer()->P.getView();
    //auto dtview = itsBunch_m->getParticleContainer()->dt.getView();

    Kokkos::parallel_for(
        "pushParticles", ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const size_t i) {
            /** 
             * \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} +
             * \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
             *
             * \code
             * R[i] += 0.5 * P[i] * recpgamma;
             * \endcode
             */
            /// \todo check +-
            
            Vector_t<double, 3> x = Rview(i); 
            pusher.push(x, Pview(i), 0); // this 0 is "dt" that is not used with unitless positions!
            Rview(i) = x;
        });

    itsBunch_m->switchOffUnitlessPositions(true);
    /// \todo update gives different results on one rank?
    //itsBunch_m->getParticleContainer()->update();
    ippl::Comm->barrier();
}

/**
 * @brief Kicks particles
 * 
 * @param pusher The BorisPusher
 */
void ParallelTracker::kickParticles(const BorisPusher& pusher) {

    // auto Rview  = itsBunch_m->getParticleContainer()->R.getView();
    auto Pview  = itsBunch_m->getParticleContainer()->P.getView();
    auto dtview = itsBunch_m->getParticleContainer()->dt.getView(); 
    auto Efview = itsBunch_m->getParticleContainer()->E.getView();
    auto Bfview = itsBunch_m->getParticleContainer()->B.getView();

    /// \todo Apparently, we want mass in eV and charge in elementary charges here to match OPAL's BorisPusher
    //double mass = itsBunch_m->getMassPerParticle(); // itsReference.getM();
    //double charge = itsBunch_m->getChargePerParticle();  // itsReference.getQ();
        // Get reference particle mass and charge
    const double mass = itsReference.getM();
    const double charge = itsReference.getQ();


    std::cout << charge << " " << mass << std::endl;

    Kokkos::parallel_for(
        "kickParticles", ippl::getRangePolicy(Pview),
        KOKKOS_LAMBDA(const size_t i) {
            /**
             *
             * Implementation follows chapter 4-4, pp. 61–63, from:
             * Birdsall, C. K. and Langdon, A. B. (1985). Plasma Physics via Computer Simulation.
             *
             * Up to finite-precision effects, the new implementation is equivalent to the old one,
             * but uses fewer floating-point operations.
             *
             * The relativistic variant implemented below is described in chapter 15-4, pp. 356–357.
             * However, since different units are used here, small modifications are required.
             * The relativistic variant can be derived from the nonrelativistic one by replacing:
             *   mass
             * with:
             *   gamma * rest mass
             * and transforming the units accordingly.
             *
             * Parameters:
             *   R      Scaled position, R = x / (c * dt). Not used here.
             *   P      Scaled velocity, P = (v / c) * gamma.
             *   Ef     Electric field.
             *   Bf     Magnetic field.
             *   dt     Timestep.
             *   mass   Rest energy, i.e., rest mass * c^2.
             *   charge Particle charge.
             *
             */
           
            Vector_t<double, 3> p = Pview(i); 
            pusher.kick(0, p, Efview(i), Bfview(i), dtview(i), mass, charge); /// \todo might want to remove dt and R altogether from the kick!
            Pview(i) = p; 
        });
        
    /// \todo unnecessary update? kick does not modify positions
    //itsBunch_m->getParticleContainer()->update();
    ippl::Comm->barrier();
}
/* ========================================================================== */ 
/* ============================= Functions ================================== */
/**
 * @brief Sets up beamline
 */
void ParallelTracker::prepareSections() {
    // Calls ParallelTracker::visitBeamline() -> TBeamline::iterate() ->
    // For each element: 
    // FlaggedElemPtr::accept() -> DefaultVisitor::visitFlaggedElmPtr() ->
    // ElementBase::accept() -> ParallelTracker::visit[ElemName]() ->
    // OpalBeamline::visit([ElemName], bunch):
    // This initialises the FieldList elements_m object of OpalBeamline
    // with clones of all the elements
    itsBeamline_m.accept(*this);
    
    // Sorts the elements in OpalBeamline::elements_m by starting position
    itsOpalBeamline_m.prepareSections();

    // Computes the coordinate transformations for non-straight sections
    itsOpalBeamline_m.compute3DLattice();

    // Write 3D Lattice
    itsOpalBeamline_m.save3DLattice();
    itsOpalBeamline_m.save3DInput();
}
void ParallelTracker::selectDT(bool backTrack) {
    double dt = dtCurrentTrack_m;
    itsBunch_m->setdT(dt);

    if (backTrack) {
        itsBunch_m->setdT(-std::abs(itsBunch_m->getdT()));
    }
}

void ParallelTracker::changeDT(bool backTrack) {

    selectDT(backTrack);

    auto dtview  = itsBunch_m->getParticleContainer()->dt.getView();
    double newdT = itsBunch_m->getdT();

    Kokkos::parallel_for(
                         "changeDT", ippl::getRangePolicy(dtview),
                         KOKKOS_LAMBDA(const int i) {
                             dtview(i) = newdT;
                         });                     
    
}

void ParallelTracker::doBinaryRepartition() {
    if (itsBunch_m->hasFieldSolver()) {
        *ippl::Info << "*****************************************************************" << endl;
        *ippl::Info << "do repartition because of repartFreq_m" << endl;
        *ippl::Info << "*****************************************************************" << endl;
        itsBunch_m->do_binaryRepart();
        ippl::Comm->barrier();
        IpplTimings::stopTimer(BinRepartTimer_m);
        *ippl::Info << "*****************************************************************" << endl;
        *ippl::Info << "do repartition done" << endl;
        *ippl::Info << "*****************************************************************" << endl;
    }
}

void ParallelTracker::updateReference(const BorisPusher& pusher) {
    updateReferenceParticle(pusher);
    updateRefToLabCSTrafo();
}

void ParallelTracker::updateReferenceParticle(const BorisPusher& pusher) {
    const double direction = back_track ? -1 : 1;
    const double dt = direction * std::min(itsBunch_m->getT(), direction * itsBunch_m->getdT());
    const double scaleFactor = Physics::c * dt;
    Vector_t<double, 3> Ef(0.0), Bf(0.0);

    itsBunch_m->RefPartR_m /= scaleFactor;
    pusher.push(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, dt);
    itsBunch_m->RefPartR_m *= scaleFactor;

    IndexMap::value_t elements           = itsOpalBeamline_m.getElements(itsBunch_m->RefPartR_m);
    IndexMap::value_t::const_iterator it = elements.begin();
    const IndexMap::value_t::const_iterator end = elements.end();

    for (; it != end; ++it) {
        const CoordinateSystemTrafo& refToLocalCSTrafo =
            itsOpalBeamline_m.getCSTrafoLab2Local((*it));

        Vector_t<double, 3> localR = refToLocalCSTrafo.transformTo(itsBunch_m->RefPartR_m);
        Vector_t<double, 3> localP = refToLocalCSTrafo.rotateTo(itsBunch_m->RefPartP_m);
        Vector_t<double, 3> localE(0.0), localB(0.0);

        if ((*it)->applyToReferenceParticle(
                localR, localP, itsBunch_m->getT() - 0.5 * dt, localE, localB)) {
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

void ParallelTracker::transformBunch(const CoordinateSystemTrafo& trafo) 
{
    trafo.transformBunchTo(itsBunch_m->getParticleContainer()->R.getView());
    trafo.rotateBunchTo(itsBunch_m->getParticleContainer()->P.getView());
    trafo.rotateBunchTo(itsBunch_m->getParticleContainer()->E.getView());
    trafo.rotateBunchTo(itsBunch_m->getParticleContainer()->B.getView());
}

void ParallelTracker::updateRefToLabCSTrafo() {
    /*
    Inform m("updateRefToLabCSTrafo ");
    m << " initial RefPartR: " << itsBunch_m->RefPartR_m << endl;
    m << " RefPartR after transformFrom( " << itsBunch_m->RefPartR_m << endl;
    itsBunch_m->toLabTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);
    Vector_t<double, 3> P = itsBunch_m->RefPartP_m;
    m << " CS " << itsBunch_m->toLabTrafo_m << endl;
    m << " rotMat= " << itsBunch_m->toLabTrafo_m.getRotationMatrix() << endl;
    m << " initial: path Lenght= " << pathLength_m << endl;
    m << " dt= " << itsBunch_m->getdT() << " R=" << R << endl;
    */

    Vector_t<double, 3> R = itsBunch_m->toLabTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
    Vector_t<double, 3> P = itsBunch_m->toLabTrafo_m.transformFrom(itsBunch_m->RefPartP_m);

    pathLength_m += std::copysign(1, itsBunch_m->getdT()) * euclidean_norm(R);

    CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t<double, 3>(0, 0, 1)));

    transformBunch(update);

    itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();
}

void ParallelTracker::applyFractionalStep(const BorisPusher& pusher, double tau) {
    double t = itsBunch_m->getT();
    t += tau;
    itsBunch_m->setT(t);

    // the push method below pushes for half a time step. Hence the ref particle
    // should be pushed for 2 * tau
    itsBunch_m->RefPartR_m /= (Physics::c * 2 * tau);
    pusher.push(itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m, tau);
    itsBunch_m->RefPartR_m *= (Physics::c * 2 * tau);

    pathLength_m = zstart_m;
    itsBunch_m->toLabTrafo_m.transformFrom(itsBunch_m->RefPartR_m);
    Vector_t<double, 3> R = itsBunch_m->RefPartR_m;
    itsBunch_m->toLabTrafo_m.rotateFrom(itsBunch_m->RefPartP_m);
    Vector_t<double, 3> P = itsBunch_m->RefPartP_m;
    CoordinateSystemTrafo update(R, getQuaternion(P, Vector_t<double, 3>(0, 0, 1)));
    itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();
}

void ParallelTracker::findStartPosition(const BorisPusher& pusher) {
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

        Vector_t<double, 3> oldR = itsBunch_m->RefPartR_m;
        updateReferenceParticle(pusher);

        // \todo should not need the tmp
        Vector_t<double, 3> tmp = itsBunch_m->RefPartR_m - oldR;
        pathLength_m += euclidean_norm(tmp);

        tmp = itsBunch_m->RefPartP_m * Physics::c / Util::getGamma(itsBunch_m->RefPartP_m);
        double speed = euclidean_norm(tmp);
                                      

        if (pathLength_m > stepSizesCopy.getZStop()) {
            ++stepSizesCopy;

            if (stepSizesCopy.reachedEnd()) {
                --stepSizesCopy;
                double tau = (stepSizesCopy.getZStop() - pathLength_m) / speed;
                applyFractionalStep(pusher, tau);

                break;
            }

            dtCurrentTrack_m = stepSizesCopy.getdT();
            selectDT();
        }

        if (std::abs(pathLength_m - zstart_m) <= 0.5 * itsBunch_m->getdT() * speed) {
            double tau = (zstart_m - pathLength_m) / speed;
            applyFractionalStep(pusher, tau);

            break;
        }
    }

    changeDT();
}
void ParallelTracker::dumpStats(long long step, bool psDump, bool statDump) {
    OPALTimer::Timer myt2;

    /*
    if (itsBunch_m->getGlobalTrackStep() % 1000 + 1 == 1000) {
        *gmsg << level1;
    } else if (itsBunch_m->getGlobalTrackStep() % 100 + 1 == 100) {
        *gmsg << level2;
    } else {
        *gmsg << level3;
    }
    */
    
    if (numParticlesInSimulation_m == 0) {
        *gmsg << "* " << myt2.time() << " "
            << "Step " << std::setw(6) << itsBunch_m->getGlobalTrackStep() << "; "
            << "   -- no emission yet --     "
            << "t= " << Util::getTimeString(itsBunch_m->getT()) << endl;
        return;
    }

    // \todo itsBunch_m->calcEMean();
    //    size_t totalParticles_f = numParticlesInSimulation_m;
    if (std::isnan(pathLength_m) || std::isinf(pathLength_m)) {
        throw OpalException(
            "ParallelTracker::dumpStats()",
            "there seems to be something wrong with the position of the bunch!");
    } else {
        *gmsg << "* " << myt2.time() << " "
            << "Step " << std::setw(6) << itsBunch_m->getGlobalTrackStep() << " "
            << "at " << Util::getLengthString(pathLength_m) << ", "
            << "t= " << Util::getTimeString(itsBunch_m->getT()) << ", "
            << "E=" << Util::getEnergyString(itsBunch_m->get_meanKineticEnergy()) << endl;

        writePhaseSpace(step, psDump, statDump);
    }
}

void ParallelTracker::setOptionalVariables() {

    /*
    minStepforReBin_m = Options::minStepForRebin;
    RealVariable* br =
        dynamic_cast<RealVariable*>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if (br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << level2 << "MINSTEPFORREBIN " << minStepforReBin_m << endl;
    */

    // there is no point to do repartitioning with one node
    if (ippl::Comm->size() == 1) {
        repartFreq_m = std::numeric_limits<unsigned long long>::max();
    } else {
        repartFreq_m = Options::repartFreq * 100;
        RealVariable* rep =
            dynamic_cast<RealVariable*>(OpalData::getInstance()->find("REPARTFREQ"));
        if (rep)
            repartFreq_m = static_cast<int>(rep->getReal());
        *gmsg  << "* REPARTFREQ " << repartFreq_m << endl;
    }
}

bool ParallelTracker::hasEndOfLineReached(const BoundingBox& globalBoundingBox) {
    // \todo check in IPPL 1.0 it was OpBitwiseAndAssign()
    ippl::Comm->reduce(globalEOL_m,  globalEOL_m, 1, std::logical_and<bool>());
    globalEOL_m = globalEOL_m || globalBoundingBox.isOutside(itsBunch_m->RefPartR_m);
    return globalEOL_m;
}

void ParallelTracker::setTime() {

    auto dtview  = itsBunch_m->getParticleContainer()->dt.getView();
    double newdT = itsBunch_m->getdT();

    Kokkos::parallel_for(
                         "changeDT", ippl::getRangePolicy(dtview),
                         KOKKOS_LAMBDA(const int i) {
                             dtview(i) = newdT;
                         });                     
}

void ParallelTracker::writePhaseSpace(const long long /*step*/, bool psDump, bool statDump) {
    Vector_t<double, 3> externalE, externalB;
    Vector_t<double, 3> FDext[2];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t<double, 3> rmin, rmax;
    itsBunch_m->get_bounds(rmin, rmax);

    if (psDump || statDump) {
        externalB = Vector_t<double, 3>(0.0);
        externalE = Vector_t<double, 3>(0.0);
        itsOpalBeamline_m.getFieldAt(
            itsBunch_m->RefPartR_m, itsBunch_m->RefPartP_m,
            itsBunch_m->getT() - 0.5 * itsBunch_m->getdT(), externalE, externalB);
        FDext[0] = externalB;  // \todo itsBunch_m->toLabTrafo_m.rotateFrom(externalB);
        FDext[1] =
            (externalE * Units::Vpm2MVpm);  // \todo itsBunch_m->toLabTrafo_m.rotateFrom(externalE *
                                            // Units::Vpm2MVpm);
    }

    if (statDump) {
        itsDataSink_m->dumpSDDS(itsBunch_m, FDext, -1.0);
        *gmsg << "* Wrote beam statistics." << endl;
    }

    if (psDump && (itsBunch_m->getTotalNum() > 0)) {
        // Write bunch to .h5 file.
        itsDataSink_m->dumpH5(itsBunch_m, FDext);

        /*
        // Write fields to .h5 file.
        const size_t localNum    = itsBunch_m->getLocalNum();
        double distToLastStop    = stepSizes_m.getFinalZStop() - pathLength_m;
        Vector_t<double, 3> beta = itsBunch_m->RefPartP_m / Util::getGamma(itsBunch_m->RefPartP_m);
        Vector_t<double, 3> driftPerTimeStep =
            itsBunch_m->getdT()
            * Physics::c;  // \todo  * itsBunch_m->toLabTrafo_m.rotateFrom(beta);
        bool driftToCorrectPosition =
            std::abs(distToLastStop) < 0.5 * euclidean_norm(driftPerTimeStep);
        // \todo Ppos_t stashedR;
        Vector_t<double, 3> stashedR;
        Vector_t<double, 3> stashedRefPartR;

        if (driftToCorrectPosition) {
            const double tau =
                distToLastStop / euclidean_norm(driftPerTimeStep) * itsBunch_m->getdT();

            if (localNum > 0) {

                stashedR.create(localNum);
                stashedR        = itsBunch_m->R;
                stashedRefPartR = itsBunch_m->RefPartR_m;

                for (size_t i = 0; i < localNum; ++i) {
                    itsBunch_m->R[i] +=
                        tau
                        * (Physics::c * itsBunch_m->P[i] / Util::getGamma(itsBunch_m->P[i])
                           - driftPerTimeStep / itsBunch_m->getdT());
                }

            }

            driftPerTimeStep = itsBunch_m->toLabTrafo_m.rotateTo(driftPerTimeStep);
            itsBunch_m->RefPartR_m =
                itsBunch_m->RefPartR_m + tau * driftPerTimeStep / itsBunch_m->getdT();
            CoordinateSystemTrafo update(
                tau * driftPerTimeStep / itsBunch_m->getdT(), Quaternion(1.0, 0.0, 0.0, 0.0));
            itsBunch_m->toLabTrafo_m = itsBunch_m->toLabTrafo_m * update.inverted();

            itsBunch_m->set_sPos(stepSizes_m.getFinalZStop());

            itsBunch_m->calcBeamParameters();
        }
        if (!statDump && !driftToCorrectPosition)
            itsBunch_m->calcBeamParameters();

            msg << *itsBunch_m << endl;
        
        itsDataSink_m->dumpH5(itsBunch_m, FDext);
        */
        
        /*
        
        if (driftToCorrectPosition) {
            if (localNum > 0) {
                itsBunch_m->R = stashedR;
                stashedR.destroy(localNum, 0);
            }

            itsBunch_m->RefPartR_m = stashedRefPartR;
            itsBunch_m->set_sPos(pathLength_m);

            itsBunch_m->calcBeamParameters();
        }

        */
        
        *gmsg << level2 << "* Wrote beam phase space." << endl;
    }
}
/* ========================================================================== */
/* ============================ Autophasing ================================= */
void ParallelTracker::updateRFElement(std::string elName, double maxPhase) {
    FieldList cavities       = 
        itsOpalBeamline_m.getElementByType(ElementType::RFCAVITY);
    FieldList travelingwaves = 
        itsOpalBeamline_m.getElementByType(ElementType::TRAVELINGWAVE);
    cavities.insert(cavities.end(), travelingwaves.begin(), travelingwaves.end());

    for (FieldList::iterator fit = cavities.begin(); fit != cavities.end(); ++fit) {
        if ((*fit).getElement()->getName() == elName) {
            RFCavity* element = static_cast<RFCavity*>((*fit).getElement().get());

            element->setPhasem(maxPhase);
            element->setAutophaseVeto();

            *ippl::Info << "Restored cavity phase from the h5 file. Name: " << 
                element->getName() << ", phase: " << maxPhase << " rad" << endl;
            return;
        }
    }
}
void ParallelTracker::saveCavityPhases() {
    itsDataSink_m->storeCavityInformation();
}

void ParallelTracker::restoreCavityPhases() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if (OpalData::getInstance()->hasPriorTrack() || OpalData::getInstance()->inRestartRun()) {
        iterator_t it  = OpalData::getInstance()->getFirstMaxPhases();
        iterator_t end = OpalData::getInstance()->getLastMaxPhases();
        for (; it < end; ++it) {
            updateRFElement((*it).first, (*it).second);
        }
    }
}

void ParallelTracker::autophaseCavities(const BorisPusher& pusher) {
    double t                  = itsBunch_m->getT();
    Vector_t<double, 3> nextR = itsBunch_m->RefPartR_m / (Physics::c * itsBunch_m->getdT());
    pusher.push(nextR, itsBunch_m->RefPartP_m, itsBunch_m->getdT());
    nextR *= Physics::c * itsBunch_m->getdT();

    auto elementSet = itsOpalBeamline_m.getElements(nextR);
    for (auto element : elementSet) {
        if (element->getType() == ElementType::TRAVELINGWAVE) {
            const TravelingWave* TWelement = static_cast<const TravelingWave*>(element.get());
            if (!TWelement->getAutophaseVeto()) {
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(
                    itsOpalBeamline_m.transformToLocalCS(element, itsBunch_m->RefPartR_m),
                    itsOpalBeamline_m.rotateToLocalCS(element, itsBunch_m->RefPartP_m), t,
                    itsBunch_m->getdT());
            }

        } else if (element->getType() == ElementType::RFCAVITY) {
            const RFCavity* RFelement = static_cast<const RFCavity*>(element.get());
            if (!RFelement->getAutophaseVeto()) {
                CavityAutophaser ap(itsReference, element);
                ap.getPhaseAtMaxEnergy(
                    itsOpalBeamline_m.transformToLocalCS(element, itsBunch_m->RefPartR_m),
                    itsOpalBeamline_m.rotateToLocalCS(element, itsBunch_m->RefPartP_m), t,
                    itsBunch_m->getdT());
            }
        }
    }
}

/* ========================================================================== */
/* ============================ RING FUNCTIONS ============================== */
void ParallelTracker::buildupFieldList(
    double BcParameter[], ElementType elementType, Component* elptr) {
    beamline_list::iterator sindex;

    type_pair* localpair = new type_pair();
    localpair->first     = elementType;

    for (int i = 0; i < 8; i++)
        *(((localpair->second).first) + i) = *(BcParameter + i);

    (localpair->second).second = elptr;

    // always put cyclotron as the first element in the list.
    if (elementType == ElementType::RING) {
        sindex = FieldDimensions.begin();
    } else {
        sindex = FieldDimensions.end();
    }
    FieldDimensions.insert(sindex, localpair);
}
bool ParallelTracker::applyPluginElements(const double dt) {
    IpplTimings::startTimer(PluginElemTimer_m);

    bool flag = false;
    for (PluginElement* element : pluginElements_m) {
        bool tmp = element->check(itsBunch_m, turnnumber_m, itsBunch_m->getT(), dt);
        flag |= tmp;

        if (tmp) {
            itsBunch_m->updateNumTotal();
            *gmsg << "* Total number of particles after PluginElement= "
                  << itsBunch_m->getTotalNum() << endl;
        }
    }

    IpplTimings::stopTimer(PluginElemTimer_m);
    return flag;
}
/* ========================================================================== */

struct DistributionInfo {
    unsigned int who;
    unsigned int whom;
    unsigned int howMany;
};
