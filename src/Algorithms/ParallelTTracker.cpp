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

#include "Algorithms/ParallelTTracker.h"

#include <cfloat>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

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
#include "Utilities/Options.h"
#include "Utilities/Options.h"

#include "Distribution/Distribution.h"
#include "ValueDefinitions/RealVariable.h"
#include "Utilities/Timer.h"
#include "Utilities/OpalException.h"
#include "Solvers/SurfacePhysicsHandler.hh"
#include "Structure/BoundaryGeometry.h"
#define EPS 10e-10
class PartData;

using namespace std;

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
				   size_t N):
Tracker(beamline, reference, revBeam, revTrack),
itsBunch(NULL),
itsDataSink_m(NULL),
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
zStop_m(),
scaleFactor_m(1.0),
vscaleFactor_m(scaleFactor_m),
recpGamma_m(1.0),
rescale_coeff_m(1.0),
dtCurrentTrack_m(0.0),
dtAllTracks_m(),
surfaceEmissionStop_m(-1),
specifiedNPart_m(N),
minStepforReBin_m(-1),
minBinEmitted_m(std::numeric_limits<size_t>::max()),
repartFreq_m(-1),
lastVisited_m(-1),
numRefs_m(-1),
gunSubTimeSteps_m(-1),
emissionSteps_m(std::numeric_limits<unsigned int>::max()),
localTrackSteps_m(),
maxNparts_m(0),
numberOfFieldEmittedParticles_m(std::numeric_limits<size_t>::max()),
bends_m(0),
numParticlesInSimulation_m(0),
totalParticlesInSimulation_m(0),
space_orientation_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
timeFieldEvaluation_m(IpplTimings::getTimer("Fieldeval")),
BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
Nimpact_m(0),
SeyNum_m(0.0),
timeIntegrationTimer1Push_m(IpplTimings::getTimer("TIntegration1Push")),
timeIntegrationTimer2Push_m(IpplTimings::getTimer("TIntegration2Push"))
{
}

ParallelTTracker::ParallelTTracker(const Beamline &beamline,
                                   PartBunch &bunch,
                                   DataSink &ds,
                                   const PartData &reference,
                                   bool revBeam,
                                   bool revTrack,
                                   const std::vector<unsigned long long> &maxSteps,
                                   const std::vector<double> &zstop,
                                   int timeIntegrator,
                                   const std::vector<double> &dt,
				   size_t N):
Tracker(beamline, reference, revBeam, revTrack),
itsBunch(&bunch),
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
scaleFactor_m(itsBunch->getdT() * Physics::c),
vscaleFactor_m(scaleFactor_m),
recpGamma_m(1.0),
rescale_coeff_m(1.0),
dtCurrentTrack_m(0.0),
surfaceEmissionStop_m(-1),
specifiedNPart_m(N),
minStepforReBin_m(-1),
minBinEmitted_m(std::numeric_limits<size_t>::max()),
repartFreq_m(-1),
lastVisited_m(-1),
numRefs_m(-1),
gunSubTimeSteps_m(-1),
emissionSteps_m(numeric_limits<unsigned int>::max()),
maxNparts_m(0),
numberOfFieldEmittedParticles_m(numeric_limits<size_t>::max()),
bends_m(0),
numParticlesInSimulation_m(0),
totalParticlesInSimulation_m(0),
space_orientation_m(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
timeIntegrationTimer1_m(IpplTimings::getTimer("TIntegration1")),
timeIntegrationTimer2_m(IpplTimings::getTimer("TIntegration2")),
timeFieldEvaluation_m(IpplTimings::getTimer("Fieldeval")),
BinRepartTimer_m(IpplTimings::getTimer("Binaryrepart")),
WakeFieldTimer_m(IpplTimings::getTimer("WakeField")),
timeIntegrator_m(timeIntegrator),
Nimpact_m(0),
SeyNum_m(0.0),
timeIntegrationTimer1Push_m(IpplTimings::getTimer("TIntegration1Push")),
timeIntegrationTimer2Push_m(IpplTimings::getTimer("TIntegration2Push"))
{

    for (std::vector<unsigned long long>::const_iterator it = maxSteps.begin(); it != maxSteps.end(); ++ it) {
        localTrackSteps_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = dt.begin(); it != dt.end(); ++ it) {
        dtAllTracks_m.push(*it);
    }
    for (std::vector<double>::const_iterator it = zstop.begin(); it != zstop.end(); ++ it) {
        zStop_m.push(*it);
    }

    //    itsBeamline = dynamic_cast<Beamline*>(beamline.clone());

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        dksbase.setAPI("Cuda", 4);
        dksbase.setDevice("-gpu", 4);
        dksbase.initDevice();
    }
#endif

}


ParallelTTracker::~ParallelTTracker() {

}

void ParallelTTracker::writeTrackOrbitFile(double spos) {

    /// note this is needs to be changed adapting
    /// the OPAL-cyc idea

    int locIdofPart1 = -1;
    int locIdofPart2 = -1;

    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG4, IPPL_APP_CYCLE);

    /// find the particles
    for(size_t ii = 0; ii < (itsBunch->getLocalNum()); ii++) {
        if(itsBunch->ID[ii] == 1)
            locIdofPart1 = ii;
        if(itsBunch->ID[ii] == 2)
            locIdofPart2 = ii;
    }

    if(Ippl::myNode() == 0) {
        /// write my one particle first
        if (locIdofPart1 >= 0)
            outfTrackOrbit_m << "ID" << itsBunch->ID[locIdofPart1] << "\t "
                             << spos << "\t "
                             << itsBunch->R[locIdofPart1](0) << "\t "
                             << itsBunch->R[locIdofPart1](1) << "\t "
                             << itsBunch->R[locIdofPart1](2) << "\t "
                             << itsBunch->P[locIdofPart1](0) << "\t "
                             << itsBunch->P[locIdofPart1](1) << "\t "
                             << itsBunch->P[locIdofPart1](2) << "\t "
                             << std::endl;

        if (locIdofPart2 >= 0)
            outfTrackOrbit_m << "ID" << itsBunch->ID[locIdofPart2] << "\t "
                             << spos << "\t "
                             << itsBunch->R[locIdofPart2](0) << "\t "
                             << itsBunch->R[locIdofPart2](1) << "\t "
                             << itsBunch->R[locIdofPart2](2) << "\t "
                             << itsBunch->P[locIdofPart2](0) << "\t "
                             << itsBunch->P[locIdofPart2](1) << "\t "
                             << itsBunch->P[locIdofPart2](2) << "\t "
                             << std::endl;

        outfTrackOrbit_m.flush();


        int notReceived = Ippl::getNodes() - 1;
        int dataBlocks = 0;
        Vector_t x,p;
        int id;
        while(notReceived > 0) {
            int node = COMM_ANY_NODE;
            std::unique_ptr<Message> rmsg(Ippl::Comm->receive_block(node, tag));
            if(!bool(rmsg)) {
                ERRORMSG("Could not receive from client nodes in writeTrackOrbitFile." << endl);
            }
            notReceived--;
            rmsg->get(&dataBlocks);
            for(int i = 0; i < dataBlocks; ++i) {
                rmsg->get(&id);
                rmsg->get(x);
                rmsg->get(p);

                outfTrackOrbit_m << "ID" << id << " \t "
                                 << spos << " \t "
                                 << x(0) << " \t "
                                 << x(1) << " \t "
                                 << x(2) << " \t "
                                 << p(0) << " \t "
                                 << p(1) << " \t "
                                 << p(2) << " \t "
                                 << std::endl;
                outfTrackOrbit_m.flush();
            }
        }
    }
    else {
        Message *smsg = new Message();
        int dataToSend = 0;
        if (locIdofPart1>=0)
            dataToSend++;
        if (locIdofPart2>=0)
            dataToSend++;

        smsg->put(dataToSend);

        if (locIdofPart1>=0){
            smsg->put(itsBunch->ID[locIdofPart1]);
            smsg->put(itsBunch->R[locIdofPart1]);
            smsg->put(itsBunch->P[locIdofPart1]);
        }

        if (locIdofPart2>=0){
            smsg->put(itsBunch->ID[locIdofPart2]);
            smsg->put(itsBunch->R[locIdofPart2]);
            smsg->put(itsBunch->P[locIdofPart2]);
        }
        bool res = Ippl::Comm->send(smsg, 0, tag);
        if(! res)
            ERRORMSG("Ippl::Comm->send(smsg, 0, tag) failed in writeTrackOrbitFile" << endl);
    }
}

void ParallelTTracker::initTrackOrbitFile() {

    if(Ippl::myNode() == 0) {
        std::string f = std::string("data/") + OpalData::getInstance()->getInputBasename() + std::string("-trackOrbit.dat");
        outfTrackOrbit_m.setf(std::ios::scientific, std::ios::floatfield);
        outfTrackOrbit_m.precision(8);

        if(OpalData::getInstance()->inRestartRun()) {

            outfTrackOrbit_m.open(f.c_str(), std::ios::app);
            outfTrackOrbit_m << "# Restart at integration step " << itsBunch->getLocalTrackStep() << std::endl;
        } else {

            outfTrackOrbit_m.open(f.c_str());
            outfTrackOrbit_m << "# The six-dimensional phase space data in the global Cartesian coordinates" << std::endl;
            outfTrackOrbit_m << "# Part. ID \t  SPOS [m] \t  x [m] \t  beta_x*gamma \t y [m] \t beta_y*gamma \t z [m] \t beta_z*gamma" << std::endl;
        }
    }
}

void ParallelTTracker::applyEntranceFringe(double angle, double curve,
                                           const BMultipoleField &field, double scale) {
}


void ParallelTTracker::applyExitFringe(double angle, double curve,
                                       const BMultipoleField &field, double scale) {
}

void ParallelTTracker::updateRFElement(std::string elName, double maxPhase) {
    /**
     The maximum phase is added to the nominal phase of
     the element. This is done on all nodes except node 0 where
     the Autophase took place.
     */
    double phase = 0.0;
    double frequency = 0.0;
    double globalTimeShift = OpalData::getInstance()->getGlobalPhaseShift();
    for (FieldList::iterator fit = cavities_m.begin(); fit != cavities_m.end(); ++fit) {
        if ((*fit).getElement()->getName() == elName) {
            if ((*fit).getElement()->getType() == ElementBase::TRAVELINGWAVE) {
                phase  =  static_cast<TravelingWave *>((*fit).getElement().get())->getPhasem();
                frequency = static_cast<TravelingWave *>((*fit).getElement().get())->getFrequencym();
                maxPhase -= frequency * globalTimeShift;

                static_cast<TravelingWave *>((*fit).getElement().get())->updatePhasem(phase + maxPhase);
            } else {
                phase  = static_cast<RFCavity *>((*fit).getElement().get())->getPhasem();
                frequency = static_cast<RFCavity *>((*fit).getElement().get())->getFrequencym();
                maxPhase -= frequency * globalTimeShift;

                static_cast<RFCavity *>((*fit).getElement().get())->updatePhasem(phase + maxPhase);
            }

            break;
        }
    }
}

void ParallelTTracker::handleAutoPhasing() {
    typedef std::vector<MaxPhasesT>::iterator iterator_t;

    if(Options::autoPhase == 0) return;

    if(!OpalData::getInstance()->inRestartRun()) {
        itsDataSink_m->storeCavityInformation();
    }

    iterator_t it = OpalData::getInstance()->getFirstMaxPhases();
    iterator_t end = OpalData::getInstance()->getLastMaxPhases();
    for(; it < end; ++ it) {
        updateRFElement((*it).first, (*it).second);
    }
}

void ParallelTTracker::execute() {

    initTrackOrbitFile() ;

    if(timeIntegrator_m == 3) {
        executeAMTSTracker();
    } else {
        executeDefaultTracker();
    }

    if(Ippl::myNode() == 0) outfTrackOrbit_m.close();
}

void ParallelTTracker::executeDefaultTracker() {
    Inform msg("ParallelTTracker ", *gmsg);
    const Vector_t vscaleFactor_m = Vector_t(scaleFactor_m);
    BorisPusher pusher(itsReference);
    secondaryFlg_m = false;
    dtCurrentTrack_m = itsBunch->getdT();

    // upper limit of particle number when we do field emission and secondary emission
    // simulation. Could be reset to another value in input file with MAXPARTSNUM.
    maxNparts_m = 100000000;
    nEmissionMode_m = true;

    prepareSections();

    if (OpalData::getInstance()->hasBunchAllocated() &&
        (itsBunch->getGlobalTrackStep() - 1) % Options::statDumpFreq != 0) {
        // delete last entry of sdds file and load balance file
        // if we are in a follow-up track
        itsDataSink_m->rewindLinesSDDS(1);
        itsDataSink_m->rewindLinesLBal(1);
    }

    handleAutoPhasing();

    numParticlesInSimulation_m = itsBunch->getTotalNum();
    totalParticlesInSimulation_m = itsBunch->getTotalNum();

    OPALTimer::Timer myt1;

    setTime();

    double t = itsBunch->getT();

    unsigned long long step = itsBunch->getLocalTrackStep();

    *gmsg << "Track start at: " << myt1.time() << ", t= " << t << "; zstop at: " << zStop_m.front() << " [m]" << endl;

    gunSubTimeSteps_m = 10;
    prepareEmission();

    doSchottyRenormalization();

    *gmsg << level1
          << "Executing ParallelTTracker, initial DT " << itsBunch->getdT() << " [s];\n"
          << "max integration steps " << localTrackSteps_m.front() << ", next step= " << step << "\n";
    *gmsg << "Using default (Boris-Buneman) integrator" << endl;

    itsOpalBeamline_m.print(*gmsg);

    setupSUV(!(OpalData::getInstance()->inRestartRun() || (OpalData::getInstance()->hasBunchAllocated() && !Options::scan)));

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

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        //get number of elements in the bunch
        numDeviceElements = itsBunch->getLocalNum();

        //allocate memory on device
        r_ptr = dksbase.allocateMemory<Vector_t>(numDeviceElements, ierr);
        p_ptr = dksbase.allocateMemory<Vector_t>(numDeviceElements, ierr);
        x_ptr = dksbase.allocateMemory<Vector_t>(numDeviceElements, ierr);

        lastSec_ptr = dksbase.allocateMemory<long>(numDeviceElements, ierr);
        dt_ptr = dksbase.allocateMemory<double>(numDeviceElements, ierr);

        orient_ptr = dksbase.allocateMemory<Vector_t>(itsOpalBeamline_m.sections_m.size(), ierr);

        //get all the section orientations
        int nsec = itsOpalBeamline_m.sections_m.size();
        Vector_t *orientation = new Vector_t[nsec];
        for (long i = 0; i < nsec; i++) {
            orientation[i] = itsOpalBeamline_m.getOrientation(i);
        }

        //write orientations to device
        dksbase.writeData<Vector_t>(orient_ptr, orientation, nsec);

        //free local orientation memory
        delete[] orientation;

        //allocate memory on device for particle
        allocateDeviceMemory();

        //page lock itsBunch->X, itsBunch->R, itsBunch-P
        registerHostMemory();

        //write R, P and X data to device
        dksbase.writeDataAsync<Vector_t>(r_ptr, &itsBunch->R[0], itsBunch->getLocalNum());
        dksbase.writeDataAsync<Vector_t>(p_ptr, &itsBunch->P[0], itsBunch->getLocalNum());
        dksbase.writeDataAsync<Vector_t>(x_ptr, &itsBunch->X[0], itsBunch->getLocalNum());

        //create two new streams
        dksbase.createStream(stream1);
        dksbase.createStream(stream2);
    }
#endif

    while (localTrackSteps_m.size() > 0) {
        localTrackSteps_m.front() += step;
        dtCurrentTrack_m = dtAllTracks_m.front();
        changeDT();

        for(; step < localTrackSteps_m.front(); ++step) {
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

            switchElements(10.0);

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

            bool const psDump = itsBunch->getGlobalTrackStep() % Options::psDumpFreq == 0;
            bool const statDump = itsBunch->getGlobalTrackStep() % Options::statDumpFreq == 0;

            dumpStats(step, psDump, statDump);

            if(hasEndOfLineReached()) break;

            itsBunch->incTrackSteps();

        }

        dtAllTracks_m.pop();
        localTrackSteps_m.pop();
        zStop_m.pop();
    }

    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
        numParticlesInSimulation_m = itsBunch->getTotalNum();
    }

    bool const psDump = (itsBunch->getGlobalTrackStep() - 1) % Options::psDumpFreq != 0;
    bool const statDump = (itsBunch->getGlobalTrackStep() - 1) % Options::statDumpFreq != 0;
    writePhaseSpace((step + 1), itsBunch->get_sPos(), psDump, statDump);
    msg << level2 << "Dump phase space of last step" << endl;
    OPALTimer::Timer myt3;
    itsOpalBeamline_m.switchElementsOff();

    *gmsg << "done executing ParallelTTracker at " << myt3.time() << endl;

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        //free device memory
        freeDeviceMemory();
        dksbase.freeMemory<Vector_t>(orient_ptr, itsOpalBeamline_m.sections_m.size());
        //unregister page lock itsBunch->X, itsBunch->R, itsBunch-P
        unregisterHostMemory();
    }
#endif
}

void ParallelTTracker::executeAMTSTracker() {
    Inform msg("ParallelTTracker ", *gmsg);
    const Vector_t vscaleFactor_m = Vector_t(scaleFactor_m);
    dtCurrentTrack_m = itsBunch->getdT();

    // upper limit of particle number when we do field emission and secondary emission
    // simulation. Could be reset to another value in input file with MAXPARTSNUM.
    maxNparts_m = 100000000;

    prepareSections();

    handleAutoPhasing();

    numParticlesInSimulation_m = itsBunch->getTotalNum();
    setTime();
    unsigned long long step = itsBunch->getLocalTrackStep();
    msg << "Track start at: " << OPALTimer::Timer().time() << ", t = " << itsBunch->getT() << "; zstop at: " << zStop_m.front() << " [m]" << endl;
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
    double const tEnd = itsBunch->getT() + double(localTrackSteps_m.front() - step) * dt_inner_target;
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

void ParallelTTracker::doSchottyRenormalization() {
    Inform msg("ParallelTTracker ", *gmsg);
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

double ParallelTTracker::schottkyLoop(double rescale_coeff) {

    Inform msg("ParallelTTracker ", *gmsg);

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
        itsOpalBeamline_m.switchElements(rmin(2) - margin, rmax(2) + margin, getEnergyMeV(RefPartP_suv_m));
    }

    double minBinEmitted  = 10.0;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar) {
        minBinEmitted = ar->getReal();  // the space charge solver crashes if we use less than ~10 particles.
        // This variable controls the number of particles to be emitted before we use
        // the space charge solver.

        msg << level3 << "MINBINEMITTED " << minBinEmitted << endl;
    }


    double minStepforReBin  = 10000.0;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br) {
        minStepforReBin = br->getReal();  // this variable controls the minimal number of steps of emission (using bins)
        // before we can merge the bins
        msg << level3 << "MINSTEPFORREBIN " << minStepforReBin << endl;
    }

    int repartFreq = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep) {
        repartFreq = static_cast<int>(rep->getReal());  // this variable controls the minimal number of steps until we repartition the particles
        msg << level3 << "REPARTFREQ " << repartFreq << endl;
    }

    // there is no point to do repartitioning with one node
    if(Ippl::getNodes() == 1)
        repartFreq = 1000000;

    size_t totalParticles_f = 0;

    for(; step < localTrackSteps_m.front(); ++step) {
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

        // if(itsBunch->getLocalNum() == 0)
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
            msg << "* Deleted in Shotky " << ne << " particles, remaining "
                << totalParticles_f << " particles" << endl; //benchmark output

        kickParticles(pusher);

        if(totalParticles_f > 0) {
            // none of the particles is in a bending element
            updateReferenceParticle();
        }

        itsBunch->RefPart_R = RefPartR_zxy_m;
        itsBunch->RefPart_P = RefPartP_zxy_m;

        /*
          calculate the dimensions of the bunch and add a small margin to them;
          then decide which elements have to be triggered when an element is
          triggered memory is allocated and the field map is read in
        */
        itsBunch->get_bounds(rmin, rmax);

        // trigger the elements
        margin = 3. * RefPartP_suv_m(2) * recpgamma;
        margin = 0.01 > margin ? 0.01 : margin;
        itsOpalBeamline_m.switchElements(
            (rmin(2) - margin)*scaleFactor_m,
            (rmax(2) + margin)*scaleFactor_m,
            getEnergyMeV(RefPartP_suv_m));

        // start normal particle loop part 2 for simulation without boundary geometry.
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            //and scale back to dimensions
            itsBunch->R[i] *= Vector_t (
                Physics::c * itsBunch->dt[i],
                Physics::c * itsBunch->dt[i],
                Physics::c * itsBunch->dt[i]);
            // update local coordinate system
            itsBunch->X[i] /= vscaleFactor;
            pusher.push(
                itsBunch->X[i],
                TransformTo (
                    itsBunch->P[i],
                    itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
                itsBunch->getdT());
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
                INFOMSG(myt2.time() << " Step " << step
                        << "; there seems to be something wrong with the position of the bunch!" << endl);
            } else {
                INFOMSG(myt2.time() << " Step " << step
                        << " at " << sposRef << " [m] t= " << t << " [s] E=" << itsBunch->get_meanEnergy() << " [MeV] " << endl);
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
            if(sposRef > zStop_m.front()) {
                localTrackSteps_m.front() = step;
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


void ParallelTTracker::applySchottkyCorrection(PartBunch &itsBunch, int ne, double t, double rescale_coeff) {

    Inform msg("ParallelTTracker ", *gmsg);
    const long ls = 0;
    /*
     Now I can calculate E_{rf} at each position
     of the newely generated particles and rescale Q

     Note:
     For now I only sample the field of the last emitted particles.
     Space charge is not yet included
     */

    double laser_erg = itsBunch.getLaserEnergy();
    double workFunction = itsBunch.getWorkFunctionRf();
    // coeffecient for calculate schottky potenial from E field [eV/(MV^0.5)]
    const double schottky_coeff = 0.037947;

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

        /*
          fabs(Ez): if the field of cathode surface is in the right direction,
          it will increase the energy which electron obtain. If the field is in
          the wrong direction, this particle will be back to the cathode surface
          and then be deleted automaticly by OPAL,  we don't add another logical
          branch to handle this. So fabs is the simplest way to handle this
        */
        obtain_erg = laser_erg - workFunction + schottky_coeff * sqrt(fabs(Ez) / 1E6);
        double schottkyScale = obtain_erg * obtain_erg * rescale_coeff;

        itsBunch.Q[n] *= schottkyScale;
        itsBunch.R[n] /= Vector_t(Physics::c * itsBunch.dt[n]);
    }
}

void ParallelTTracker::bgf_main_collision_test() {
    if(!bgf_m) return;

    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();
    /**
     Here we check if a particles is
     outside the geometry, flag it for
     deletion and create secondaries
     */
    size_t Inc_num = itsBunch->getLocalNum();
    for(size_t i = 0; i < Inc_num; i++) {
        double dtime = 0.5 * itsBunch->getdT();
        double seyNum = 0;
        if(secondaryFlg_m != 0) {//chuan: only for the secondary emission, change the particle type to be the old secondaries.
            if(itsBunch->PType[i] == ParticleType::NEWSECONDARY)
                itsBunch->PType[i] = ParticleType::SECONDARY;
        }

        /*
          for primary bunch, primary dark current particles,
          old secondaries in previous time steps and newly
          generated secondaries which have no collision with
          boundary in both first and second half step, do main
          collision test and emit the secondaries.
        */
        int res = 0;
        int triId = itsBunch->TriID[i];
        Vector_t position = itsBunch->R[i];
        if (triId == 0) {
            // no previous collision
            Vector_t intersection_pt = outr;
            res = bgf_m->PartInside (
                itsBunch->R[i],
                itsBunch->P[i],
                dtime,
                itsBunch->PType[i],
                itsBunch->Q[i],
                intersection_pt,
                itsBunch->TriID[i]);//Chuan: should be itsBunch->TriID[i];
            if (res < 0) continue;
            position = intersection_pt;
        }

        /*
          Attention: New secondaries have not been kicked and are without new momentum.
        */
        if (secondaryFlg_m == 1) {
            res += bgf_m->emitSecondaryFurmanPivi (
                position,
                i,
                itsBunch, seyNum);

        } else if (secondaryFlg_m != 0) {
            res += bgf_m->emitSecondaryVaughan (
                position,
                i,
                itsBunch, seyNum);

        } else {
            res += bgf_m->emitSecondaryNone (
                position, triId);
        }
        if(res >= 0) {
            itsBunch->Bin[i] = -1;
            Nimpact_m++;
            SeyNum_m += seyNum;
        }
    } // end for()

    for(size_t i = 0; i < itsBunch->getLocalNum(); i++) {
        if (bgf_m->isOutsideApperture(itsBunch->R[i])) {
            itsBunch->Bin[i] = -1;
        }
    }

    // Now we do field emission
    if(itsBunch->getT() < surfaceEmissionStop_m)
        numberOfFieldEmittedParticles_m +=
            bgf_m->doFNemission (itsOpalBeamline_m, itsBunch, itsBunch->getT());

    itsBunch->boundp();
    numParticlesInSimulation_m = itsBunch->getTotalNum();
}

void ParallelTTracker::handleOverlappingMonitors() {
    // make sure that no monitor has overlap with two tracks
    FieldList monitors = itsOpalBeamline_m.getElementByType(ElementBase::MONITOR);
    for(FieldList::iterator it = monitors.begin(); it != monitors.end(); ++ it) {
        double zbegin, zend;
        it->getElement()->getDimensions(zbegin, zend);
        if(zbegin < zStop_m.front() && zend >= zStop_m.back()) {
            ERRORMSG(
                "\033[0;31m"
                << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
                << "% Removing '" << it->getElement()->getName() << "' since it resides in two tracks.   %\n"
                << "% Please adjust zstop or place your monitor at a different position to prevent this. %\n "
                << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"
                << "\033[0m"
                << endl);
            static_cast<Monitor *>(it->getElement().get())->moveBy(-zend - 0.001);
            itsOpalBeamline_m.removeElement(it->getElement()->getName());
        }
    }
}

void ParallelTTracker::prepareSections() {

    itsBeamline_m.accept(*this);
    handleOverlappingMonitors();
    itsOpalBeamline_m.prepareSections();

    cavities_m = itsOpalBeamline_m.getElementByType(ElementBase::RFCAVITY);
    FieldList travelingwaves = itsOpalBeamline_m.getElementByType(ElementBase::TRAVELINGWAVE);
    cavities_m.merge(travelingwaves, ClassicField::SortAsc);
}

void ParallelTTracker::timeIntegration1(BorisPusher & pusher) {

    IpplTimings::startTimer(timeIntegrationTimer1_m);
    if(bgf_m != NULL && secondaryFlg_m > 0) return;

    IpplTimings::startTimer(timeIntegrationTimer1Push_m);

    if (IpplInfo::DKSEnabled) {
#ifdef OPAL_DKS
        //if bunch is largen than before reallocate memory
        if (itsBunch->getLocalNum() > numDeviceElements) {
            freeDeviceMemory();
            unregisterHostMemory();
            allocateDeviceMemory();
            registerHostMemory();
        }

        //get all the section orientations
        int nsec = itsOpalBeamline_m.sections_m.size();
        Vector_t *orientation = new Vector_t[nsec];
        for (long i = 0; i < nsec; i++) {
            orientation[i] = itsOpalBeamline_m.getOrientation(i);
        }

        //write orientations to device
        dksbase.writeData<Vector_t>(orient_ptr, orientation, nsec);

        //free local orientation memory
        delete[] orientation;

        //write R to device
        dksbase.writeDataAsync<Vector_t>(r_ptr, &itsBunch->R[0],
                                         itsBunch->getLocalNum(), stream1);
        //write P to device
        dksbase.writeDataAsync<Vector_t>(p_ptr, &itsBunch->P[0],
                                         itsBunch->getLocalNum(), stream1);
        //write X to device
        dksbase.writeDataAsync<Vector_t>(x_ptr, &itsBunch->X[0],
                                         itsBunch->getLocalNum(), stream2);

        //write dt to device
        dksbase.writeDataAsync<double>(dt_ptr, &itsBunch->dt[0],
                                       itsBunch->getLocalNum(), stream1);

        //calc push
        dksbase.callParallelTTrackerPush (r_ptr, p_ptr, itsBunch->getLocalNum(), dt_ptr,
                                          itsBunch->getdT(), Physics::c, true, stream1);
        //read R from device
        dksbase.readDataAsync<Vector_t> (r_ptr, &itsBunch->R[0],
                                         itsBunch->getLocalNum(), stream1);

        //write LastSection to device
        dksbase.writeDataAsync<long> (lastSec_ptr, &itsBunch->LastSection[0],
                                      itsBunch->getLocalNum(), stream2);

        //sync the device to guarantee that p_ptr has been transfered
        dksbase.syncDevice();

        //calc push
        dksbase.callParallelTTrackerPushTransform(x_ptr, p_ptr, lastSec_ptr, orient_ptr,
                                                  itsBunch->getLocalNum(),
                                                  itsOpalBeamline_m.sections_m.size(),
                                                  dt_ptr, itsBunch->getdT(), Physics::c,
                                                  false, stream2);
        //read R from device
        dksbase.readDataAsync<Vector_t>(x_ptr, &itsBunch->X[0], itsBunch->getLocalNum(), stream2);

        //sync and wait till all tasks and reads are complete
        dksbase.syncDevice();
#endif
    } else {
        itsBunch->switchToUnitlessPositions();
        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
            //scale each particle with c*dt
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);

            // update local coordinate system of particleInform &PartBunc
            pusher.push(
                itsBunch->X[i],
                TransformTo(
                    itsBunch->P[i],
                    itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
                itsBunch->getdT());
        }
        itsBunch->switchOffUnitlessPositions();
    }
    IpplTimings::stopTimer(timeIntegrationTimer1Push_m);

    if(numParticlesInSimulation_m > minBinEmitted_m) {
        itsBunch->boundp();
    }
    IpplTimings::stopTimer(timeIntegrationTimer1_m);
}

void ParallelTTracker::timeIntegration1_bgf(BorisPusher & pusher) {
    if(bgf_m == NULL || secondaryFlg_m == 0) return;
    IpplTimings::startTimer(timeIntegrationTimer1_m);

    /*
      We do collision test for newly generated secondaries before integration
      in the first half step of each time step.  This is because only secondary
      emission model yield non zero inital momenta. The initial momenta of
      field emitted particles are zero.  If hit, we set itsBunch->R[i] to
      intersection points, else we do normal integration.
    */
    Nimpact_m = 0;  // Initial parallel plate benchmark variable.
    SeyNum_m = 0;   // Initial parallel plate benchmark variable.

    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();
    double dt = itsBunch->getdT();
    double bgf_scaleFactor = dt * Physics::c;
    Vector_t bgf_vscaleFactor = Vector_t(bgf_scaleFactor);

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        bool particleHitBoundary = false;
        Vector_t intecoords = outr;
        int triId = 0;
        if (itsBunch->PType[i] == ParticleType::NEWSECONDARY) {
            particleHitBoundary = bgf_m->PartInside(
                itsBunch->R[i],
                itsBunch->P[i],
                0.5 * itsBunch->dt[i],
                itsBunch->PType[i],
                itsBunch->Q[i],
                intecoords,
                triId) == 0;
        }
        Vector_t rtmp = itsBunch->R[i];
        if(particleHitBoundary && (std::fabs(rtmp[0]-intecoords[0])>EPS) && (std::fabs(rtmp[1]-intecoords[1])>EPS) && (std::fabs(rtmp[2]-intecoords[2])>EPS)) {
            /*
              set particle position to intersection points coordinates and
              scale the position;
              no scaling required
            */
            itsBunch->R[i] = intecoords/bgf_vscaleFactor;
            itsBunch->TriID[i] = triId;
        } else {
            itsBunch->R[i] /= bgf_vscaleFactor;
            pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
        }
        /*
          :FIXME: is the local update necessary here?
          update local coordinate system for particle
        */
        itsBunch->X[i] /= bgf_vscaleFactor;
        pusher.push (
            itsBunch->X[i],
            TransformTo(
                itsBunch->P[i],
                itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
            itsBunch->getdT());

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
    // push the reference particle by a half step
    double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;
    kickParticles(pusher);
    handleBends();

    //switchElements();
    IpplTimings::startTimer(timeIntegrationTimer2Push_m);
    if (IpplInfo::DKSEnabled) {
#ifdef OPAL_DKS

        //if bunch is largen than before reallocate memory
        if (itsBunch->getLocalNum() > numDeviceElements) {
            freeDeviceMemory();
            unregisterHostMemory();
            allocateDeviceMemory();
            registerHostMemory();
        }

        //get all the section orientations
        int nsec = itsOpalBeamline_m.sections_m.size();
        Vector_t *orientation = new Vector_t[nsec];
        for (long i = 0; i < nsec; i++) {
            orientation[i] = itsOpalBeamline_m.getOrientation(i);
        }

        //write orientations to device
        dksbase.writeData<Vector_t>(orient_ptr, orientation, nsec);

        //free local orientation memory
        delete[] orientation;

        //write R to device
        dksbase.writeDataAsync<Vector_t>(r_ptr, &itsBunch->R[0],
                                         itsBunch->getLocalNum(), stream1);
        //write P to device
        dksbase.writeDataAsync<Vector_t>(p_ptr, &itsBunch->P[0],
                                         itsBunch->getLocalNum(), stream1);
        //write X to device
        dksbase.writeDataAsync<Vector_t>(x_ptr, &itsBunch->X[0],
                                         itsBunch->getLocalNum(), stream2);

        //wrote dt to device
        dksbase.writeDataAsync<double>(dt_ptr, &itsBunch->dt[0],
                                       itsBunch->getLocalNum(), stream1);
        //calc push
        dksbase.callParallelTTrackerPush(r_ptr, p_ptr, itsBunch->getLocalNum(), dt_ptr,
                                         itsBunch->getdT(), Physics::c, true, stream1);
        //read R from device
        dksbase.readDataAsync<Vector_t>(r_ptr, &itsBunch->R[0], itsBunch->getLocalNum(), stream1);

        //write LastSection to device
        dksbase.writeDataAsync<long>(lastSec_ptr, &itsBunch->LastSection[0],
                                     itsBunch->getLocalNum(), stream2);

        //sync the device to guarantee that p_ptr has been transfered
        dksbase.syncDevice();

        //calc push
        dksbase.callParallelTTrackerPushTransform(x_ptr, p_ptr, lastSec_ptr, orient_ptr,
                                                  itsBunch->getLocalNum(),
                                                  itsOpalBeamline_m.sections_m.size(),
                                                  dt_ptr, itsBunch->getdT(), Physics::c,
                                                  false, stream2);
        //read R from device
        dksbase.readDataAsync<Vector_t>(x_ptr, &itsBunch->X[0], itsBunch->getLocalNum(), stream2);

        //sync and wait till all tasks and reads are complete
        dksbase.syncDevice();
#endif
    } else {
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
            pusher.push (
                itsBunch->X[i],
                TransformTo(
                    itsBunch->P[i],
                    itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
                itsBunch->getdT());
            //reset time step if particle was emitted in the first half-step
            //the particle is now in sync with the simulation timestep

        }
        itsBunch->switchOffUnitlessPositions(true);
    }
    IpplTimings::stopTimer(timeIntegrationTimer2Push_m);

    //std::fill(itsBunch->dt.begin(), itsBunch->dt.end(), itsBunch->getdT());
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->dt[i] = itsBunch->getdT();
    }
    IpplTimings::stopTimer(timeIntegrationTimer2_m);

}

void ParallelTTracker::timeIntegration2_bgf(BorisPusher & pusher) {

    if(!bgf_m) return;
    /*
      After kick, we do collision test before integration in second half step
      with new momentum, if hit, then move collision particles to the position
      where collision occurs.
    */
    IpplTimings::startTimer(timeIntegrationTimer2_m);

    // push the reference particle by a half step
    double recpgamma = 1.0 / sqrt(1.0 + dot(RefPartP_suv_m, RefPartP_suv_m));
    RefPartR_zxy_m += RefPartP_zxy_m * recpgamma / 2. * scaleFactor_m;
    kickParticles(pusher, 0);
    handleBends();

    const Vector_t outr = bgf_m->getmaxcoords() + bgf_m->gethr();
    double dtime = 0.5 * itsBunch->getdT();

    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        bool particleHitBoundary = false;
        Vector_t intecoords = outr;
        int triId = 0;

        /*
          test all particles except those already have collided the boundary in the first half step.
          particles which already have collided the boundary in the first half step will stay at
          the positions that they collide on the boundaries.
         */
        if(itsBunch->TriID[i] == 0) {
            particleHitBoundary =  bgf_m->PartInside(
                itsBunch->R[i],
                itsBunch->P[i],
                dtime,
                itsBunch->PType[i],
                itsBunch->Q[i],
                intecoords,
                triId) == 0;
            if(particleHitBoundary) {
                itsBunch->R[i] = intecoords / Vector_t (Physics::c * itsBunch->dt[i]);
                itsBunch->TriID[i] = triId;
            } else {
                //if no collision do normal push in the second half-step
                itsBunch->R[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
                pusher.push(itsBunch->R[i], itsBunch->P[i], itsBunch->dt[i]);
            }
            itsBunch->X[i] /= Vector_t(Physics::c * itsBunch->dt[i]);
            pusher.push(
                itsBunch->X[i],
                TransformTo(
                    itsBunch->P[i],
                    itsOpalBeamline_m.getOrientation(itsBunch->LastSection[i])),
                itsBunch->getdT());
            itsBunch->R[i] *= Vector_t(Physics::c * itsBunch->dt[i]);
            itsBunch->X[i] *= Vector_t(Physics::c * itsBunch->dt[i]);
        }
    }

    //fill(itsBunch->dt.begin(), itsBunch->dt.end(), itsBunch->getdT());
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->dt[i] = itsBunch->getdT();
    }

    IpplTimings::stopTimer(timeIntegrationTimer2_m);
}

void ParallelTTracker::selectDT() {

    if(itsBunch->GetIfBeamEmitting()) {
        double dt = itsBunch->GetEmissionDeltaT();
        itsBunch->setdT(dt);
        scaleFactor_m = dt * Physics::c;
        vscaleFactor_m = Vector_t(scaleFactor_m);
    } else {
        double dt = dtCurrentTrack_m;
        itsBunch->setdT(dt);
        scaleFactor_m = dt * Physics::c;
        vscaleFactor_m = Vector_t(scaleFactor_m);
    }
}

void ParallelTTracker::changeDT() {
    selectDT();
    for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
        itsBunch->dt[i] = itsBunch->getdT();
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
            itsBunch->setGlobalMeanR(itsBunch->get_centroid());
            itsBunch->computeSelfFields(binNumber);
            itsBunch->Q = Q_back;

        }

    } else {
        itsBunch->setGlobalMeanR(itsBunch->get_centroid());
        itsBunch->computeSelfFields();
    }

    // itsBunch->switchOffUnitlessPositions(true);
    // // FIXME! (see above)
    // itsBunch->boundp();
}


void ParallelTTracker::computeExternalFields() {
    IpplTimings::startTimer(timeFieldEvaluation_m);
    Inform msg("ParallelTTracker ", *gmsg);

    unsigned long hasWake = 0;
    int hasSurfacePhysics = 0;
    long wfSection = 0;
    std::set<long> sphysSections;

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
            sphysSections.insert(ls);
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

    //if any new surfacephysics sections enabled, check if there are other surface physics
    //elements fallowing directy after it, then handle these elements as well, to process particles
    //moving from one degrader to another
    std::set<long> attachedSphysSections;
    for (long it: sphysSections) {

        //get the section it + 1 if it has surface physics
        long newit = it + 1;
        while ( itsOpalBeamline_m.getSection(newit).hasSurfacePhysics() ) {

            //check if two sections following one another dont share the same surface physics handler
            if (itsOpalBeamline_m.getSection(newit).getSurfacePhysicsHandler() !=
                itsOpalBeamline_m.getSection(newit - 1).getSurfacePhysicsHandler())
            {
                attachedSphysSections.insert(newit);
            }
            newit++;
        }

    }
    //merge two sets
    sphysSections.insert(attachedSphysSections.begin(), attachedSphysSections.end());

    IpplTimings::stopTimer(timeFieldEvaluation_m);

    reduce(hasWake, hasWake, OpAddAssign());
    reduce(bends_m, bends_m, OpAddAssign());

    hasSurfacePhysics = sphysSections.size();
    reduce(hasSurfacePhysics, hasSurfacePhysics, OpMaxAssign());

    if(hasWake > 0) {
        IpplTimings::startTimer(WakeFieldTimer_m);
        reduce(wfSection, wfSection, OpMaxAssign());
        WakeFunction *wf = itsOpalBeamline_m.getWakeFunction(wfSection);
        /*--------- Added by Xiaoying Pang 04/22/2014 ---------------
         * If the CSR is turned on for a dipole, save its pointer to the CSRWakeFunction
         * and reuse it in the following drift.*/
        // xpang: start
        std::shared_ptr<const ElementBase> element = itsOpalBeamline_m.getWakeFunctionOwner(wfSection);
        if(dynamic_cast<CSRWakeFunction*>(wf))
            {
                if(dynamic_cast<const RBend*>(element.get()) || dynamic_cast<const SBend*>(element.get()))
                    wakeFunction_m = wf;
                if(dynamic_cast<const Drift*>(element.get()))
                    wf = wakeFunction_m;
            }
        // xpang: end
        if(!wakeStatus_m) {
            msg << level2 << "============== START WAKE CALCULATION =============" << endl;
            std::shared_ptr<const ElementBase> element = itsOpalBeamline_m.getWakeFunctionOwner(wfSection);
            wf->initialize(element.get());
            wakeStatus_m = true;
        }
        if(wf == NULL) {
            INFOMSG("no wakefunction attached" << endl);
        } else {
            wf->apply(*itsBunch);
        }
        IpplTimings::stopTimer(WakeFieldTimer_m);

    } else if(wakeStatus_m) {
        msg << level2 << "=============== END WAKE CALCULATION ==============" << endl;
        wakeStatus_m = false;
    }

    if(hasSurfacePhysics > 0) {    // in a section we have an element with surface physics

        // this is a workaround fighting the zero partilce problem on a core
        // noy yet understood
        // with setMimumNumberOfParticlesPerCore(2) degrader simulations for the PSI gantry are running
        itsBunch->setMimumNumberOfParticlesPerCore(0);

        std::set<long> old_sphysSections;
        std::vector<long> leftBehindSections, newSections;
        for (auto it: sphys_m) {
            old_sphysSections.insert(it.first);
        }

        std::vector<long> sectionsWithSurfacePhysics(Ippl::getNodes() * hasSurfacePhysics, -1);
        unsigned int insPosition = Ippl::myNode() * hasSurfacePhysics;
        for (auto it :sphysSections)
            sectionsWithSurfacePhysics[insPosition ++] = it;

        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                      &sectionsWithSurfacePhysics[0], hasSurfacePhysics, MPI_LONG,
                      Ippl::getComm());
        sphysSections.insert(sectionsWithSurfacePhysics.begin(),
                             sectionsWithSurfacePhysics.end());
        auto dummySection = sphysSections.find(-1);
        if (dummySection != sphysSections.end()) {
            sphysSections.erase(dummySection);
        }

        newSections.resize(leftBehindSections.size());
        {
            leftBehindSections.resize(std::max(old_sphysSections.size(),
                                               sphysSections.size()));
            auto last = std::set_difference(old_sphysSections.begin(), old_sphysSections.end(),
                                            sphysSections.begin(), sphysSections.end(),
                                            leftBehindSections.begin());
            leftBehindSections.resize(last - leftBehindSections.begin());
        }

        for (long it: leftBehindSections) {
            if (!sphys_m[it]->stillActive()) {
                sphys_m[it] = NULL;
                sphys_m.erase(it);
            }
        }

        {
            newSections.resize(std::max(old_sphysSections.size(),
                                        sphysSections.size()));
            auto last = std::set_difference(sphysSections.begin(), sphysSections.end(),
                                            old_sphysSections.begin(), old_sphysSections.end(),
                                            newSections.begin());
            newSections.resize(last - newSections.begin());
        }

        for (long it: newSections) {
            sphys_m.insert(std::make_pair(it, itsOpalBeamline_m.getSurfacePhysicsHandler(it)));
        }

        if(!surfaceStatus_m) {
            msg << level2 << "============== START SURFACE PHYSICS CALCULATION =============" << endl;
            surfaceStatus_m = true;
        }

    }

    if (surfaceStatus_m) {
        do {
            ///all particles in material if max per node is 2 and other degraders have 0 particles
            //check if more than one degrader has particles
            long degraderWithParticles;
            int degradersWithParticlesCount = 0;
            for (auto it: sphys_m) {
                it.second->AllParticlesIn(false);
                if (it.second->getParticlesInMat() > 0) {
                    degraderWithParticles = it.first;
                    degradersWithParticlesCount++;
                }
            }

            //if max particles per node is 2, and only one degrader has particles set
            //AllParticlesIn for this degrader to true
            int maxPerNode = itsBunch->getLocalNum();
            reduce(maxPerNode, maxPerNode, OpMaxAssign());
            bool allParticlesInMat = ( (unsigned)maxPerNode <= itsBunch->getMinimumNumberOfParticlesPerCore() && degradersWithParticlesCount == 1 );

            if (allParticlesInMat) {
                sphys_m[degraderWithParticles]->AllParticlesIn(true);
                msg << "All particles in degrader" << endl;
            }

            unsigned redifusedParticles = 0;
            for (auto it: sphys_m) {
                it.second->apply(*itsBunch, totalParticlesInSimulation_m);
                it.second->print(msg);

                redifusedParticles += it.second->getRedifused();
                //if all particles where in material update time to time in degrader
                if (allParticlesInMat)
                    itsBunch->setT(it.second->getTime() - itsBunch->getdT());
            }

            //perform boundp only if there are particles in the bunch, or there are particles
            //comming out of the degrader
            if ((unsigned)maxPerNode > itsBunch->getMinimumNumberOfParticlesPerCore() ||
                redifusedParticles > 0)
            {
                itsBunch->boundp();
            }

            //if bunch has no particles update time to time in degrader
            if (itsBunch->getTotalNum() == 0)
                itsBunch->setT(itsBunch->getT() + itsBunch->getdT());

        } while (itsBunch->getTotalNum() == 0);


        if (sphys_m.size() == 0) {
            msg << level2 << "============== END SURFACE PHYSICS CALCULATION =============" << endl;
            surfaceStatus_m = false;
            itsBunch->setMimumNumberOfParticlesPerCore(0);
        }
    }

    size_t ne = 0;
    bool globPartOutOfBounds = (min(itsBunch->Bin) < 0) && (itsBunch->getTotalNum() > itsBunch->getMinimumNumberOfParticlesPerCore());
    if(globPartOutOfBounds) {
        ne = itsBunch->boundp_destroyT();
        numParticlesInSimulation_m  = itsBunch->getTotalNum();
    }

    if(numParticlesInSimulation_m > minBinEmitted_m || itsBunch->getTotalNum() > minBinEmitted_m) {
        itsBunch->update();
        numParticlesInSimulation_m = itsBunch->getTotalNum();
    }

    //remove surface physics elements from start of the list till hit one which is still active
    if (surfaceStatus_m) {
        for (auto it: sphys_m) {
            if (it.second->stillActive()) {
                break;
            } else {
                //check if degrader is empty and bunch has moved past it (needed to free GPU memory)
                if ( !it.second->stillAlive(*itsBunch) )
                    msg << "Surface physiscs at element " << it.second->getName() << " has become inactive" << endl;

                it.second = NULL;
                sphys_m.erase(it.first);
            }
        }
    }

    if(ne > 0 || surfaceStatus_m) {
        msg << level1 << "* Deleted " << ne << " particles, "
          << "remaining " << itsBunch->getTotalNum() << " particles" << endl;
        totalParticlesInSimulation_m -= ne;
    }

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
    double margin = scaleMargin * RefPartP_suv_m(2) * recpgamma * Physics::c * itsBunch->getdT();
    margin = 0.01 > margin ? 0.01 : margin;
    itsOpalBeamline_m.switchElements((rmin(2) - margin), (rmax(2) + margin), getEnergyMeV(RefPartP_suv_m));
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
    Inform msg("ParallelTTracker ", *gmsg);

    if ((itsBunch->getGlobalTrackStep() + 1) % 1000 == 0) {
        msg << level1;
    } else if ((itsBunch->getGlobalTrackStep() + 1) % 100 == 0) {
        msg << level2;
    } else {
        msg << level3;
    }

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
        itsDataSink_m->writeImpactStatistics(
            *itsBunch, ustep, Nimpact_m, SeyNum_m, numberOfFieldEmittedParticles_m, nEmissionMode_m, std::string("data/PartStatistics"));

        if(((Options::surfDumpFreq) > 0) && ((step % Options::surfDumpFreq) == 0)) {
            itsDataSink_m->writeSurfaceInteraction(*itsBunch, ustep, *bgf_m, std::string("SurfaceInteraction"));
        }

        // If we are dealing with field emission and secondary emission, set upper
        // limit of particle number in simulation to prevent memory overflow.
        if(numParticlesInSimulation_m > maxNparts_m)
            localTrackSteps_m.front() = step;

        // ada reset Nimpact_m, does not make sense to integrate this we obtain a rediculus large number !!
        Nimpact_m = 0;
    }

    if(sposRef > zStop_m.front())
        localTrackSteps_m.front() = step;

}


void ParallelTTracker::setOptionalVariables() {
    Inform msg("ParallelTTracker ", *gmsg);

    minBinEmitted_m  = 10;
    RealVariable *ar = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINBINEMITTED"));
    if(ar)
        minBinEmitted_m = static_cast<size_t>(ar->getReal());
    msg << level2 << "MINBINEMITTED " << minBinEmitted_m << endl;

    minStepforReBin_m  = 200;
    RealVariable *br = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("MINSTEPFORREBIN"));
    if(br)
        minStepforReBin_m = static_cast<int>(br->getReal());
    msg << level2 << "MINSTEPFORREBIN " << minStepforReBin_m << endl;

    surfaceEmissionStop_m  = 1000.0;
    RealVariable *cr = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("SURFACEEMISSIONSTOP"));
    if(cr)
        surfaceEmissionStop_m = static_cast<double>(cr->getReal());
    msg << level2 << "SURFACEEMISSIONSTOP after " << surfaceEmissionStop_m << " seconds" <<  endl;

    repartFreq_m = 1000;
    RealVariable *rep = dynamic_cast<RealVariable *>(OpalData::getInstance()->find("REPARTFREQ"));
    if(rep)
        repartFreq_m = static_cast<int>(rep->getReal());
    msg << level2 << "REPARTFREQ " << repartFreq_m << endl;

}


bool ParallelTTracker::hasEndOfLineReached() {
    reduce(&globalEOL_m, &globalEOL_m + 1, &globalEOL_m, OpBitwiseAndAssign());
    return globalEOL_m;
}

void ParallelTTracker::setupSUV(bool updateReference) {

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

    if (updateReference) {
        RefPartP_zxy_m = RefPartP_suv_m;

        updateSpaceOrientation(false);
        RefPartP_suv_m = itsBunch->get_pmean();
    }
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
    Inform msg("ParallelTTracker ", *gmsg);

    if(mpacflg_m || !itsBunch->doEmission()) return;

    emissionSteps_m = itsBunch->GetNumberOfEmissionSteps();
    msg << level2 << "Do emission for " << itsBunch->getTEmission() << " [s] using "
        << itsBunch->GetNumberOfEnergyBins() << " energy bins " << endl
        << "Change dT from " <<  itsBunch->getdT() << " [s] to "
        <<  itsBunch->GetEmissionDeltaT() << " [s] during emission " << endl;;

}

void ParallelTTracker::initializeBoundaryGeometry() {

    Inform msg("ParallelTTracker ", *gmsg);

    // for the time being, the Boundary geomentry must be attachedto the first element
    bgf_m = itsOpalBeamline_m.getBoundaryGeometry(0);
    if (bgf_m) {
        Distribution *dist = NULL;
        Distribution *distrand = NULL;
        std::vector<std::string> distr_str = bgf_m->getDistributionArray();

        if(distr_str.size() == 0) {
            std::string distr = bgf_m->getDistribution();
            if(!distr.empty()) {
                msg << "* Find boundary geometry, start at: " << bgf_m->getS()
                    << " (m) Distribution= " << bgf_m->getDistribution() << endl;
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
                    /*
                      here nbparts should be non zero as these particles will
                      be the initialization of primary bunch.
                    */
                    size_t nbparts = distrand->GetNumberOfDarkCurrentParticles();
                    double darkinwardmargin = distrand->GetDarkCurrentParticlesInwardMargin();
                    double einitthreshold = distrand->GetEInitThreshold();
                    /*
                      make sure that the elements have already been switched
                      on before initializing particles in position where the
                      electric field > einitthreshold.
                    */
                    bgf_m->setEInitThreshold(einitthreshold);
                    if(!mpacflg_m) {
                        bgf_m->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
                        distrand->CreatePriPart(itsBunch, *bgf_m);
                        numParticlesInSimulation_m = itsBunch->getTotalNum();
                    } else {
                        /*
                          Multipacting flag set true. Generate primary particles.
                          Activate all elements (switch on the field map of elements
                          in multipacting) in multipacting simulation
                        */

                        itsOpalBeamline_m.switchAllElements();
                        /*
                          it is possible to generate initial particles according
                          to E field, since all elements switched on before we
                          create particles.
                        */
                        bgf_m->createPriPart(nbparts, darkinwardmargin, itsOpalBeamline_m, itsBunch);
                        /*
                          for Parallel Plate benchmark, Vw should be defined in
                          input file and will be invoked by getVw method in
                          createPriPart().  For other multipacting simulation
                          no need to define the Vw in SURFACERANDCREATE in input
                          file.
                        */
                        distrand->CreatePriPart(itsBunch, *bgf_m);
                        numParticlesInSimulation_m = itsBunch->getTotalNum();
                        itsBunch->calcBeamParameters();
                        for(unsigned int i = 0; i < itsBunch->getLocalNum(); ++i) {
                            long &l = itsBunch->LastSection[i];
                            l = -1;
                            itsOpalBeamline_m.getSectionIndexAt (itsBunch->R[i], l);
                            itsBunch->ResetLocalCoordinateSystem (
                                i,
                                itsOpalBeamline_m.getOrientation(l),
                                itsOpalBeamline_m.getSectionStart(l));
                        }

                        /*
                          Check if there are any particles in simulation. If there are,
                          as in a restart, use the usual function to calculate beam
                          parameters. If not, calculate beam parameters of the initial
                          beam distribution.
                        */
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
                                        d->GetTypeofDistribution() +
                                        "\". Need to check the input file... ...");
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
                    /*
                      set thermal velocity of Maxwellian distribution of
                      secondaries for benchmark
                    */
                    double vVThermal = dist->GetvVThermal();
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
                // sey_0 in Vaughan's model
                double vSeyZero = dist->GetvSeyZero();
                // energy related to sey_0 in Vaughan's model
                double vEZero = dist->GetvEZero();
                // sey max in Vaughan's model
                double vSeyMax = dist->GetvSeyMax();
                // Emax in Vaughan's model
                double vEmax = dist->GetvEmax();
                // fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
                double vKenergy = dist->GetvKenergy();
                // fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
                double vKtheta = dist->GetvKtheta();
                // return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
                double vVThermal = dist->GetvVThermal();
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
            /*
              :FIXME: maybe need to be called in each time step for modeling
              creating darkcurrent in each time step
            */
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
            // set upper limit of particle number in simulation
            maxNparts_m = static_cast<size_t>(maxnp->getReal());
        }

        msg << "Boundary geometry initialized " << endl;
    }
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
    Inform msg("ParallelTTracker ", *gmsg);
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

void ParallelTTracker::writePhaseSpace(const long long step, const double &sposRef, bool psDump, bool statDump) {
    extern Inform *gmsg;
    Inform msg("OPAL ", *gmsg);
    Vector_t externalE, externalB;
    Vector_t FDext[6];  // FDext = {BHead, EHead, BRef, ERef, BTail, ETail}.

    // Sample fields at (xmin, ymin, zmin), (xmax, ymax, zmax) and the centroid location. We
    // are sampling the electric and magnetic fields at the back, front and
    // center of the beam.
    Vector_t rmin, rmax;
    itsBunch->get_bounds(rmin, rmax);

    Vector_t pos[3];
    pos[0] = Vector_t(rmax(0), rmax(1), rmax(2));
    pos[1] = Vector_t(0.0, 0.0, sposRef);
    pos[2] = Vector_t(rmin(0), rmin(1), rmin(2));

    if (psDump || statDump) {
        for(int k = 0; k < 3; ++k) {
            externalB = Vector_t(0.0);
            externalE = Vector_t(0.0);
            itsOpalBeamline_m.getFieldAt(pos[k], itsBunch->get_rmean(), itsBunch->getT() - 0.5 * itsBunch->getdT(), externalE, externalB);
            FDext[2 * k]   = externalB;
            FDext[2 * k + 1] = externalE * 1e-6;
        }
    }

    if (itsBunch->getTotalNum() > 0) {
        if(psDump) {
            // Write fields to .h5 file.
            itsDataSink_m->writePhaseSpace(*itsBunch, FDext, rmax(2), sposRef, rmin(2));
            msg << level3 << "* Wrote beam phase space." << endl;
            msg << *itsBunch << endl;
        } else {
            itsBunch->calcBeamParameters();
        }
    }

    if(statDump) {
        std::vector<std::pair<std::string, unsigned int> > collimatorLosses;
        FieldList collimators = itsOpalBeamline_m.getElementByType(ElementBase::COLLIMATOR);
	if (collimators.size() != 0) {
	  for (FieldList::iterator it = collimators.begin(); it != collimators.end(); ++ it) {
  	    Collimator* coll = static_cast<Collimator*>(it->getElement().get());
            std::string name = coll->getName();
            unsigned int losses = coll->getLosses();
            collimatorLosses.push_back(std::make_pair(name, losses));
	  }
	  std::sort(collimatorLosses.begin(), collimatorLosses.end(),
		    [](const std::pair<std::string, unsigned int>& a, const std::pair<std::string, unsigned int>& b) ->bool {
                      return a.first < b.first;
		    });
	  std::vector<unsigned int> bareLosses(collimatorLosses.size(),0);
	  for (size_t i = 0; i < collimatorLosses.size(); ++ i){
            bareLosses[i] = collimatorLosses[i].second;
	  }

	  reduce(&bareLosses[0], &bareLosses[0] + bareLosses.size(), &bareLosses[0], OpAddAssign());

	  for (size_t i = 0; i < collimatorLosses.size(); ++ i){
            collimatorLosses[i].second = bareLosses[i];
	  }
	}
        // Write statistical data.
        msg << level3 << "* Wrote beam statistics." << endl;
        itsDataSink_m->writeStatData(*itsBunch, FDext, rmax(2), sposRef, rmin(2), collimatorLosses);

        writeTrackOrbitFile(sposRef);
    }
}

#ifdef OPAL_DKS

void ParallelTTracker::registerHostMemory() {

    if (Ippl::getNodes() == 1) {
        dksbase.registerHostMemory(&itsBunch->R[0], itsBunch->getLocalNum());
        dksbase.registerHostMemory(&itsBunch->P[0], itsBunch->getLocalNum());
        dksbase.registerHostMemory(&itsBunch->X[0], itsBunch->getLocalNum());
        dksbase.registerHostMemory(&itsBunch->dt[0], itsBunch->getLocalNum());
        dksbase.registerHostMemory(&itsBunch->LastSection[0], itsBunch->getLocalNum());
    }
    numDeviceElements = itsBunch->getLocalNum();

}

void ParallelTTracker::unregisterHostMemory() {

    if (Ippl::getNodes() == 1) {
        dksbase.unregisterHostMemory(&itsBunch->R[0]);
        dksbase.unregisterHostMemory(&itsBunch->P[0]);
        dksbase.unregisterHostMemory(&itsBunch->X[0]);
        dksbase.unregisterHostMemory(&itsBunch->dt[0]);
        dksbase.unregisterHostMemory(&itsBunch->LastSection[0]);
    }

}

void ParallelTTracker::allocateDeviceMemory() {

    //allocate memory on device
    r_ptr = dksbase.allocateMemory<Vector_t>(itsBunch->getLocalNum(), ierr);
    p_ptr = dksbase.allocateMemory<Vector_t>(itsBunch->getLocalNum(), ierr);
    x_ptr = dksbase.allocateMemory<Vector_t>(itsBunch->getLocalNum(), ierr);
    lastSec_ptr = dksbase.allocateMemory<long>(itsBunch->getLocalNum(), ierr);
    dt_ptr = dksbase.allocateMemory<double>(itsBunch->getLocalNum(), ierr);

    numDeviceElements = itsBunch->getLocalNum();

}

void ParallelTTracker::freeDeviceMemory() {

    //free device memory
    dksbase.freeMemory<Vector_t>(r_ptr, numDeviceElements);
    dksbase.freeMemory<Vector_t>(p_ptr, numDeviceElements);
    dksbase.freeMemory<Vector_t>(x_ptr, numDeviceElements);
    dksbase.freeMemory<long>(lastSec_ptr, numDeviceElements);
    dksbase.freeMemory<double>(dt_ptr, numDeviceElements);

}

#endif

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c++
// c-basic-offset: 4
// indent-tabs-mode:nil
// End: