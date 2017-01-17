// -----------------------------------------------------------------------
// /*$RCSfile*/: TrackRun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackRun
//   The class for the OPAL RUN command.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackRun.h"
#include "AbstractObjects/BeamSequence.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include "Algorithms/Tracker.h"
#include "Algorithms/ThinTracker.h"
#include "Algorithms/ThickTracker.h"

#include "Algorithms/bet/EnvelopeBunch.h"

#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelSliceTracker.h"
#include "Algorithms/ParallelCyclotronTracker.h"
#include "Algorithms/AutophaseTracker.h"

#include "Attributes/Attributes.h"
#include "Beamlines/TBeamline.h"

#include "BasicActions/Option.h"

#include "Elements/OpalBeamBeam3D.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include "Utilities/Round.h"
#include "Structure/Beam.h"
#include "Structure/FieldSolver.h"
#include "Structure/DataSink.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPT.h"
#include "Structure/H5PartWrapperForPC.h"
#include "Structure/H5PartWrapperForPS.h"
#include "Distribution/Distribution.h"
#include "Structure/BoundaryGeometry.h"

#ifdef HAVE_AMR_SOLVER
#define DIM 3
#include <ParallelDescriptor.H>
#endif

#include <fstream>
#include <iomanip>

extern Inform *gmsg;

using namespace Physics;

// ------------------------------------------------------------------------

namespace {

    // The attributes of class TrackRun.
    enum {
        METHOD,       // Tracking method to use.
        FNAME,        // The name of file to be written.
        TURNS,        // The number of turns to be tracked.
        MBMODE,       // The working way for multi-bunch mode for OPAL-cycl: "FORCE" or "AUTO"
        PARAMB,       // The control parameter for "AUTO" mode of multi-bunch
        BEAM,         // The beam to track
        FIELDSOLVER,  // The field solver attached
        BOUNDARYGEOMETRY, // The boundary geometry
        DISTRIBUTION, // The particle distribution
        MULTIPACTING, // MULTIPACTING flag
        // THE INTEGRATION TIMESTEP IN SEC
        SIZE
    };
}

const std::string TrackRun::defaultDistribution("DISTRIBUTION");

TrackRun::TrackRun():
    Action(SIZE, "RUN",
           "The \"RUN\" sub-command tracks the defined particles through "
           "the given lattice."),
    itsTracker(NULL),
    dist(NULL),
    distrs_m(),
    fs(NULL),
    ds(NULL),
    phaseSpaceSink_m(NULL) {
    itsAttr[METHOD] = Attributes::makeString
                      ("METHOD", "Name of tracking algorithm to use:\n"
                       "\t\t\t\"THIN\" (default) or \"THICK,PARALLEL-T,CYCLOTRON-T,PARALLEL-SLICE,AUTOPHASE\".", "THIN");
    itsAttr[TURNS] = Attributes::makeReal
                     ("TURNS", "Number of turns to be tracked; Number of neighboring bunches to be tracked in cyclotron", 1.0);

    itsAttr[MBMODE] = Attributes::makeString
                      ("MBMODE", "The working way for multi-bunch mode for OPAL-cycl: FORCE or AUTO ", "FORCE");

    itsAttr[PARAMB] = Attributes::makeReal
                      ("PARAMB", " Control parameter to define when to start multi-bunch mode, only available in \"AUTO\" mode ", 5.0);

    itsAttr[FNAME] = Attributes::makeString
                     ("FILE", "Name of file to be written", "TRACK");

    itsAttr[BEAM] = Attributes::makeString
                    ("BEAM", "Name of beam ", "BEAM");
    itsAttr[FIELDSOLVER] = Attributes::makeString
                           ("FIELDSOLVER", "Field solver to be used ", "FIELDSOLVER");
    itsAttr[BOUNDARYGEOMETRY] = Attributes::makeString
                           ("BOUNDARYGEOMETRY", "Boundary geometry to be used NONE (default)", "NONE");
    itsAttr[DISTRIBUTION] = Attributes::makeStringArray
                             ("DISTRIBUTION", "List of particle distributions to be used ");
    itsAttr[MULTIPACTING] = Attributes::makeBool
                            ("MULTIPACTING", "Multipacting flag, default: false. Set true to initialize primary particles according to BoundaryGeometry", false);
    opal = OpalData::getInstance();
}


TrackRun::TrackRun(const std::string &name, TrackRun *parent):
    Action(name, parent),
    itsTracker(NULL),
    dist(NULL),
    distrs_m(),
    fs(NULL),
    ds(NULL),
    phaseSpaceSink_m(NULL) {
    opal = OpalData::getInstance();
}


TrackRun::~TrackRun()
{
    delete phaseSpaceSink_m;
    phaseSpaceSink_m = NULL;
}


TrackRun *TrackRun::clone(const std::string &name) {
    return new TrackRun(name, this);
}


void TrackRun::execute() {

    // Get algorithm to use.
    std::string method = Attributes::getString(itsAttr[METHOD]);
    if(method == "THIN") {
        //std::cerr << "  method == \"THIN\"" << std::endl;
        itsTracker = new ThinTracker(*Track::block->use->fetchLine(),
                                     *Track::block->bunch, Track::block->reference,
                                     false, false);
    } else if(method == "THICK") {
        //std::cerr << "  method == \"THICK\"" << std::endl;
        itsTracker = new ThickTracker(*Track::block->use->fetchLine(),
                                      *Track::block->bunch, Track::block->reference,
                                      false, false);
    } else if(method == "PARALLEL-SLICE") {
        setupSliceTracker();
    } else if(method == "PARALLEL-T") {
        setupTTracker();
    } else if(method == "PARALLEL-Z") {
        *gmsg << "  method == \"PARALLEL-Z\"" << endl;

    } else if(method == "CYCLOTRON-T") {
        setupCyclotronTracker();
    } else if(method == "AUTOPHASE") {
        executeAutophaseTracker();

        return;
    } else {
        throw OpalException("TrackRun::execute()",
                            "Method name \"" + method + "\" unknown.");
    }

    if(method == "THIN" || method == "THICK") {
        /*
          OLD SERIAL STUFF
        */
        // Open output file.
        std::string file = Attributes::getString(itsAttr[FNAME]);
        std::ofstream os(file.c_str());
        if(os.bad()) {
            throw OpalException("TrackRun::execute()",
                                "Unable to open output file \"" + file + "\".");
        }

        // Print initial conditions.
        os << "\nInitial particle positions:\n"
           << itsTracker->getBunch() << std::endl;

        int turns = int(Round(Attributes::getReal(itsAttr[TURNS])));
        // Track for the all but last turn.
        for(int turn = 1; turn < turns; ++turn) {
            itsTracker->execute();
            os << "\nParticle positions after turn " << turn << ":\n"
               << itsTracker->getBunch() << std::endl;
        }
        // Track last turn, with statistics.
        itsTracker->execute();

        // Print final conditions.
        os << "Particle positions after turn " << turns << ":\n"
           << itsTracker->getBunch() << std::endl;
        //    Track::block->bunch = itsTracker->getBunch();
    } else {
        itsTracker->execute();
        opal->setRestartRun(false);
    }

    opal->bunchIsAllocated();
    if(method == "PARALLEL-SLICE")
        opal->slbunchIsAllocated();

    delete itsTracker;
}

void TrackRun::setupSliceTracker() {
    OpalData::getInstance()->setInOPALEnvMode();
    bool isFollowupTrack = opal->hasSLBunchAllocated() && !Options::scan;
    if(!opal->hasSLBunchAllocated()) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-SLICE, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else if(isFollowupTrack) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-SLICE, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else if(opal->hasSLBunchAllocated() && Options::scan) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == PARALLEL-SLICE, SCAN TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam   *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack) {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPS(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    const size_t numberOfDistributions = distr_str.size();
    if (numberOfDistributions == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distr_str.at(0));
        dist->setNumberOfDistributions(numberOfDistributions);

        if(numberOfDistributions > 1) {
            *gmsg << "Found more than one distribution: ";
            for(size_t i = 1; i < numberOfDistributions; ++ i) {
                Distribution *d = Distribution::find(distr_str.at(i));

                d->setNumberOfDistributions(numberOfDistributions);
                distrs_m.push_back(d);

                *gmsg << " " << distr_str.at(i);
            }
            *gmsg << endl;
        }
    }

    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
    fs->initCartesianFields();

    double charge = 0.0;

    if(!opal->hasSLBunchAllocated()) {
        if(!opal->inRestartRun()) {

            dist->CreateOpalE(beam, distrs_m, Track::block->slbunch, 0.0, 0.0);
            opal->setGlobalPhaseShift(0.5 * dist->GetTEmission());

        } else {
            /***
                reload slice distribution
            */

            dist->DoRestartOpalE(*Track::block->slbunch, beam->getNumberOfParticles(), opal->getRestartStep(), phaseSpaceSink_m);
        }
    } else {
        charge = 1.0;
    }

    Track::block->slbunch->setdT(Track::block->dT.front());
    // set the total charge
    charge = beam->getCharge() * beam->getCurrent() / beam->getFrequency();
    Track::block->slbunch->setCharge(charge);
    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->slbunch->setCouplingConstant(coefE);
    //Track::block->slbunch->calcBeamParameters();


    if(!opal->inRestartRun()) {
        if(!opal->hasDataSinkAllocated()) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, -1));
    }

    ds = opal->getDataSink();

    if(!opal->hasBunchAllocated())
        *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *Track::block->slbunch  << endl;
    *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq
          << " Inputfile is " << opal->getInputFn() << endl;


    findPhasesForMaxEnergy();

    itsTracker = new ParallelSliceTracker(*Track::block->use->fetchLine(),
                                          dynamic_cast<EnvelopeBunch &>(*Track::block->slbunch),
                                          *ds,
                                          Track::block->reference,
                                          false, false,
                                          Track::block->localTimeSteps.front(),
                                          Track::block->zstop.front());
}

void TrackRun::setupTTracker(){
    OpalData::getInstance()->setInOPALTMode();
    bool isFollowupTrack = (opal->hasBunchAllocated() && !Options::scan);

    if(opal->hasBunchAllocated()) {
        if (Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "* Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK in SCAN MODE" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
            Track::block->bunch->setLocalTrackStep(0);
            Track::block->bunch->setGlobalTrackStep(0);
        } else {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "* Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
            Track::block->bunch->setLocalTrackStep(0);
        }
    } else {
        if (Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "* Selected Tracking Method == PARALLEL-T, NEW TRACK in SCAN MODE" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "* Selected Tracking Method == PARALLEL-T, NEW TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        }
    }

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            BoundaryGeometry *bg = BoundaryGeometry::find(Attributes::getString(itsAttr[BOUNDARYGEOMETRY]))->
                clone(getOpalName() + std::string("_geometry"));
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }

    setupFieldsolver();

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (isFollowupTrack) {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPT(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    double charge = setDistributionParallelT(beam);

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;

    if (!isFollowupTrack && !opal->inRestartRun())
        Track::block->bunch->setT(Track::block->t0_m);

    bool mpacflg = Attributes::getBool(itsAttr[MULTIPACTING]);
    if(!mpacflg) {
        Track::block->bunch->setCharge(charge);
        // set coupling constant
        double coefE = 1.0 / (4 * pi * epsilon_0);
        Track::block->bunch->setCouplingConstant(coefE);


        // statistical data are calculated (rms, eps etc.)
        Track::block->bunch->calcBeamParameters();
    } else {
        Track::block->bunch->setChargeZeroPart(charge);// just set bunch->qi_m=charge, don't set bunch->Q[] as we have not initialized any particle yet.
        Track::block->bunch->calcBeamParametersInitial();// we have not initialized any particle yet.
    }

#ifdef HAVE_AMR_SOLVER
    setupAMRSolver();
#endif

    if(!opal->inRestartRun()) {
        if(!opal->hasDataSinkAllocated() && !Options::scan) {
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else if(Options::scan) {
            ds = opal->getDataSink();
            if(ds)
                delete ds;
            opal->setDataSink(new DataSink(phaseSpaceSink_m));
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    } else {
        opal->setDataSink(new DataSink(phaseSpaceSink_m, -1));//opal->getRestartStep()));
    }

    ds = opal->getDataSink();

    if(opal->hasBunchAllocated() && Options::scan)
        ds->reset();

    if(!opal->hasBunchAllocated() || Options::scan) {
        if(!mpacflg) {
            *gmsg << *dist << endl;
        } else {
            *gmsg << "* Multipacting flag is true. The particle distribution in the run command will be ignored " << endl;
        }
    }

    if (Track::block->bunch->getTotalNum() > 0) {
        double spos = Track::block->bunch->get_sPos();
        auto &zstop = Track::block->zstop;
        auto &timeStep = Track::block->localTimeSteps;
        auto &dT = Track::block->dT;

        unsigned int i = 0;
        while (i + 1 < zstop.size() && zstop[i + 1] < spos) {
            ++ i;
        }

        zstop.erase(zstop.begin(), zstop.begin() + i);
        timeStep.erase(timeStep.begin(), timeStep.begin() + i);
        dT.erase(dT.begin(), dT.begin() + i);

        Track::block->bunch->setdT(dT.front());
    }

    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    *gmsg << *Track::block->bunch  << endl;

    findPhasesForMaxEnergy();

    *gmsg << level2
          << "Phase space dump frequency " << Options::psDumpFreq << " and "
          << "statistics dump frequency " << Options::statDumpFreq << " w.r.t. the time step." << endl;
    
#ifdef P3M_TEST

    Track::block->bunch->runTests();

#else
    itsTracker = new ParallelTTracker(*Track::block->use->fetchLine(),
                                      dynamic_cast<PartBunch &>(*Track::block->bunch), *ds,
                                      Track::block->reference, false, false, Track::block->localTimeSteps,
                                      Track::block->zstop, Track::block->timeIntegrator, Track::block->dT,  
				      beam->getNumberOfParticles());
#endif
    itsTracker->setMpacflg(mpacflg); // set multipacting flag in ParallelTTracker
}

void TrackRun::setupCyclotronTracker(){
    OpalData::getInstance()->setInOPALCyclMode();
    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    if (Attributes::getString(itsAttr[BOUNDARYGEOMETRY]) != "NONE") {
        // Ask the dictionary if BoundaryGeometry is allocated.
        // If it is allocated use the allocated BoundaryGeometry
        if (!OpalData::getInstance()->hasGlobalGeometry()) {
            BoundaryGeometry *bg = BoundaryGeometry::find(Attributes::getString(itsAttr[BOUNDARYGEOMETRY]))->
                clone(getOpalName() + std::string("_geometry"));
            OpalData::getInstance()->setGlobalGeometry(bg);
        }
    }

    setupFieldsolver();

    Track::block->bunch->PType = ParticleType::REGULAR;

    std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
    if (distr_str.size() == 0) {
        dist = Distribution::find(defaultDistribution);
    } else {
        dist = Distribution::find(distr_str.at(0));
    }

    // set macromass and charge for simulation particles
    double macromass = 0.0;
    double macrocharge = 0.0;

    const int specifiedNumBunch = int(std::abs(Round(Attributes::getReal(itsAttr[TURNS]))));

    if(opal->inRestartRun()) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  opal->getRestartStep(),
                                                  OpalData::getInstance()->getRestartFileName(),
                                                  H5_O_WRONLY);
    } else if (opal->hasBunchAllocated() && !Options::scan) {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  -1,
                                                  opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    } else {
        phaseSpaceSink_m = new H5PartWrapperForPC(opal->getInputBasename() + std::string(".h5"),
                                                  H5_O_WRONLY);
    }

    if(beam->getNumberOfParticles() < 3 || beam->getCurrent() == 0.0) {
        macrocharge = beam->getCharge() * q_e;
        macromass = beam->getMass();
        dist->CreateOpalCycl(*Track::block->bunch,
                             beam->getNumberOfParticles(),
                             beam->getCurrent(),*Track::block->use->fetchLine(),
                             Options::scan);

    } else {

        /**
           getFrequency() gets RF frequency [MHz], NOT isochronous  revolution frequency of particle!
           getCurrent() gets beamcurrent [A]

        */
        macrocharge = beam->getCurrent() / (beam->getFrequency() * 1.0e6); // [MHz]-->[Hz]

        if(!opal->hasBunchAllocated()) {
            if(!opal->inRestartRun()) {
                macrocharge /= beam->getNumberOfParticles();
                macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
                dist->CreateOpalCycl(*Track::block->bunch,
                                     beam->getNumberOfParticles(),
                                     beam->getCurrent(),
                                     *Track::block->use->fetchLine(),
                                     Options::scan);

            } else {
                dist->DoRestartOpalCycl(*Track::block->bunch,
                                        beam->getNumberOfParticles(),
                                        opal->getRestartStep(),
                                        specifiedNumBunch,
                                        phaseSpaceSink_m);
                macrocharge /= beam->getNumberOfParticles();
                macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
            }
        } else if(opal->hasBunchAllocated() && Options::scan) {
            macrocharge /= beam->getNumberOfParticles();
            macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
            dist->CreateOpalCycl(*Track::block->bunch,
                                 beam->getNumberOfParticles(),
                                 beam->getCurrent(),
                                 *Track::block->use->fetchLine(),
                                 Options::scan);
        }
    }
    Track::block->bunch->setMass(macromass); // set the Mass per macro-particle, [GeV/c^2]
    Track::block->bunch->setCharge(macrocharge);  // set the charge per macro-particle, [C]

    *gmsg << "* Mass of simulation particle= " << macromass << " GeV/c^2" << endl;
    *gmsg << "* Charge of simulation particle= " << macrocharge << " [C]" << endl;

    Track::block->bunch->setdT(1.0 / (Track::block->stepsPerTurn * beam->getFrequency() * 1.0e6));
    Track::block->bunch->setStepsPerTurn(Track::block->stepsPerTurn);

    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    if(!opal->inRestartRun())
        if(!opal->hasDataSinkAllocated()) {
            ds = new DataSink(phaseSpaceSink_m);
            opal->setDataSink(ds);
        } else {
            ds = opal->getDataSink();
            ds->changeH5Wrapper(phaseSpaceSink_m);
        }
    else {
        ds = new DataSink(phaseSpaceSink_m, -1);
        opal->setDataSink(ds);
    }

    if(opal->hasBunchAllocated() && Options::scan)
        ds->reset();

    if(!opal->hasBunchAllocated() && !Options::scan) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == CYCLOTRON-T, NEW TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else if(opal->hasBunchAllocated() && !Options::scan) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == CYCLOTRON-T, FOLLOWUP TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    } else if(opal->hasBunchAllocated() && Options::scan) {
        *gmsg << "* ********************************************************************************** " << endl;
        *gmsg << "* Selected Tracking Method == CYCLOTRON-T, SCAN TRACK" << endl;
        *gmsg << "* ********************************************************************************** " << endl;
    }
    *gmsg << "* Number of neighbour bunches= " << specifiedNumBunch << endl;
    *gmsg << "* DT                         = " << Track::block->dT.front() << endl;
    *gmsg << "* MAXSTEPS                   = " << Track::block->localTimeSteps.front() << endl;
    *gmsg << "* Phase space dump frequency = " << Options::psDumpFreq << endl;
    *gmsg << "* Statistics dump frequency  = " << Options::statDumpFreq << " w.r.t. the time step." << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    itsTracker = new ParallelCyclotronTracker(*Track::block->use->fetchLine(),
                                              dynamic_cast<PartBunch &>(*Track::block->bunch), *ds, Track::block->reference,
                                              false, false, Track::block->localTimeSteps.front(), Track::block->timeIntegrator);

    itsTracker->setNumBunch(specifiedNumBunch);

    if(opal->inRestartRun()) {

        H5PartWrapperForPC *h5pw = static_cast<H5PartWrapperForPC*>(phaseSpaceSink_m);
        itsTracker->setBeGa(h5pw->getMeanMomentum());

        itsTracker->setPr(h5pw->getReferencePr());
        itsTracker->setPt(h5pw->getReferencePt());
        itsTracker->setPz(h5pw->getReferencePz());

        itsTracker->setR(h5pw->getReferenceR());
        itsTracker->setTheta(h5pw->getReferenceT());
        itsTracker->setZ(h5pw->getReferenceZ());

        // The following is for restarts in local frame
        itsTracker->setPhi(h5pw->getAzimuth());
        itsTracker->setPsi(h5pw->getElevation());
        itsTracker->setPreviousH5Local(h5pw->getPreviousH5Local());
    }

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    // *gmsg << *Track::block->bunch  << endl;

    if(specifiedNumBunch > 1) {

        // only for regular  run of multi bunches, instantiate the  PartBins class
        // note that for restart run of multi bunches, PartBins class is instantiated in function doRestart_cycl()
        if(!opal->inRestartRun()) {

            // already exist bins number initially
            const int BinCount = 1;
            //allowed maximal bin
            const int MaxBinNum = 1000;

            // initialize particles number for each bin (both existed and not yet emmitted)
            size_t partInBin[MaxBinNum];
            for(int ii = 0; ii < MaxBinNum; ii++) partInBin[ii] = 0;
            partInBin[0] =  beam->getNumberOfParticles();

            Track::block->bunch->setPBins(new PartBinsCyc(MaxBinNum, BinCount, partInBin));
            // the allowed maximal bin number is set to 100
            //Track::block->bunch->setPBins(new PartBins(100));

        }

        // mode of generating new bunches:
        // "FORCE" means generating one bunch after each revolution, until get "TURNS" bunches.
        // "AUTO" means only when the distance between two neighbor bunches is bellow the limitation,
        //        then starts to generate new bunches after each revolution,until get "TURNS" bunches;
        //        otherwise, run single bunch track

        *gmsg << "***---------------------------- MULTI-BUNCHES MULTI-ENERGY-BINS MODE------ ----------------------------*** " << endl;

        double paraMb = Attributes::getReal(itsAttr[PARAMB]);
        itsTracker->setParaAutoMode(paraMb);

        if(opal->inRestartRun()) {

            itsTracker->setLastDumpedStep(opal->getRestartStep());

            if(Track::block->bunch->pbin_m->getLastemittedBin() < 2) {
                itsTracker->setMultiBunchMode(2);
                *gmsg << "In this restart job, the multi-bunches mode is forcely set to AUTO mode." << endl;
            } else {
                itsTracker->setMultiBunchMode(1);
                *gmsg << "In this restart job, the multi-bunches mode is forcely set to FORCE mode." << endl
                      << "If the existing bunch number is less than the specified number of TURN, readin the phase space of STEP#0 from h5 file consecutively" << endl;
            }
        } else {
            //////
            if((Attributes::getString(itsAttr[MBMODE])) == std::string("FORCE")) {
                itsTracker->setMultiBunchMode(1);
                *gmsg << "FORCE mode: The multi bunches will be injected consecutively after each revolution, until get \"TURNS\" bunches." << endl;


            }
            //////
            else if((Attributes::getString(itsAttr[MBMODE])) == std::string("AUTO")) {


                itsTracker->setMultiBunchMode(2);

                *gmsg << "AUTO mode: The multi bunches will be injected only when the distance between two neighborring bunches " << endl
                      << "is bellow the limitation. The control parameter is set to " << paraMb << endl;
            }
            //////
            else
                throw OpalException("TrackRun::execute()",
                                    "MBMODE name \"" + Attributes::getString(itsAttr[MBMODE]) + "\" unknown.");
        }

    }
}

void TrackRun::setupFieldsolver() {
    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
    fs->initCartesianFields();
    Track::block->bunch->setSolver(fs);
    if (fs->hasPeriodicZ())
        Track::block->bunch->setBCForDCBeam();
    else
        Track::block->bunch->setBCAllOpen();
}

double TrackRun::setDistributionParallelT(Beam *beam) {

    // If multipacting flag is not set, get distribution(s).
    if (!Attributes::getBool(itsAttr[MULTIPACTING])) {
        /*
         * Distribution(s) can be set via a single distribution or a list
         * (array) of distributions. If an array is defined the first in the
         * list is treated as the primary distribution. All others are added to
         * it to create the full distribution.
         */
        std::vector<std::string> distributionArray
            = Attributes::getStringArray(itsAttr[DISTRIBUTION]);
        const size_t numberOfDistributions = distributionArray.size();

        if (numberOfDistributions == 0) {
            dist = Distribution::find(defaultDistribution);
        } else {
            dist = Distribution::find(distributionArray.at(0));
            dist->setNumberOfDistributions(numberOfDistributions);

            if (numberOfDistributions > 1) {
                *gmsg << endl
                      << "---------------------------------" << endl
                      << "Found more than one distribution:" << endl << endl;
                *gmsg << "Main Distribution" << endl
                      << "---------------------------------" << endl
                      << distributionArray.at(0) << endl << endl
                      << "Secondary Distribution(s)" << endl
                      << "---------------------------------" << endl;

                for (size_t i = 1; i < numberOfDistributions; ++ i) {
                    Distribution *distribution = Distribution::find(distributionArray.at(i));
                    distribution->setNumberOfDistributions(numberOfDistributions);
                    distrs_m.push_back(distribution);

                    *gmsg << distributionArray.at(i) << endl;
                }
                *gmsg << endl
                      << "---------------------------------" << endl << endl;
            }
        }
    }

    /*
     * Initialize distributions.
     */
    size_t numberOfParticles = beam->getNumberOfParticles();
    if (!opal->hasBunchAllocated()) {
        if (!opal->inRestartRun()) {
            if (!Attributes::getBool(itsAttr[MULTIPACTING])) {
                /*
                 * Here we are not doing a restart or doing a mulitpactor run
                 * and we do not have a bunch already allocated.
                 */
                Track::block->bunch->setDistribution(dist,
                                                     distrs_m,
                                                     numberOfParticles,
                                                     Options::scan);

                /*
                 * If this is an injected beam (rather than an emitted beam), we
                 * make sure it doesn't have any particles at z < 0.
                 */
                Vector_t rMin;
                Vector_t rMax;
                Track::block->bunch->get_bounds(rMin, rMax);

                opal->setGlobalPhaseShift(0.5 * dist->GetTEmission());
            }
        } else {
            /*
             * Read in beam from restart file.
             */
            dist->DoRestartOpalT(*Track::block->bunch, numberOfParticles, opal->getRestartStep(), phaseSpaceSink_m);
        }
    } else if (opal->hasBunchAllocated() && Options::scan) {
        /*
         * We are in scan mode and have already allocated a bunch. So, we need to
         * do a couple more things.
         */
        Track::block->bunch->setDistribution(dist,
                                             distrs_m,
                                             numberOfParticles,
                                             Options::scan);
        Track::block->bunch->resetIfScan();
        Track::block->bunch->LastSection = 1;

        opal->setGlobalPhaseShift(0.5 * dist->GetTEmission());
    }

    // Return charge per macroparticle.
    return beam->getCharge() * beam->getCurrent() / beam->getFrequency() / numberOfParticles;

}

void TrackRun::findPhasesForMaxEnergy(bool writeToFile) const {
    if (Options::autoPhase == 0 ||
        OpalData::getInstance()->inRestartRun() ||
        OpalData::getInstance()->hasBunchAllocated()) return;

    std::queue<unsigned long long> localTrackSteps;
    std::queue<double> dtAllTracks;
    std::queue<double> zStop;

    const PartBunch *bunch = Track::block->bunch;
    if (bunch == NULL) {
        bunch = Track::block->slbunch;
    }
    if (bunch == NULL) {
        throw OpalException("findPhasesForMaxEnergy()",
                            "no valid PartBunch object found");
    }

    Vector_t meanP = bunch->get_pmean();
    Vector_t meanR = bunch->get_rmean();
    if (bunch->getTotalNum() == 0) {
        meanP = bunch->get_pmean_Distribution();
        meanR = 0.0;
    }
    AutophaseTracker apTracker(*Track::block->use->fetchLine(),
                               Track::block->reference,
                               bunch->getT(),
                               meanR(2),
                               meanP(2));

    {
        std::vector<unsigned long long>::const_iterator it = Track::block->localTimeSteps.begin();
        std::vector<unsigned long long>::const_iterator end = Track::block->localTimeSteps.end();
        for (; it != end; ++ it) {
            localTrackSteps.push(*it);
        }
    }

    {
        std::vector<double>::const_iterator it = Track::block->dT.begin();
        std::vector<double>::const_iterator end = Track::block->dT.end();
        for (; it != end; ++ it) {
            dtAllTracks.push(*it);
        }
    }

    {
        std::vector<double>::const_iterator it = Track::block->zstop.begin();
        std::vector<double>::const_iterator end = Track::block->zstop.end();
        for (; it != end; ++ it) {
            zStop.push(*it);
        }
    }

    apTracker.execute(dtAllTracks, zStop, localTrackSteps);

    if (writeToFile) {
        std::string fileName = OpalData::getInstance()->getInputBasename() + "_phases.in";
        apTracker.save(fileName);
    }
}

void TrackRun::executeAutophaseTracker() {

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
    /*double charge =*/ setDistributionParallelT(beam);

    Track::block->bunch->setdT(Track::block->dT.front());
    Track::block->bunch->dtScInit_m = Track::block->dtScInit;
    Track::block->bunch->deltaTau_m = Track::block->deltaTau;
    Track::block->bunch->setT(Track::block->t0_m);

    // Track::block->bunch->setCharge(charge);

    double couplingConstant = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(couplingConstant);


    Track::block->bunch->calcBeamParameters();

    findPhasesForMaxEnergy(true);

}

#ifdef HAVE_AMR_SOLVER
void TrackRun::setupAMRSolver() {
    /*
    if ( fs->isAMRSolver() ) {
        *gmsg << *Track::block->bunch  << endl;
        *gmsg << *fs   << endl;
        
        std::vector<double> x(3);
        std::vector<double> attr(11);

        for (size_t i=0; i<Track::block->bunch->getLocalNum(); i++) {
            // X, Y, Z are stored separately from the other attributes
            for (unsigned int k=0; k<3; k++)
                x[k] = Track::block->bunch->R[i](k);

            //  This allocates 11 spots -- 1 for Q, 3 for v, 3 for E, 3 for B, 1 for ID.
            //  IMPORTANT: THIS ORDERING IS ASSUMED WHEN WE FILL E AT THE PARTICLE LOCATION
            //             IN THE MOVEKICK CALL -- if Evec no longer starts at component 4
            //             then you must change "start_comp_for_e" in Accel_advance.cpp
	    //
            // Q      : 0
            //  Vvec   : 1, 2, 3 the velocity
            //  Evec   : 4, 5, 6 the electric field at the particle location
            //  Bvec   : 7, 8, 9 the electric field at the particle location
            //  id+1   : 10 (we add 1 to make the particle ID > 0)

            // This is the charge
            attr[0] = Track::block->bunch->Q[i];

            // These are the velocity components
            double gamma=sqrt(1+ dot(Track::block->bunch->P[i],Track::block->bunch->P[i]));
            for (unsigned int k=0; k<DIM; k++)
                attr[k+1] = Track::block->bunch->P[i](k) * Physics::c /gamma;

            // These are E and B
            for (unsigned int k=4; k<10; k++)
                attr[k]= 0.0;

            //
            // The Particle stuff in AMR requires ids > 0
            //   (because we flip the sign to make them invalid)
            // So we just make id->id+1 here.
            int particle_id = Track::block->bunch->ID[i] + 1;
            attr[3*DIM + 1] = particle_id;

//             fs->getAmrPtr()->addOneParticle(particle_id, Ippl::myNode(), x, attr);
        }

        // It is essential that we call this routine since the particles
        //    may not currently be defined on the same processor as the grid
        //    that will hold them in the AMR world.
        fs->getAmrPtr()->RedistributeParticles();
    
        // This part of the call must come after we add the particles
        // since this one calls post_init which does the field solve.
        double start_time = 0.0;
        double stop_time = 1.0;
        
        fs->getAmrPtr()->FinalizeInit(start_time, stop_time);
        fs->getAmrPtr()->writePlotFile();
    
        *gmsg << "A M R Initialization DONE" << endl;
    }
    */
}

#endif
