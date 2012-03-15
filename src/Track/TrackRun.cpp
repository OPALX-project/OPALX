// ------------------------------------------------------------------------
// $RCSfile: TrackRun.cpp,v $
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
#include "Algorithms/ThinTracker.h"
#include "Algorithms/ThickTracker.h"
#include "Algorithms/ParallelTTracker.h"
#ifdef HAVE_ENVELOPE_SOLVER
#include "Algorithms/ParallelSliceTracker.h"
#endif
#include "Algorithms/ParallelCyclotronTracker.h"

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
#include "Distribution/Distribution.h"

#include <fstream>
#include <iomanip>


using namespace Options;
using namespace Physics;


// Class TrackRun
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
    DISTRIBUTION, // The particle distribution
    // THE INTEGRATION TIMESTEP IN SEC
    SIZE
  };
}

TrackRun::TrackRun():
  Action(SIZE, "RUN",
	 "The \"RUN\" sub-command tracks the defined particles through "
	 "the given lattice.")
{
  itsAttr[METHOD] = Attributes::makeString
    ("METHOD", "Name of tracking algorithm to use:\n"
               "\t\t\t\"THIN\" (default) or \"THICK,PARALLEL-T,PARALLEL-Z,PARALLEL-SLICE\".", "THIN");
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
  itsAttr[DISTRIBUTION] = Attributes::makeString
    ("DISTRIBUTION", "Particle distribution to be used ", "DISTRIBUTION");
}


TrackRun::TrackRun(const string &name, TrackRun *parent):
  Action(name, parent)
{}


TrackRun::~TrackRun()
{}


TrackRun *TrackRun::clone(const string &name)
{
  return new TrackRun(name, this);
}

void TrackRun::execute()
{

  // Get algorithm to use.
  string method = Attributes::getString(itsAttr[METHOD]);

  if (method == "THIN") {
    //std::cerr << "  method == \"THIN\"" << std::endl;
    itsTracker = new ThinTracker(*Track::block->use->fetchLine(),
                                 *Track::block->bunch, Track::block->reference,
                                 false, false);
  } else if (method == "THICK") {
    //std::cerr << "  method == \"THICK\"" << std::endl;
    itsTracker = new ThickTracker(*Track::block->use->fetchLine(),
                                 *Track::block->bunch, Track::block->reference,
                                 false, false);
  } else if (method == "PARALLEL-SLICE") {
#ifdef HAVE_ENVELOPE_SOLVER
    if (!OPAL.hasSLBunchAllocated()) {
      *gmsg << "* ********************************************************************************** " << endl;
      *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, NEW TRACK" << endl;
      *gmsg << "* ********************************************************************************** " << endl;
    }
    else {
      *gmsg << "* ********************************************************************************** " << endl;
      *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, FOLLOWUP TRACK" << endl;
      *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam   *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));

    double charge = 0.0;

    if (!OPAL.hasSLBunchAllocated()) {
      if (!OPAL.inRestartRun()) {

	/**
	   Create slice distribution
	*/
	dist->createSlicedBunch(*Track::block->slbunch);
      }
      else {

	/**
	   reload slice distribution
	*/
	// dist->doRestart(*Track::block->slbunch, beam->getNumberOfParticles(), OPAL.getRestartStep());
      }
    }
    else {
      charge = 1.0;
    }

    Track::block->slbunch->setdT(Track::block->dT);

    // set the total charge
    Track::block->slbunch->setCharge(charge);

    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->slbunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->slbunch->calcBeamParameters();

    if (!OPAL.inRestartRun()) {
      if (!OPAL.hasSLDataSinkAllocated())
	OPAL.setSLDataSink(new SLDataSink());
    }
    else
      OPAL.setSLDataSink(new SLDataSink(OPAL.getRestartStep()+1));

    slds = OPAL.getSLDataSink();

    if (!OPAL.hasBunchAllocated())
      *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *Track::block->slbunch  << endl;

    *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq
	  << " Inputfile is " << OPAL.getInputFn() << endl;

    // write into text file
    string fn = OPAL.getInputFn();
    int pos=fn.find(string("."),0);
    fn.erase(pos,fn.size()-pos);
    fn += string(".stat");

    ds->writeStatData(*Track::block->bunch, fn, string("IC"));

    itsTracker = new ParallelSliceTracker(*Track::block->use->fetchLine(),
					  dynamic_cast<SLPartBunch&>(*Track::block->slbunch), *slds, Track::block->reference,
					  false, false,Track::block->maxTSteps);
#else
      *gmsg << "* ********************************************************************************** " << endl;
      *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, but the code is not compiled for the " << endl;
      *gmsg << "  Use --envelope-solver in the configure make process   " << endl;
      *gmsg << "* ********************************************************************************** " << endl;
      exit(1);
#endif

  } else if (method == "PARALLEL-T") {

    if (!OPAL.hasBunchAllocated()) {
    *gmsg << "* ********************************************************************************** " << endl;
    *gmsg << "  Selected Tracking Method == PARALLEL-T, NEW TRACK" << endl;
    *gmsg << "* ********************************************************************************** " << endl;
    }
    else {
    *gmsg << "* ********************************************************************************** " << endl;
    *gmsg << "  Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK" << endl;
    *gmsg << "* ********************************************************************************** " << endl;
    }

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));

    fs->initCartesianFields();

    Track::block->bunch->setSolver(fs);

    Track::block->bunch->setBCAllOpen();

    dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));

    double charge = beam->getCharge() * beam->getCurrent() / beam->getFrequency();

    if (!OPAL.hasBunchAllocated()) {
      if (!OPAL.inRestartRun()) {
    	  dist->create(*Track::block->bunch, beam->getNumberOfParticles());
    	  charge /= beam->getNumberOfParticles();
      }
      else {
    	  dist->doRestart(*Track::block->bunch, beam->getNumberOfParticles(), OPAL.getRestartStep());
    	  charge /= Track::block->bunch->getTotalNum();
      }
    }
    else {
      charge /= Track::block->bunch->getTotalNum();
    }
    // fixMe: put boundp in distribution
    // Track::block->bunch.boundp();
    Track::block->bunch->setdT(Track::block->dT);

    // put bin zero out of the cathode
    // so that we have particles to track
    //if (Track::block->bunch.doEmission())
    //  Track::block->bunch.pbin_m->setBinEmitted(0);

    // set the total charge
    Track::block->bunch->setCharge(charge);

    // set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters();

    if (!OPAL.inRestartRun()) {
      if (!OPAL.hasDataSinkAllocated())
    	  OPAL.setDataSink(new DataSink());
    } else
      OPAL.setDataSink(new DataSink(OPAL.getRestartStep()+1));

    ds = OPAL.getDataSink();

    if (!OPAL.hasBunchAllocated())
      *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    *gmsg << *Track::block->bunch  << endl;

    *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq << endl
		<< " Inputfile is " << OPAL.getInputFn() << endl;

    // write into text file
    string fn = OPAL.getInputFn();
    int pos=fn.find(string("."),0);
    fn.erase(pos,fn.size()-pos);
    fn += string(".stat");

    // fort.xx
    ds->writeStatData(*Track::block->bunch, fn, string("IC"));

    // write into H5part container
    // ds->writePhaseSpace(Track::block->bunch, string("IC"), 0);

    itsTracker = new ParallelTTracker(*Track::block->use->fetchLine(),
									  dynamic_cast<PartBunch&>(*Track::block->bunch), *ds,
									  Track::block->reference, false, false, Track::block->maxTSteps);
  } else if (method == "PARALLEL-Z") {
    *gmsg << "  method == \"PARALLEL-Z\"" << endl;

  } else if (method == "CYCLOTRON-T") {

    *gmsg << "* ********************************************************************************** " << endl;
    *gmsg << "  Selected Tracking Method == CYCLOTRON-T" << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    fs   = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
    fs->initCartesianFields();
    Track::block->bunch->setSolver(fs);
    Track::block->bunch->setBCAllOpen();

    dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));

    // getCharge() gets number of charge, for proton it is 1.0
    // getFrequency() gets RF frequency [MHz], NOT isochronous  revolution frequency of particle!
    // getCurrent() gets beamcurrent [A]
    double charge = beam->getCurrent() / (beam->getFrequency()*1.0e6);// [MHz]-->[Hz]

    *gmsg <<"Beam current ="<<beam->getCurrent()*1000.0<<" [mA]"<<endl;
    *gmsg <<"Total charge in a bunch ="<<charge*1.0E12<<" [pC]"<<endl;

    const int specifiedNumBunch = int(abs(Round(Attributes::getReal(itsAttr[TURNS]))));
    *gmsg <<"Number of tracked neighbour bunches is set to "<<specifiedNumBunch<<endl;

    if (!OPAL.inRestartRun()) {

      dist->create(*Track::block->bunch, beam->getNumberOfParticles());
      charge /= beam->getNumberOfParticles();
    }
    else{

      dist->doRestart_cycl(*Track::block->bunch, beam->getNumberOfParticles(), OPAL.getRestartStep(), specifiedNumBunch);
      charge /= Track::block->bunch->getTotalNum();
    }

    Track::block->bunch->setdT(1.0/(Track::block->stepsPerTurn*beam->getFrequency()*1.0e6));
    Track::block->bunch->setStepsPerTurn(Track::block->stepsPerTurn);

    /// set the charge per macro-particle
    Track::block->bunch->setCharge(charge);

    /// set coupling constant
    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);

    // statistical data are calculated (rms, eps etc.)
    Track::block->bunch->calcBeamParameters_cycl();

    *gmsg << "* DT          " << Track::block->dT << endl;
    *gmsg << "* MAXSTEPS    " << Track::block->maxTSteps << endl;

    if (!OPAL.inRestartRun())
      if (!OPAL.hasDataSinkAllocated()) {
	ds = new DataSink();
	OPAL.setDataSink(ds);
      }
      else
	ds = OPAL.getDataSink();
    else {
      ds = new DataSink(OPAL.getRestartStep()+1);
      OPAL.setDataSink(ds);
    }

    *gmsg << *dist << endl;
    *gmsg << *beam << endl;
    *gmsg << *fs   << endl;
    *gmsg << *Track::block->bunch  << endl;

    *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq<<endl
          << "Inputfile is " << OPAL.getInputFn() << endl;

    // write into text file
    string fn = OPAL.getInputFn();
    int pos=fn.find(string("."),0);
    fn.erase(pos,fn.size()-pos);
    fn += string(".stat");

    ds->writeStatData(*Track::block->bunch, fn, string("IC"));

    itsTracker = new ParallelCyclotronTracker(*Track::block->use->fetchLine(),
                                              dynamic_cast<PartBunch&>(*Track::block->bunch), *ds, Track::block->reference,
                                              false, false,Track::block->maxTSteps);

    itsTracker->setNumBunch(specifiedNumBunch);

    string ParticleName = beam->getParticleName();
    *gmsg<<"The particle species: "<<ParticleName<<endl;

    const double ratioCh_M = beam->getCharge()/beam->getMass(); // chargeNumber/Mass[GeV]
    *gmsg<<"Charge number = "<<beam->getCharge()<<endl;
    *gmsg<<"Rest mass = "<<beam->getMass()<< " [GeV]"<<endl;
    if( ParticleName != string("PROTON"))
      itsTracker->setTrackCoeff(ratioCh_M);

    if ( specifiedNumBunch > 1 ){

      // only for regular  run of multi bunches, instantiate the  PartBins class
      // note that for restart run of multi bunches, PartBins class is instantiated in function doRestart_cycl()
      if (!OPAL.inRestartRun()) {

        // initialize particles number for each bin (both existed and not yet emmitted)
        size_t partInBin[specifiedNumBunch];
        for (int ii=0; ii<specifiedNumBunch; ii++) partInBin[ii] = 0;

        // already exist bins number
        int BinCount = 0;

        partInBin[0] =  beam->getNumberOfParticles();
        BinCount=1;
        Track::block->bunch->setPBins( new PartBins(specifiedNumBunch, BinCount, partInBin ));
      }

      // mode of generating new bunches:
      // "FORCE" means generating one bunch after each revolution, until get "TURNS" bunches.
      // "AUTO" means only when the distance between two neighbor bunches is bellow the limitation,
      //        then starts to generate new bunches after each revolution,until get "TURNS" bunches;
      //        otherwise, run single bunch track

      *gmsg << "***---------------------------- MULTI-BUNCHES MULTI-ENERGY-BINS MODE------ ----------------------------*** "<<endl;

      double paraMb = Attributes::getReal(itsAttr[PARAMB]);
      itsTracker->setParaAutoMode(paraMb);

      if (OPAL.inRestartRun()) {

        itsTracker->setLastDumpedStep(OPAL.getRestartStep());

        if ( Track::block->bunch->pbin_m->getLastemittedBin() < 2 ) {
          itsTracker->setMultiBunchMode(2);

          *gmsg<<"In this restart job, the multi-bunches mode is forcely set to AUTO mode."<<endl;
        }
        else {
          itsTracker->setMultiBunchMode(1);
          *gmsg<<"In this restart job, the multi-bunches mode is forcely set to FORCE mode."<<endl
               <<"If the existing bunch number is less than the specified number of TURN, readin the phase space of STEP#0 from h5 file consecutively"<<endl;
        }
      }
      else {
        //////
        if ( (Attributes::getString(itsAttr[MBMODE])) == string("FORCE") )
        {
          itsTracker->setMultiBunchMode(1);
          *gmsg<< "FORCE mode: The multi bunches will be injected consecutively after each revolution, until get \"TURNS\" bunches."<<endl;


        }
        //////
        else if ( (Attributes::getString(itsAttr[MBMODE])) == string("AUTO") )
        {


          itsTracker->setMultiBunchMode(2);

          *gmsg << "AUTO mode: The multi bunches will be injected only when the distance between two neighborring bunches " <<endl
                << "is bellow the limitation. The control parameter is set to "<< paraMb <<endl;
        }
        //////
        else
          throw OpalException("TrackRun::execute()",
                              "MBMODE name \"" +Attributes::getString(itsAttr[MBMODE])+ "\" unknown.");
      }

    }

  } else {
    throw OpalException("TrackRun::execute()",
                        "Method name \"" + method + "\" unknown.");
  }

  if ((method != "PARALLEL-Z") && (method != "PARALLEL-T") && (method != "CYCLOTRON-T")) {
    /*
      OLD SERIAL STUFF
    */
    // Open output file.
    string file = Attributes::getString(itsAttr[FNAME]);
    std::ofstream os(file.c_str());
    if (os.bad()) {
      throw OpalException("TrackRun::execute()",
			 "Unable to open output file \"" + file + "\".");
    }

    // Print initial conditions.
    os << "\nInitial particle positions:\n"
       << itsTracker->getBunch() << std::endl;

    int turns = int(Round(Attributes::getReal(itsAttr[TURNS])));
    // Track for the all but last turn.
    for (int turn = 1; turn < turns; ++turn) {
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
  }
  else {
   itsTracker->execute();
  }
  //  delete [] fs;
  //  delete [] dist;
  OPAL.bunchIsAllocated();
}
