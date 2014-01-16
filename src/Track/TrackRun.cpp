// -----------------------------------------------------------------------
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
#include "Algorithms/Tracker.h"
#include "Algorithms/ThinTracker.h"
#include "Algorithms/ThickTracker.h"

#include "Algorithms/bet/EnvelopeBunch.h"

#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelSliceTracker.h"
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

#ifdef HAVE_AMR_SOLVER
#define DIM 3
#include <Amr.H>
#include <ParallelDescriptor.H>
#include "Solvers/amr/BoundaryDomain.h"
#endif

#include <fstream>
#include <iomanip>

extern Inform *gmsg;

using namespace Options;
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
        DISTRIBUTION, // The particle distribution
        DISTRIBUTIONS, // A list of  particle distributions
        MULTIPACTING, // MULTIPACTING flag
        // THE INTEGRATION TIMESTEP IN SEC
        SIZE
    };
}

TrackRun::TrackRun():
    Action(SIZE, "RUN",
           "The \"RUN\" sub-command tracks the defined particles through "
           "the given lattice.") {
    itsAttr[METHOD] = Attributes::makeString
                      ("METHOD", "Name of tracking algorithm to use:\n"
                       "\t\t\t\"THIN\" (default) or \"THICK,PARALLEL-T,PARALLEL-TA,PARALLEL-Z,PARALLEL-SLICE\".", "THIN");
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
    itsAttr[DISTRIBUTIONS] = Attributes::makeStringArray
                             ("DISTRIBUTIONS", "List of particle distributions to be used ");
    itsAttr[MULTIPACTING] = Attributes::makeBool
                            ("MULTIPACTING", "Multipacting flag, default: false. Set true to initialize primary particles according to BoundaryGeometry", false);
    OPAL = OpalData::getInstance();
}


TrackRun::TrackRun(const string &name, TrackRun *parent):
    Action(name, parent) {
    OPAL = OpalData::getInstance();
}


TrackRun::~TrackRun()
{}


TrackRun *TrackRun::clone(const string &name) {
    return new TrackRun(name, this);
}


void TrackRun::execute() {

    // Get algorithm to use.
    std::string method = Attributes::getString(itsAttr[METHOD]);
    bool mpacflg = Attributes::getBool(itsAttr[MULTIPACTING]);
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
        OpalData::getInstance()->setInOPALEnvMode();
        if(!OPAL->hasSLBunchAllocated()) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, NEW TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else if(OPAL->hasSLBunchAllocated() && !Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, FOLLOWUP TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else if(OPAL->hasSLBunchAllocated() && Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-SLICE, SCAN TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        }

        Beam   *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
        dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));


        fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
        fs->initCartesianFields();

        double charge = 0.0;

        std::vector<std::string> distr_str = Attributes::getStringArray(itsAttr[DISTRIBUTIONS]);
        if(distr_str.size() > 0) {
            *gmsg << "Found more than one distribution: ";
            for(std::vector<std::string>::const_iterator dit = distr_str.begin(); dit != distr_str.end(); ++ dit) {
                Distribution *d = Distribution::find(*dit);
                *gmsg << " " << *dit;
                distrs_m.push_back(d);
            }
            *gmsg << endl;
        }

        if(!OPAL->hasSLBunchAllocated()) {
            if(!OPAL->inRestartRun()) {

                dist->CreateOpalE(beam, distrs_m, Track::block->slbunch, 0.0, 0.0);
                OPAL->setGlobalPhaseShift(0.5 * dist->GetTEmission());

            } else {
                /***
                    reload slice distribution
                */

                dist->DoRestartOpalE(*Track::block->slbunch, beam->getNumberOfParticles(), OPAL->getRestartStep());
            }
        } else {
            charge = 1.0;
        }

        Track::block->slbunch->setdT(Track::block->dT);
        // set the total charge
        charge = beam->getCharge() * beam->getCurrent() / beam->getFrequency();
        Track::block->slbunch->setCharge(charge);
        // set coupling constant
        double coefE = 1.0 / (4 * pi * epsilon_0);
        Track::block->slbunch->setCouplingConstant(coefE);
        //Track::block->slbunch->calcBeamParameters();


        if(!OPAL->inRestartRun()) {
            if(!OPAL->hasDataSinkAllocated())
                OPAL->setDataSink(new DataSink());
        } else
            OPAL->setDataSink(new DataSink(OPAL->getRestartStep() + 1));

        ds = OPAL->getDataSink();

        if(!OPAL->hasBunchAllocated())
            *gmsg << *dist << endl;
        *gmsg << *beam << endl;
        *gmsg << *Track::block->slbunch  << endl;
        *gmsg << "Phase space dump frequency is set to " << Options::psDumpFreq
              << " Inputfile is " << OPAL->getInputFn() << endl;


        ParallelTTracker *mySlApTracker = setupForAutophase();
        itsTracker = new ParallelSliceTracker(*Track::block->use->fetchLine(),
                                              dynamic_cast<EnvelopeBunch &>(*Track::block->slbunch),
                                              *ds,
                                              Track::block->reference,
                                              false, false,
                                              Track::block->localTimeSteps,
                                              Track::block->zstop,
                                              *mySlApTracker);
    } else if(method == "PARALLEL-T") {
        OpalData::getInstance()->setInOPALTMode();
        bool isFollowupTrack = (OPAL->hasBunchAllocated() && !Options::scan);

        if(!OPAL->hasBunchAllocated() && !Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-T, NEW TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else if(isFollowupTrack) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
            Track::block->bunch->setLocalTrackStep(0);
        } else if(OPAL->hasBunchAllocated() && Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-T, FOLLOWUP TRACK in SCAN MODE" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
            Track::block->bunch->setLocalTrackStep(0);
            Track::block->bunch->setGlobalTrackStep(0);
        } else if(!OPAL->hasBunchAllocated() && Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == PARALLEL-T, NEW TRACK in SCAN MODE" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else
            *gmsg << "  Selected Tracking Method is NOT implemented, good luck ..." << endl;

        Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));
        fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
        fs->initCartesianFields();
        Track::block->bunch->setSolver(fs);
	if (fs->hasPeriodicZ())
	  Track::block->bunch->setBCForDCBeam();
	else
	  Track::block->bunch->setBCAllOpen();


        double charge = SetDistributionParallelT(beam);

        Track::block->bunch->setdT(Track::block->dT);
        Track::block->bunch->dtScInit_m = Track::block->dtScInit;
        Track::block->bunch->deltaTau_m = Track::block->deltaTau;

        if (!isFollowupTrack && !OPAL->inRestartRun())
            Track::block->bunch->setT(Track::block->t0_m);

        if(!mpacflg) {
            Track::block->bunch->setCharge(charge);
        } else {
            Track::block->bunch->setChargeZeroPart(charge);// just set bunch->qi_m=charge, don't set bunch->Q[] as we have not initialized any particle yet.
        }
        if(!mpacflg) {
            // set coupling constant
            double coefE = 1.0 / (4 * pi * epsilon_0);
            Track::block->bunch->setCouplingConstant(coefE);


            // statistical data are calculated (rms, eps etc.)
            Track::block->bunch->calcBeamParameters();
        } else {
            Track::block->bunch->calcBeamParametersInitial();// we have not initialized any particle yet.
        }
#ifdef HAVE_AMR_SOLVER
	Amr* amrptr;
//	if (fs->isAMRSolver()) {
	  *gmsg << "A M R Initialization " << endl;
	  *gmsg << *Track::block->bunch  << endl;
	  *gmsg << *fs   << endl;

	  FieldLayout<3>::iterator_iv locDomBegin = Track::block->bunch->getFieldLayout().begin_iv();
	  FieldLayout<3>::iterator_iv locDomEnd = Track::block->bunch->getFieldLayout().end_iv();
	  FieldLayout<3>::iterator_dv globDomBegin = Track::block->bunch->getFieldLayout().begin_rdv();
	  FieldLayout<3>::iterator_dv globDomEnd = Track::block->bunch->getFieldLayout().end_rdv();
	  
	  NDIndex<3> ipplDom = Track::block->bunch->getFieldLayout().getDomain();
	
	  BoxArray lev0_grids(Ippl::getNodes());
	
	  Array<int> procMap;
	  procMap.resize(lev0_grids.size()+1); // +1 is a historical thing, do not ask	

	  // first iterate over the local owned domain(s)
	  for(FieldLayout<3>::const_iterator_iv v_i = locDomBegin ; v_i != locDomEnd; ++v_i) {
	    std::ostringstream stream;
	    stream << *((*v_i).second);	 

	    std::pair<Box,unsigned int> res = getBlGrids(stream.str());
	    lev0_grids.set(res.second,res.first);
	    procMap[res.second] = Ippl::myNode();
	  }

	  // then iterate over the non-local domain(s)
	  for(FieldLayout<3>::iterator_dv v_i = globDomBegin ; v_i != globDomEnd; ++v_i) {
	    std::ostringstream stream;
	    stream << *((*v_i).second);
	    
	    std::pair<Box,unsigned int> res = getBlGrids(stream.str());
	    lev0_grids.set(res.second,res.first);
	    procMap[res.second] = res.second;
	  }
	  procMap[lev0_grids.size()] = Ippl::myNode();
	  
	  // This init call will cache the distribution map as determined by procMap
	  // so that all data will end up on the right processor
	  RealBox rb;
	  Array<Real> prob_lo(3);
	  Array<Real> prob_hi(3);

	  Vector_t rmin;
	  Vector_t rmax;
	  Track::block->bunch->get_bounds(rmin, rmax);

          // Set up a problem size with dx = dy = dz and large enough
          //  to hold the geometric boundary.
	  prob_lo.set(0,-0.025);
	  prob_lo.set(1,-0.025);
	  prob_lo.set(2,-0.050);
	  prob_hi.set(0, 0.025);
	  prob_hi.set(1, 0.025);
	  prob_hi.set(2, 0.050);

	  rb.setLo(prob_lo);
	  rb.setHi(prob_hi);
	  
	  int coord_sys = 0;

	  Array<int> ncell(3);
	  ncell[0] = ipplDom[0].length();
	  ncell[1] = ipplDom[1].length();
	  ncell[2] = ipplDom[2].length();

          std::vector<int   > nr(3);
          std::vector<double> hr(3);
          std::vector<double> prob_lo_in(3);
          for (int i = 0; i < 3; i++)
          {
             nr[i] = ncell[i];
             hr[i] = (prob_hi[i] - prob_lo[i]) / ncell[i];
             prob_lo_in[i] = prob_lo[i];
          }
	 
	  std::cout << "#TK ncell" << ncell[0] << " " << ncell[1] << " " << ncell[2] << std::endl; 
  
          // We set this to -1 so that we can now control max_lev from the inputs file
	  int maxLevel = -1;

	  amrptr = new Amr(&rb,maxLevel,ncell,coord_sys);

	  Real strt_time = 0.0; 
	  Real stop_time = 1.0;

	  // This init call will cache the distribution map as determined by procMap                                                                                            // so that all data will end up on the right processor
	  amrptr->InitializeInit(strt_time, stop_time, &lev0_grids, &procMap);

          BoundaryDomain* bd = new BoundaryDomain(nr,hr);
	  amrptr->setBoundaryGeometry(bd->GetIntersectLoX(), bd->GetIntersectHiX(),
                                      bd->GetIntersectLoY(), bd->GetIntersectHiY());

          std::vector<double> x(3);
	  std::vector<double> attr(11);

	  for (size_t i=0; i<Track::block->bunch->getLocalNum(); i++) 
          {
            // X, Y, Z are stored separately from the other attributes
	    for (unsigned int k=0; k<3; k++)
	        x[k] = Track::block->bunch->R[i](k);

	    //  This allocates 11 spots -- 1 for Q, 3 for v, 3 for E, 3 for B, 1 for ID.
	    //  IMPORTANT: THIS ORDERING IS ASSUMED WHEN WE FILL E AT THE PARTICLE LOCATION
	    //             IN THE MOVEKICK CALL -- if Evec no longer starts at component 4
	    //             then you must change "start_comp_for_e" in Accel_advance.cpp
	    /*
	      Q      : 0
	      Vvec   : 1, 2, 3 the velocity 
	      Evec   : 4, 5, 6 the electric field at the particle location
	      Bvec   : 7, 8, 9 the electric field at the particle location
	      id+1   : 10 (we add 1 to make the particle ID > 0)
	    */
	
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

	    amrptr->addOneParticle(particle_id, Ippl::myNode(), x, attr);
	  }

          // It is essential that we call this routine since the particles
          //    may not currently be defined on the same processor as the grid
          //    that will hold them in the AMR world.
          amrptr->RedistributeParticles();

	  // This part of the call must come after we add the particles
	  // since this one calls post_init which does the field solve. 
	  amrptr->FinalizeInit(strt_time, stop_time);

	  amrptr->writePlotFile();
	
	  *gmsg << "A M R Initialization DONE" << endl;
	
	// }
//      }
#endif

        if(!OPAL->inRestartRun()) {
            if(!OPAL->hasDataSinkAllocated() && !Options::scan) {
                OPAL->setDataSink(new DataSink());
            } else if(Options::scan) {
                ds = OPAL->getDataSink();
                if(ds)
                    delete ds;
                OPAL->setDataSink(new DataSink());
            }
        } else {
            OPAL->setDataSink(new DataSink(OPAL->getRestartStep() + 1));
        }

        ds = OPAL->getDataSink();

        if(OPAL->hasBunchAllocated() && Options::scan)
            ds->reset();

        if(!OPAL->hasBunchAllocated() || Options::scan) {
            if(!mpacflg) {
                *gmsg << *dist << endl;
            } else {
                *gmsg << "* Multipacting flag is true. The particle distribution in the run command will be ignored " << endl;
            }
        }

        *gmsg << *beam << endl;
        *gmsg << *fs   << endl;
        *gmsg << *Track::block->bunch  << endl;

        *gmsg << "Phase space dump frequency " << Options::psDumpFreq << " and "
              << "statistics dump frequency " << Options::statDumpFreq << " w.r.t. the time step." << endl;
#ifdef HAVE_AMR_SOLVER
        itsTracker = new ParallelTTracker(*Track::block->use->fetchLine(),
                                          dynamic_cast<PartBunch &>(*Track::block->bunch), *ds,
                                          Track::block->reference, false, false, Track::block->localTimeSteps,
                                          Track::block->zstop, Track::block->timeIntegrator,
					  beam->getNumberOfParticles(),amrptr);
#else
        itsTracker = new ParallelTTracker(*Track::block->use->fetchLine(),
                                          dynamic_cast<PartBunch &>(*Track::block->bunch), *ds,
                                          Track::block->reference, false, false, Track::block->localTimeSteps,
                                          Track::block->zstop, Track::block->timeIntegrator,
					  beam->getNumberOfParticles());
#endif
        itsTracker->setMpacflg(mpacflg); // set multipacting flag in ParallelTTracker

    } else if(method == "PARALLEL-Z") {
        *gmsg << "  method == \"PARALLEL-Z\"" << endl;

    } else if(method == "CYCLOTRON-T") {
        OpalData::getInstance()->setInOPALCyclMode();
        Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

        fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));
        fs->initCartesianFields();
        Track::block->bunch->setSolver(fs);

	if (fs->hasPeriodicZ())
	  Track::block->bunch->setBCForDCBeam();
	else
	  Track::block->bunch->setBCAllOpen();

        Track::block->bunch->PType = 0;

        dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));

        // set macromass and charge for simulation particles
        double macromass = 0.0;
        double macrocharge = 0.0;

        const int specifiedNumBunch = int(std::abs(Round(Attributes::getReal(itsAttr[TURNS]))));

        if(beam->getNumberOfParticles() < 3 || beam->getCurrent() == 0.0) {
            macrocharge = beam->getCharge() * q_e;
            macromass = beam->getMass();
            dist->CreateOpalCycl(*Track::block->bunch, beam->getNumberOfParticles(), Options::scan);

        } else {

            /**
                   getFrequency() gets RF frequency [MHz], NOT isochronous  revolution frequency of particle!
                   getCurrent() gets beamcurrent [A]

	    */
            macrocharge = beam->getCurrent() / (beam->getFrequency() * 1.0e6); // [MHz]-->[Hz]

            if(!OPAL->hasBunchAllocated()) {
                if(!OPAL->inRestartRun()) {
                    macrocharge /= beam->getNumberOfParticles();
                    macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
                    dist->CreateOpalCycl(*Track::block->bunch, beam->getNumberOfParticles(), Options::scan);

                } else {
		    dist->DoRestartOpalCycl(*Track::block->bunch, beam->getNumberOfParticles(),
					    OPAL->getRestartStep(), specifiedNumBunch);
                    macrocharge /= beam->getNumberOfParticles();
                    macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
                }
            } else if(OPAL->hasBunchAllocated() && Options::scan) {
                macrocharge /= beam->getNumberOfParticles();
                macromass = beam->getMass() * macrocharge / (beam->getCharge() * q_e);
                dist->CreateOpalCycl(*Track::block->bunch, beam->getNumberOfParticles(), Options::scan);
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
        Track::block->bunch->calcBeamParameters_cycl();

        if(!OPAL->inRestartRun())
            if(!OPAL->hasDataSinkAllocated()) {
                ds = new DataSink();
                OPAL->setDataSink(ds);
            } else
                ds = OPAL->getDataSink();
        else {
            ds = new DataSink(OPAL->getRestartStep() + 1);
            OPAL->setDataSink(ds);
        }

        if(OPAL->hasBunchAllocated() && Options::scan)
            ds->reset();

        if(!OPAL->hasBunchAllocated() && !Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == CYCLOTRON-T, NEW TRACK" << endl;

            *gmsg << "* ********************************************************************************** " << endl;
        } else if(OPAL->hasBunchAllocated() && !Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == CYCLOTRON-T, FOLLOWUP TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        } else if(OPAL->hasBunchAllocated() && Options::scan) {
            *gmsg << "* ********************************************************************************** " << endl;
            *gmsg << "  Selected Tracking Method == CYCLOTRON-T, SCAN TRACK" << endl;
            *gmsg << "* ********************************************************************************** " << endl;
        }
        *gmsg << "* Number of neighbour bunches= " << specifiedNumBunch << endl;
        *gmsg << "* DT                         = " << Track::block->dT << endl;
        *gmsg << "* MAXSTEPS                   = " << Track::block->localTimeSteps << endl;
        *gmsg << "* Phase space dump frequency = " << Options::psDumpFreq << endl;
        *gmsg << "* Statistics dump frequency  = " << Options::statDumpFreq << " w.r.t. the time step." << endl;
        *gmsg << "* ********************************************************************************** " << endl;

        itsTracker = new ParallelCyclotronTracker(*Track::block->use->fetchLine(),
                 dynamic_cast<PartBunch &>(*Track::block->bunch), *ds, Track::block->reference,
                false, false, Track::block->localTimeSteps, Track::block->timeIntegrator);

        itsTracker->setNumBunch(specifiedNumBunch);

	if(OPAL->inRestartRun()) {
	  itsTracker->setPr(dist->GetPr());
	  itsTracker->setR(dist->GetR());
	  itsTracker->setTheta(dist->GetTheta());
	  itsTracker->setBeGa(dist->GetBeGa());
	}

        // statistical data are calculated (rms, eps etc.)
        Track::block->bunch->calcBeamParameters_cycl();

        *gmsg << *dist << endl;
        *gmsg << *beam << endl;
        *gmsg << *fs   << endl;
	//       *gmsg << *Track::block->bunch  << endl;

        if(specifiedNumBunch > 1) {

            // only for regular  run of multi bunches, instantiate the  PartBins class
            // note that for restart run of multi bunches, PartBins class is instantiated in function doRestart_cycl()
            if(!OPAL->inRestartRun()) {

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

            if(OPAL->inRestartRun()) {

                itsTracker->setLastDumpedStep(OPAL->getRestartStep());

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

    } else {
        throw OpalException("TrackRun::execute()",
                            "Method name \"" + method + "\" unknown.");
    }

    if((method != "PARALLEL-Z") && (method != "PARALLEL-T") && (method != "CYCLOTRON-T") && (method != "PARALLEL-SLICE")) {
        /*
          OLD SERIAL STUFF
        */
        // Open output file.
        string file = Attributes::getString(itsAttr[FNAME]);
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
        OPAL->setRestartRun(false);
    }
    //  delete [] fs;
    //  delete [] dist;

    OPAL->bunchIsAllocated();
    if(method == "PARALLEL-SLICE")
        OPAL->slbunchIsAllocated();

    delete itsTracker;

}

double TrackRun::SetDistributionParallelT(Beam *beam) {

    // If multipacting flag is not set, get distribution(s).
    if (!Attributes::getBool(itsAttr[MULTIPACTING])) {
        /*
         * Distribution(s) can be set via a single distribution or a list
         * (array) of distributions. If an array is defined the first in the
         * list is treated as the primary distribution. All others are added to
         * it to create the full distribution.
         */
        std::vector<std::string> distributionArray
            = Attributes::getStringArray(itsAttr[DISTRIBUTIONS]);

        if (distributionArray.size() > 0) {
            *gmsg << endl
                  << "---------------------------------" << endl
                  << "Found more than one distribution:" << endl << endl;

            for (std::vector<std::string>::const_iterator distIterator
                 = distributionArray.begin();
                 distIterator != distributionArray.end(); ++distIterator) {

                if (distIterator == distributionArray.begin()) {
                    dist = Distribution::find(*distIterator);
                    *gmsg << "Main Distribution" << endl
                          << "-----------------" << endl
                          << *distIterator << endl << endl
                          << "Secondary Distribution(s)" << endl
                          << "-------------------------" << endl;
                } else {
                    Distribution *distribution = Distribution::find(*distIterator);
                    *gmsg << *distIterator << endl;
                    distrs_m.push_back(distribution);
                }
            }
            *gmsg << endl
                  << "---------------------------------" << endl << endl;
        } else
            dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));
    }


    /*
     * Initialize distributions.
     */
    size_t numberOfParticles = beam->getNumberOfParticles();
    if (!OPAL->hasBunchAllocated()) {
        if (!OPAL->inRestartRun()) {
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

                OPAL->setGlobalPhaseShift(0.5 * dist->GetTEmission());
            }
        } else {
            /*
             * Read in beam from restart file.
             */
            dist->DoRestartOpalT(*Track::block->bunch, numberOfParticles, OPAL->getRestartStep());
        }
    } else if (OPAL->hasBunchAllocated() && Options::scan) {
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

        OPAL->setGlobalPhaseShift(0.5 * dist->GetTEmission());
    }

    // Return charge per macroparticle.
    return beam->getCharge() * beam->getCurrent() / beam->getFrequency() / numberOfParticles;

}

ParallelTTracker *TrackRun::setupForAutophase() {

    Inform m("setupForAutophase() ");

    Beam *beam = Beam::find(Attributes::getString(itsAttr[BEAM]));

    DataSink *ds = NULL;


    fs = FieldSolver::find(Attributes::getString(itsAttr[FIELDSOLVER]));

    fs->initCartesianFields();

    Track::block->bunch->setSolver(fs);

    dist = Distribution::find(Attributes::getString(itsAttr[DISTRIBUTION]));

    double charge = beam->getCharge() * beam->getCurrent() / beam->getFrequency();

    charge /= beam->getNumberOfParticles();

    Track::block->bunch->setdT(Track::block->dT);

    Track::block->bunch->resetIfScan();

    Track::block->bunch->LastSection = 1;

    Track::block->bunch->setCharge(charge);

    Track::block->bunch->setChargeZeroPart(charge);// just set bunch->qi_m=charge, don't set bunch->Q[] as we have not initialized any particle yet.

    double coefE = 1.0 / (4 * pi * epsilon_0);
    Track::block->bunch->setCouplingConstant(coefE);
    size_t N = 1;
#ifdef HAVE_AMR_SOLVER
    Amr* null_amrptr = 0;
    return new ParallelTTracker(*Track::block->use->fetchLine(),
                                dynamic_cast<PartBunch &>(*Track::block->bunch), *ds,
                                Track::block->reference, false, false, Track::block->localTimeSteps,
                                Track::block->zstop, Track::block->timeIntegrator,
				N,null_amrptr);
#else
    return new ParallelTTracker(*Track::block->use->fetchLine(),
                                dynamic_cast<PartBunch &>(*Track::block->bunch), *ds,
                                Track::block->reference, false, false, Track::block->localTimeSteps,
                                Track::block->zstop, Track::block->timeIntegrator,
				N);
#endif
}

#ifdef HAVE_AMR_SOLVER
std::vector<std::string> TrackRun::filterString(std::string str) {
  /*
    pid tokens[0]
    il == tokens[2], ih == tokens[3]
    jl == tokens[5], jh == tokens[6]
    kl == tokens[9], kh == tokens[9]
   */

  // charakters to remove from the string
  char chars[] = "Node=;vn_mDmain{[][][]}:,";
 	
  for (unsigned int i = 0; i < strlen(chars); ++i) 
    std::replace(str.begin(), str.end(), chars[i], ' ');

  std::vector<std::string> tokens;
  
  // filter spaces
  std::istringstream iss(str);
  copy(std::istream_iterator<std::string>(iss),
       std::istream_iterator<std::string>(),
       std::back_inserter<std::vector<std::string> >(tokens));
  return tokens;
}

std::pair<Box,unsigned int> TrackRun::getBlGrids(std::string str){

  std::vector<std::string> tokens = filterString(str);

  unsigned int theGrid; 
  std::istringstream (tokens[0]) >> theGrid;

  int ilo,ihi,jlo,jhi,klo,khi;

  std::istringstream (tokens[2]) >> ilo;
  std::istringstream (tokens[5]) >> jlo;
  std::istringstream (tokens[8]) >> klo;
  std::istringstream (tokens[3]) >> ihi;
  std::istringstream (tokens[6]) >> jhi;
  std::istringstream (tokens[9]) >> khi;

  Inform m2a("AMR ",INFORM_ALL_NODES);
  /*
  m2a << "Grid " << tokens[0]
      << " i (" << tokens[2] << " ... " << tokens[3] << ")" 
      << " j (" << tokens[5] << " ... " << tokens[6] << ")" 
      << " k (" << tokens[8] << " ... " << tokens[9] << ")" << " myNode " << Ippl::myNode() << endl;
  */
  IntVect loEnd(ilo,jlo,klo);
  IntVect hiEnd(ihi,jhi,khi);
  Box bx(loEnd,hiEnd);

  return std::pair<Box,unsigned int>(bx,theGrid);
}
#endif

