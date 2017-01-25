// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 *
 ***************************************************************************/

/***************************************************************************

This test program sets up a simple sine-wave electric field in 3D,
  creates a population of particles with random q/m values (charge-to-mass
  ratio) and velocities, and then tracks their motions in the static
  electric field using nearest-grid-point interpolation.

Usage:

 mpirun -np 4 testAmrPartBunch IPPL 32 32 32 100 10
 
 mpirun -np 4 testAmrPartBunch BOXLIB 32 32 32 100 10 0

***************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>


#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <ParmParse.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

#include "AmrParticleBase.h"
#include "ParticleAmrLayout.h"

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase<amrplayout_t> amrbase_t;

typedef std::deque<Particle<10,0> > PBox;
typedef std::map<int, std::deque<Particle<10,0> > > PMap;
typedef Array<std::unique_ptr<MultiFab> > container_t;


void createRandomParticles(amrbase_t *pbase, int N, int myNode) {

  srand(1);
  for (int i = 0; i < N; ++i) {
    pbase->createWithID(myNode * N + i + 1);
    pbase->qm[i] = (double)rand() / RAND_MAX;

    pbase->R[i][0] = (double)rand() / RAND_MAX;
    pbase->R[i][1] = (double)rand() / RAND_MAX;
    pbase->R[i][2] = (double)rand() / RAND_MAX;  
  }
    
}

void createRandomParticles(ParticleContainer<10,0> *pc, int N, int myNode) {

  srand(1);

  for (int i = 0; i < N; ++i) {

    std::vector<double> attrib;
    double r = (double)rand() / RAND_MAX;
    for (int i = 0; i < 10; i++) {
      attrib.push_back(r);
    }

    std::vector<double> xloc;
    xloc.push_back((double)rand() / RAND_MAX);
    xloc.push_back((double)rand() / RAND_MAX);
    xloc.push_back((double)rand() / RAND_MAX);
    pc->addOneParticle(myNode * N + i + 1, myNode, xloc, attrib);
  }

}

void writeAscii(AmrParticleBase<amrplayout_t> *pbase, int N, int myNode) {
  std::ofstream myfile;
  std::string fname = "Ippl-";
  fname += std::to_string(myNode);
  fname += ".dat"; 
  myfile.open(fname);

  for (size_t i = 0; i < pbase->getLocalNum(); i++) {
    myfile << std::setprecision(3) << pbase->ID[i] << "\t" << pbase->R[i][0] 
	   << "\t" << pbase->R[i][1] 
	   << "\t" << pbase->R[i][2] << "\t" << pbase->m_lev[i] 
	   << "\t" << pbase->m_grid[i] << "\n";
  }

  myfile.close();

}

void writeAscii(ParticleContainer<10,0> *pc, int N, size_t nLevels, int myNode) {
  std::ofstream myfile;
  std::string fname = "BoxLib-";
  fname += std::to_string(myNode);
  fname += ".dat"; 
  myfile.open(fname);

  for (unsigned int lev = 0; lev < nLevels; lev++) {
    const PMap& pmap = pc->GetParticles(lev);

    for (const auto& kv : pmap) {
      const PBox& pbox = kv.second;

      for (auto it = pbox.cbegin(); it != pbox.cend(); ++it) {
	myfile << std::setprecision(3) << it->m_id << "\t" << it->m_pos[0] 
	       << "\t" << it->m_pos[1] << "\t" 
	       << it->m_pos[2] << "\t" << it->m_lev << "\t" 
	       << it->m_grid  << "\n";
      }
    }
  }

  myfile.close();

}

void doIppl(Array<Geometry> &geom, Array<BoxArray> &ba, 
	    Array<DistributionMapping> &dmap, Array<int> &rr, int myNode, 
	    PArray<MultiFab> &chargeOnGrid) 
{

  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
  IpplTimings::startTimer(mainTimer);

  amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);

  AmrParticleBase<amrplayout_t>* pbase = new AmrParticleBase<amrplayout_t>();
  pbase->initialize(PL);
  pbase->initializeAmr();

  int N = 1e4;
  createRandomParticles(pbase, N, myNode);

  pbase->update();
  pbase->sort();

  chargeOnGrid[0].setVal(0.0);

  pbase->setAllowParticlesNearBoundary(true);
  pbase->AssignDensitySingleLevel(pbase->qm, chargeOnGrid[0], 0);

  std::cout << "Charge on grid: " << chargeOnGrid[0].sum() << std::endl;

  writeAscii(pbase, N, myNode);

  delete pbase;
  
  IpplTimings::stopTimer(mainTimer);

  IpplTimings::print();

}

void doBoxLib(Array<Geometry> &geom, Array<BoxArray> &ba, 
	      Array<DistributionMapping> &dmap, Array<int> &rr,
	      size_t nLevels, int myNode, PArray<MultiFab> &chargeOnGrid) 
{


  int N = 1e4;
  ParticleContainer<10,0> *pc = new ParticleContainer<10,0>(geom, dmap, ba, rr);

  createRandomParticles(pc, N, myNode);
  pc->Redistribute();

  chargeOnGrid[0].setVal(0.0);

  pc->SetAllowParticlesNearBoundary(true);
  pc->AssignDensitySingleLevel(0, chargeOnGrid[0], 0, 0);

  std::cout << "Charge on grid: " << chargeOnGrid[0].sum() << std::endl;

  writeAscii(pc, N, nLevels, myNode);

}

int main(int argc, char *argv[]) {
    
  Ippl ippl(argc, argv);
    
  Inform msg("AMRParticle");


  /* Setup BoxLib */
  BoxLib::Initialize(argc,argv, false);

  size_t nLevels = 2;
  size_t maxBoxSize = 8;

  //set up the geometry
  int n_cell = 16;
  IntVect low(0, 0, 0);
  IntVect high(n_cell - 1, n_cell - 1, n_cell - 1);    
  Box bx(low, high);

  //physical domain boundaries
  RealBox domain;
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    domain.setLo(i, 0.0);
    domain.setHi(i, 1.0);
  }

  RealBox fine_domain;
  for (int i = 0; i < BL_SPACEDIM; ++i) {
    fine_domain.setLo(i, 0.0);
    fine_domain.setHi(i, 0.5);
  }

  //periodic boundary conditions in all directions
  int bc[BL_SPACEDIM] = {0, 0, 0};

  //Container for geometry at all levels
  Array<Geometry> geom;
  geom.resize(nLevels);

  // Container for boxes at all levels
  Array<BoxArray> ba;
  ba.resize(nLevels);    

  // level 0 describes physical domain
  geom[0].define(bx, &domain, 0, bc);

  //refinement for each level
  Array<int> rr(nLevels - 1);
  for (unsigned int lev = 0; lev < rr.size(); ++lev)
    rr[lev] = 2;

  // geometries of refined levels
  for (unsigned int lev = 1; lev < nLevels; ++lev)
    geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(), rr[lev - 1]), &domain, 0, bc);
        
  // box at level 0
  ba[0].define(bx);
  ba[0].maxSize(maxBoxSize);

  //box at level 1
    
  //build boxes at finer levels
  if (nLevels > 1) {
    int n_fine = n_cell * rr[0];
    IntVect refined_lo(0, 0, 0);
    IntVect refined_hi(15, 15, 15);
    
    Box refined_box(refined_lo, refined_hi);
    ba[1].define(refined_box);
    ba[1].maxSize(maxBoxSize);
  }

  /*
   * distribution mapping
   */
  Array<DistributionMapping> dmap;
  dmap.resize(nLevels);
  dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
  if (nLevels > 1)
    dmap[1].define(ba[1], ParallelDescriptor::NProcs() /*nprocs*/);
   

  /* BoxLib geometry setup done */

  
  for (size_t i = 0; i < nLevels; ++i) {
    std::cout << "Level: " << i << std::endl;
    std::cout << geom[i] << std::endl;
    std::cout << std::endl;
    std::cout << ba[i] << std::endl;
  }


  //create a multifab
  PArray<MultiFab> chargeOnGrid;
  chargeOnGrid.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    chargeOnGrid.set(lev, new MultiFab(ba[lev], 1, 1, dmap[lev]));

    
  doIppl(geom, ba, dmap, rr, Ippl::myNode(), chargeOnGrid);

  doBoxLib(geom, ba, dmap, rr, nLevels, Ippl::myNode(), chargeOnGrid);

  return 0;
}
