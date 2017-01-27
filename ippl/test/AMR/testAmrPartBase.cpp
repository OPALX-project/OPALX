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

//typedef std::deque<Particle<1,0> > PBox;
//typedef std::map<int, std::deque<Particle<1,0> > > PMap;
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

void createRandomParticles(ParticleContainer<1,0> *pc, int N, int myNode) {

  srand(1);

  for (int i = 0; i < N; ++i) {

    std::vector<double> attrib;
    double r = (double)rand() / RAND_MAX;
    for (int i = 0; i < 1; i++) {
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

void writeAscii(ParticleContainer<1,0> *pc, int N, size_t nLevels, int myNode) {
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

void readData(std::string fname, std::vector<int> &data) {
  std::ifstream myfile;
  myfile.open(fname);

  std::string line;
  while ( getline(myfile, line) ) {
    std::istringstream ss(line);
    int id;
    ss >> id;
    data.push_back(id);
  }
}

void compareDistribution(int node) {

  //read the data files containing particle information
  std::vector<int> ippldata, bldata;
  readData("Ippl-" + std::to_string(node) + ".dat", ippldata);
  readData("BoxLib-" + std::to_string(node) + ".dat", bldata);

  //check if the size of particles per node is the same for Ippl and BoxLib versions
  if ( ippldata.size() != bldata.size() ) {
    std::cout << "===ERROR=== Particle distribution on node " << node << " doesn't match!" 
	      << std::endl; 
    return;
  }
  
  //sort the particle IDs
  std::sort(ippldata.begin(), ippldata.end());
  std::sort(bldata.begin(), bldata.end());

  //check if both nodes contain the same particles
  for (unsigned i = 0; i < ippldata.size(); ++i) {
    if (ippldata[i] != bldata[i]) {
      std::cout << "===ERROR=== Particle distribution on node " << node << " doesn't match!" 
		<< std::endl; 
      return;  
    }
  }

  std::cout << "Particle distribution for Ippl and BoxLib matches for node " << node 
	    << std::endl;
}

void compareFields(PArray<MultiFab> &field_ippl, PArray<MultiFab> &field_bl, int node) {

  bool fields_match = true;
  double ippl_sum, bl_sum;

  for (int lev = 0; lev < field_ippl.size(); ++lev) {
    //calculate the sum of all the components in multifab
    ippl_sum = field_ippl[lev].sum();
    bl_sum = field_bl[lev].sum();

    //check if the sums are the same for Ippl and BoxLib
    //only node 0 prints the error since the sum is the same on all nodes
    if ( abs( ippl_sum - bl_sum) > 1e-6 && node == 0) {
      std::cout << "===ERROR=== Fields don't match on level " << lev 
		<< ": " << ippl_sum << " != " << bl_sum << std::endl; 
      
      fields_match = false;
    }
  }

  if (fields_match && node == 0)
    std::cout << "Fields match on all levels for BoxLib and Ippl AssignDensity" << std::endl;
  
}

void doIppl(Array<Geometry> &geom, Array<BoxArray> &ba, 
	    Array<DistributionMapping> &dmap, Array<int> &rr, int myNode, 
	    PArray<MultiFab> &field,
	    int N) 
{

  static IpplTimings::TimerRef createTimer = IpplTimings::getTimer("AMR create particles");
  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
  IpplTimings::startTimer(mainTimer);

  amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);

  AmrParticleBase<amrplayout_t>* pbase = new AmrParticleBase<amrplayout_t>();
  pbase->initialize(PL);
  pbase->initializeAmr();

  IpplTimings::startTimer(createTimer);
  createRandomParticles(pbase, N, myNode);
  IpplTimings::stopTimer(createTimer);

  pbase->update();
  pbase->sort();

  pbase->setAllowParticlesNearBoundary(true);
  pbase->AssignDensitySingleLevel(pbase->qm, field[0], 0);

  pbase->AssignDensity(pbase->qm, false, field, 0, 1);

  writeAscii(pbase, N, myNode);

  delete pbase;
  
  IpplTimings::stopTimer(mainTimer);

}

void doBoxLib(Array<Geometry> &geom, Array<BoxArray> &ba, 
	      Array<DistributionMapping> &dmap, Array<int> &rr,
	      size_t nLevels, int myNode, PArray<MultiFab> &field,
	      int N) 
{


  ParticleContainer<1,0> *pc = new ParticleContainer<1,0>(geom, dmap, ba, rr);

  createRandomParticles(pc, N, myNode);
  pc->Redistribute();

  pc->SetAllowParticlesNearBoundary(true);
  pc->AssignDensitySingleLevel(0, field[0], 0, 0);

  pc->AssignDensity(0, false, field, 0, 1, 1);

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
  int bc[BL_SPACEDIM] = {1, 1, 1};

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
  PArray<MultiFab> field_ippl;
  field_ippl.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    field_ippl.set(lev, new MultiFab(ba[lev], 1, 1, dmap[lev]));

  PArray<MultiFab> field_bl;
  field_bl.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    field_bl.set(lev, new MultiFab(ba[lev], 1, 1, dmap[lev]));

   
  int N = 1e5;
  doIppl(geom, ba, dmap, rr, Ippl::myNode(), field_ippl, N);
  doBoxLib(geom, ba, dmap, rr, nLevels, Ippl::myNode(), field_bl, N);

  compareDistribution(Ippl::myNode());
  compareFields(field_ippl, field_bl, Ippl::myNode());

  IpplTimings::print();

  return 0;
}
