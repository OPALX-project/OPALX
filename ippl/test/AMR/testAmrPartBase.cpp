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
#include "PartBunchAmr.h"

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase<amrplayout_t> amrbase_t;
typedef PartBunchAmr<amrplayout_t> amrbunch_t;

typedef std::deque<Particle<4,0> > PBox4;
typedef typename std::map<int,PBox4> PMap4;

typedef Array<std::unique_ptr<MultiFab> > container_t;



void createRandomParticles(amrbunch_t *pbase, int N, int myNode, int seed = 1) {

  srand(seed);
  for (int i = 0; i < N; ++i) {
    pbase->createWithID(myNode * N + i + 1);
    pbase->qm[i] = (double)rand() / RAND_MAX;

    pbase->R[i][0] = (double)rand() / RAND_MAX;
    pbase->R[i][1] = (double)rand() / RAND_MAX;
    pbase->R[i][2] = (double)rand() / RAND_MAX;  
  }
    
}

void createRandomParticles(ParticleContainer<4,0> *pc, int N, int myNode, int seed = 1) {

  srand(seed);

  for (int i = 0; i < N; ++i) {

    std::vector<double> attrib;
    double r = (double)rand() / RAND_MAX;
    attrib.push_back(r);
    attrib.push_back(0.0);
    attrib.push_back(0.0);
    attrib.push_back(0.0);

    std::vector<double> xloc;
    xloc.push_back((double)rand() / RAND_MAX);
    xloc.push_back((double)rand() / RAND_MAX);
    xloc.push_back((double)rand() / RAND_MAX);
    pc->addOneParticle(myNode * N + i + 1, myNode, xloc, attrib);
  }

}

void writeAscii(amrbunch_t *pbase, int N, int myNode) {
  std::ofstream myfile;
  std::string fname = "Ippl-";
  fname += std::to_string(myNode);
  fname += ".dat"; 
  myfile.open(fname);

  myfile << "id\tR0\tR1\tR2\tlevel\tgrid\tE0\tE1\tE2\n";
  for (size_t i = 0; i < pbase->getLocalNum(); i++) {
    myfile << std::setprecision(3) << pbase->ID[i] << "\t" << pbase->R[i][0] 
	   << "\t" << pbase->R[i][1] << "\t" << pbase->R[i][2] 
	   << "\t" << pbase->m_lev[i] << "\t" << pbase->m_grid[i]
	   << "\t" << pbase->qm[i] << "\t" << pbase->E[i][0] 
	   << "\t" << pbase->E[i][1] << "\t" << pbase->E[i][2] << "\n";
  }

  myfile.close();

}

void writeAscii(ParticleContainer<4,0> *pc, int N, size_t nLevels, int myNode) {
  std::ofstream myfile;
  std::string fname = "BoxLib-";
  fname += std::to_string(myNode);
  fname += ".dat"; 
  myfile.open(fname);

  myfile << "id\tR0\tR1\tR2\tlevel\tgrid\tqm\tE0\tE1\tE2\n";
  for (unsigned int lev = 0; lev < nLevels; lev++) {
    const PMap4& pmap = pc->GetParticles(lev);

    for (const auto& kv : pmap) {
      const PBox4& pbox = kv.second;

      for (auto it = pbox.cbegin(); it != pbox.cend(); ++it) {
	myfile << std::setprecision(3) << it->m_id << "\t" << it->m_pos[0] 
	       << "\t" << it->m_pos[1] << "\t" 
	       << it->m_pos[2] << "\t" << it->m_lev << "\t" 
	       << it->m_grid  << "\t" << it->m_data[0] << "\t" 
	       << it->m_data[1] << "\t" << it->m_data[2] << "\t" << it->m_data[3] << "\n";
      }
    }
  }

  myfile.close();

}

void readData(std::string fname, std::vector<int> &data) {
  std::ifstream myfile;
  myfile.open(fname);

  std::string line;
  getline(myfile, line); //skip first line, since its a header
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
  int match = 0;
  for (unsigned i = 0; i < ippldata.size(); ++i) {
    if (ippldata[i] != bldata[i]) {
      std::cout << "===ERROR=== Particle distribution on node " << node << " doesn't match!" 
		<< std::endl; 
      match = 0;
      return;  
    }
  }

  int g_match = 0;
  MPI_Reduce(&match, &g_match, 1, MPI_INT, MPI_SUM, 0, Ippl::getComm());

  if (Ippl::myNode() == 0 && g_match == 0)
    std::cout << "Particle distribution for Ippl and BoxLib matches" << std::endl;
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
	    Array<DistributionMapping> &dmap, Array<int> &rr, 
	    size_t nLevels, int myNode, 
	    PArray<MultiFab> &field, PArray<MultiFab> &efield,
	    int N, int seed) 
{

  static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
  IpplTimings::startTimer(mainTimer);

  //create a new layout using ParticleAmrLayout class
  amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);

  //create a particle bunch
  PartBunchAmr<amrplayout_t>* pbase = new PartBunchAmr<amrplayout_t>();
  pbase->initialize(PL);
  pbase->initializeAmr();

  //create N random particles on each core
  createRandomParticles(pbase, N, myNode, seed);

  //update redistributes particles among the cores
  pbase->update();

  //call assign density to scatter the paarticle attribute qm on the grid
  pbase->setAllowParticlesNearBoundary(true);
  pbase->AssignDensitySingleLevel(pbase->qm, field[0], 0);
  pbase->AssignDensity(pbase->qm, false, field, 0, 1);

  //copy the valies from field to all the components of efield
  
  for (size_t lev = 0; lev < nLevels; ++lev) {
    efield[lev].setVal(0.0);
    MultiFab::Copy(efield[lev], field[lev], 0, 0, 1, 0);
    MultiFab::Copy(efield[lev], field[lev], 0, 1, 1, 0);
    MultiFab::Copy(efield[lev], field[lev], 0, 2, 1, 0);
  }

  //get values from grid to particles
  pbase->GetGravity(pbase->E, efield);

  //write the particles on the core to file - one file per core created
  writeAscii(pbase, N, myNode);

  delete pbase;
  
  IpplTimings::stopTimer(mainTimer);

}

void doBoxLib(Array<Geometry> &geom, Array<BoxArray> &ba, 
	      Array<DistributionMapping> &dmap, Array<int> &rr, 
	      size_t nLevels, int myNode, 
	      PArray<MultiFab> &field, PArray<MultiFab> &efield,
	      int N, int seed) 
{

  //create new BoxLib particle container
  ParticleContainer<4,0> *pc = new ParticleContainer<4,0>(geom, dmap, ba, rr);
  pc->SetVerbose(0);

  //create N random particles on each core
  createRandomParticles(pc, N, myNode, seed);

  //redistribute particles among the cores
  pc->Redistribute();

  //call assign density to scatter the paarticle attribute qm on the grid
  pc->SetAllowParticlesNearBoundary(true);
  pc->AssignDensitySingleLevel(0, field[0], 0, 0);
  pc->AssignDensity(0, false, field, 0, 1, 1);

  
  //copy the valies from field to all the components of efield
  for (size_t lev = 0; lev < nLevels; ++lev) {
    efield[lev].setVal(0.0);
    MultiFab::Copy(efield[lev], field[lev], 0, 0, 1, 0);
    MultiFab::Copy(efield[lev], field[lev], 0, 1, 1, 0);
    MultiFab::Copy(efield[lev], field[lev], 0, 2, 1, 0);
  }

  
  //loop trough all the levels
  for (size_t lev = 0; lev < nLevels; ++lev) {
    //get grids on the level
    PMap4& pmap = pc->GetParticles(lev);

    //loop trough grids on the level
    for (auto& kv : pmap) {
      int grid = kv.first;
      PBox4& pbox = kv.second;
      const int n = pbox.size();
      const FArrayBox& gfab = efield[lev][grid];

      //loop trough the particles in the grid and call GetGravity for each
      //assign grav to particle data components 1,2,3
      for (int i = 0; i < n; ++i) {
	Particle<4,0> &p = pbox[i];
	Real grav[BL_SPACEDIM];

	ParticleBase::GetGravity(gfab, geom[lev], p, grav);

	p.m_data[1] = grav[0];
	p.m_data[2] = grav[1];
	p.m_data[3] = grav[2];
      }

    }
  }
  
  //write the particles on the core to file - one file per core created
  writeAscii(pc, N, nLevels, myNode);

}

int main(int argc, char *argv[]) {
    
  Ippl ippl(argc, argv);   
  Inform msg("AMRParticle");

  //get command line arguments - number of particles and number of times the test is run
  int L = 1;  //how many times the test should be run
  int N = 10; //number of particles in tests
  if (argc > 1) {
    for (int arg = 1; arg < argc; ++arg) {
      if (argv[arg] == std::string("-N"))
	N = atoi(argv[arg+1]);
      
      if (argv[arg] == std::string("-loop"))
	L = atoi(argv[arg+1]);
    }
  }

  std::cout << "Start test with " << N << " particles and run " << L << " times." << std::endl;

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
  
  //print out the geomety
  /*
  for (size_t i = 0; i < nLevels; ++i) {
    std::cout << "Level: " << i << std::endl;
    std::cout << geom[i] << std::endl;
    std::cout << std::endl;
    std::cout << ba[i] << std::endl;
  }
  */

  //create a multifabs one is used with BoxLib tests, one for Ippl tests
  PArray<MultiFab> field_ippl;
  field_ippl.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    field_ippl.set(lev, new MultiFab(ba[lev], 1, 1, dmap[lev]));

  PArray<MultiFab> field_bl;
  field_bl.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    field_bl.set(lev, new MultiFab(ba[lev], 1, 1, dmap[lev]));

  PArray<MultiFab> efield;
  efield.resize(nLevels);
  for (size_t lev = 0; lev < nLevels; ++lev)
    efield.set(lev, new MultiFab(ba[lev], BL_SPACEDIM, 1, dmap[lev]));


  //Do ippl and boxlib runs multiple times.
  //At each step N particles are created at random location inside the domain.
  //Update is called to distribute the particles ammong the processes.
  //Each particle container contains additional attribute which also is assigned a random value.
  //AssignDensity is used to scatter this attribute on the grids field_ippl and field_bl.
  //The particles used for both tests are the same - same locations and attribute values.
  //compareDistribution checks if the particles are the same on each core after update
  //for BoxLib and Ippl particle containers.
  //Compare fields check if field_ippl and field_bl are the same after AssignDensity
  for (int i = 0; i < L; ++i) {

    doIppl(geom, ba, dmap, rr, nLevels, Ippl::myNode(), field_ippl, efield, N, i);
    doBoxLib(geom, ba, dmap, rr, nLevels, Ippl::myNode(), field_bl, efield, N, i);

    if (Ippl::myNode() == 0)
      std::cout << "Results for test " << i + 1 << std::endl;
    compareDistribution(Ippl::myNode());
    compareFields(field_ippl, field_bl, Ippl::myNode());

    if (Ippl::myNode() == 0)
      std::cout << std::endl;
  }

  IpplTimings::print();

  return 0;
}
