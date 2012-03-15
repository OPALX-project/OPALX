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
 * Visit http://www.
 *
 ***************************************************************************/

/***************************************************************************
 
 Tests the FEM multigrid solver
 
   test1: tests interpolation 

   test2: test poisson solver ....

 Usage:                            
  
 mpirun -np 2 ./femmgTest data1.dat 10000 OOO [test1 | test2]
 
***************************************************************************/

#include "Ippl.h"
#include "extpde.h"
#include "Physics.hh"
#include "Algorithms/PartData.h"

#include "GTConst.hh"
#include "GTConfigure.hh"
#include "GTElemData.hh"
#include "GTChargedParticles.hh"
#include "GTMGFem.hh"

#include <new>
#include <exception>
#include <vector>
#include <map>
#include <set> 
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <fstream>

#include "mpi.h"


// dimension of our positions
const unsigned Dim = 3;

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

void createDistr(ChargedParticles<playout_t> *bunch,unsigned int totalP,    
		 double globalMin[3],
		 double globalMax[3],
		 double r,
		 double R1,
		 double charge)
{
  unsigned num_my_particles = 0;
  double qi = charge/totalP;
  
  for (int i = 0; i < totalP; ++ i) {
    // Generate random position inside domain.
    // Here we rely on the fact that all processors generate the same sequence
    // of random numbers.
    
    Vektor<double,3> pos = Vektor<double,3>(2.0*r*(0.5 - IpplRandom()),
					    2.0*r*(0.5 - IpplRandom()),
					    2.0*r*(0.5 - IpplRandom()));
    
    while (!bunch->getLayout().is_in_domain(pos) || !(dot(pos,pos)<=r*r)) {
      pos = Vektor<double,3>(2.0*r*(0.5 - IpplRandom()),
			     2.0*r*(0.5 - IpplRandom()),
			     2.0*r*(0.5 - IpplRandom()));
    }
    
    // Each processor tests if pos is inside local subdomain.
    if (bunch->getLayout().is_local_pos(pos)) {
      // Create particle in container.
      bunch->create(1);
      bunch->R[num_my_particles] = pos;
      bunch->P[num_my_particles] = Vektor<double,3>(0.0);
      bunch->E[num_my_particles] = Vektor<double,3>(0.0);
      bunch->Phi[num_my_particles] = 0.0;
      bunch->Q[num_my_particles] = qi;
      num_my_particles += 1;
    }
  }

  /*
    Create some particles exactly through the space to show 
    the Phi and the E-Field
  */    
  
  int numLinePart = 10000;
  double nR = 2*R1/(numLinePart);

  Vektor<double,3> x = Vektor<double,3>(-R1, 0, 0);
  for (int i=0; i<numLinePart; i++) {
    x += Vektor<double,3>(nR, 0, 0);
    if (bunch->getLayout().is_local_pos(x)) {
      bunch->create(1);
      bunch->R[num_my_particles] = x;
      bunch->Q[num_my_particles] = 0.0;
      bunch->P[num_my_particles] = Vektor<double,3>(-1.0);
      bunch->E[num_my_particles] = Vektor<double,3>(0.0);
      bunch->Phi[num_my_particles] = 0.0;
      num_my_particles++; 
    }
  }

  x = Vektor<double,3>(0.0, 0.0, -R1);
  for (int i=0; i<numLinePart; i++) {
    x += Vektor<double,3>(0, 0, nR);
    if (bunch->getLayout().is_local_pos(x)) {
      bunch->create(1);
      bunch->R[num_my_particles] = x;
      bunch->Q[num_my_particles] = 0.0;
      bunch->P[num_my_particles] = Vektor<double,3>(-1.0);
      bunch->E[num_my_particles] = Vektor<double,3>(0.0);
      bunch->Phi[num_my_particles] = 0.0;
      num_my_particles++; 
    }
  }
}

int main(int argc, char *argv[]) {
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);
  IpplRandom.SetSeed(static_cast<unsigned long>(234244113131719.0));
  
  static IpplTimings::TimerRef mainTimer       = IpplTimings::getTimer("mainTimer");
  static IpplTimings::TimerRef partCreaTimer   = IpplTimings::getTimer("partCreaTimer");
   
  IpplTimings::startTimer(mainTimer);
    
  const unsigned int totalP = atoi(argv[2]); 

  BCT myBC;
  
  if (string(argv[2])==string("OOO")) {
    myBC = OOO; // open boundary
  } else if (string(argv[2])==string("OOP")) {
    myBC = OOP; // open boundary in x and y, periodic in z
  } else {
    myBC = PPP; // all periodic
  }
  

  bool test1 = false;
  bool test2 = false;

  if (string(argv[3])==string("test1")) 
    test1 = true;
  else if (string(argv[3])==string("test2")) 
    test2 = true;
  else 
    ERRORMSG("selected test not implemented, no test will be executed");

  double grid_stretch = 1.0;
  int r_parallel  = 1;
  double o_square = 1.0;
  double s_square = 1.0;
  double meshsize;
  int iteration;
  int iter_relax;
  int p_test;
  double a,b,l;

  double globalMin[3];
  double globalMax[3];
  double localMin[3];
  double localMax[3];

  ifstream fstr;
  fstr.open(argv[1],std::ios::in);
  fstr >> meshsize  >> iteration >> iter_relax
       >> r_parallel  >> o_square >> s_square >> p_test
       >> a >>  b >> l;
  fstr.close();

  double r = 1.0;                                             // radius of \Omega
  double R1 = 0.25;                                           // radius of charges
  double charge = 200.0e-9 / (4.0*pi*Physics::epsilon_0);         
  double dirBC = -charge/r; // (4.0*pi*Physics::epsilon_0*r);  
    
  /*
    Die Standartkugel hat radius 0.5 und den 
    Mittelpunkt bei (x,y,z) = (0.5,0.5,0.5)
    
  */      
    
  D3vector Vs(r,r,r);
  D3vector Vm(r/2,r/2,r/2);  
  D3vector Sc(-0.5,-0.5,-0.5); 
  double   st2 = grid_stretch * grid_stretch;
  D3vector Vstretch(1.0/grid_stretch,1.0/grid_stretch,1.0);
    
  /*
    Ball r=1, center=(0,0,0)
  */

  Ball ball;
  Domain*  geom_domain = new Domain((&ball+Sc)*Vstretch*2*r);
    
  Grid_gen_parameters ggen_parameter;
  ggen_parameter.Set_r_parallel(r_parallel);
  ggen_parameter.Set_offset_square(o_square);
  ggen_parameter.Set_stretch_square(s_square);
    
  Grid *grid = new Grid(meshsize, geom_domain, ggen_parameter, MPI_COMM_WORLD); // Construction of grid

  ofstream fo;
  if(Ippl::myNode()==0) fo.open("data/regionlayout.inp",ios :: out);
  grid->Print_region_processes_UCD(&fo);
  if(Ippl::myNode()==0) fo.close();

  /*
    Must initialize the solver bevore the grid 
    Grid initialize done inside PoissonSolver!
  */  
  bool printDiag = true;
  GTMGFem<double,3> *mySolver = new GTMGFem<double,3>(grid,printDiag,iter_relax,1.0E-5,dirBC);

  // FIXME: We should create a new communcator containing only the active processors, as
  //        determined by the grid generator.

  grid->getBrickInfo(); // Loose bounding box of grid
  grid->boundbox(globalMin, globalMax); 
 
  msg << "Poisson Solver definition file " << argv[1] << endl;
  msg << "Omega   = "
      << "x=(" << globalMin[0] << " ... " << globalMax[0] << ") "
      << "y=(" << globalMin[1] << " ... " << globalMax[1] << ") "
      << "z=(" << globalMin[2] << " ... " << globalMax[2] << ")" << endl;
  
  // Bounding box of (sub)grid assigned to the local processor
  grid->localbox(localMin, localMax, Ippl::myNode());
  msg2all << "Omega_" << Ippl::myNode() << " = "
	  << "x=(" << localMin[0] << " ... " << localMax[0] << ") "
	  << "y=(" << localMin[1] << " ... " << localMax[1] << ") "
	  << "z=(" << localMin[2] << " ... " << localMax[2] << ")" << endl;
  D3vector my_corner_min(localMin[0], localMin[1], localMin[2]);
  D3vector my_corner_max(localMax[0], localMax[1], localMax[2]);
  
  // Create particle container
  playout_t* layout = new ParticleLayoutFromGrid(grid, geom_domain, my_corner_min, my_corner_max);
  ChargedParticles<playout_t> *bunch = new ChargedParticles<playout_t>(layout, myBC);
  
  // Create totalP particles with random positions uniformy distributed over the grid
  IpplTimings::startTimer(partCreaTimer);
  createDistr(bunch,totalP,globalMin,globalMax,R1,r,charge);
  IpplTimings::stopTimer(partCreaTimer);
 
  bunch->update();
  bunch->writePhaseSpaceSDDS(string("data/")+string(argv[0]));

  msg << "h= " << meshsize << " a= " << a << " b= " << b << " l= " << l
      << " R1= " << R1 << " r= " << r << " Np= " << totalP 
      << " Q= " << sum(bunch->Q)*(4.0*pi*Physics::epsilon_0) << " BC= " << dirBC << endl;    
  
  msg << "Q1= " << sum(bunch->Q)*(4.0*pi*Physics::epsilon_0) << endl;
  
  if (test1)
    mySolver->test1(bunch);  
  else if (test2) {
    mySolver->computeElectricField(bunch,dirBC);
    bunch->writePhaseSpaceDiag(string("data/")+string(argv[0]));
    msg << "Q2= " << sum(bunch->Q)*(4.0*pi*Physics::epsilon_0) << endl;
    mySolver->computeElectricField(bunch,dirBC);
    bunch->writePhaseSpaceDiag(string("data/")+string(argv[0]));
  }
  
  IpplTimings::stopTimer(mainTimer); 
  IpplTimings::print();
  
  delete bunch;  // No "delete layout;" already done in ParticleBase::~ParticleBase()
  delete grid;
  delete geom_domain;
  return 0;
}
/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/
