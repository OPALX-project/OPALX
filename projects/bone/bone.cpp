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

 This test program sets up a simple sine-wave electric field in 3D, 
   creates a population of particles with random q/m values (charge-to-mass
   ratio) and velocities, and then tracks their motions in the static
   electric field using nearest-grid-point interpolation. 

Usage:

 
 $MYMPIRUN -np 2 bone 128 128 128 --processes 2 --commlib mpi

***************************************************************************/

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>

#include "mpi.h"

// dimension of our positions
const unsigned Dim = 3;

// some typedefs 
//typedef ParticleSpatialLayout<double,Dim> playout_t;

typedef UniformCartesian<Dim,double> Mesh_t;
typedef Cell                                       Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t> FieldLayout_t;


#define GUARDCELL 1


enum BC_t {OOO,OOP,PPP};
enum InterPol_t {NGP,CIC};

const double pi = acos(-1.0);
const double qmmax = 1.0;       // maximum value for particle q/m
const double dt = 1.0;          // size of timestep

int main(int argc, char *argv[]){
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);

  Vektor<int,Dim> nr(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

  msg << "Bone  test  grid = " << nr <<endl;
  
  e_dim_tag decomp[Dim];

  unsigned procs[Dim] = {1,4,4};
  NDIndex<Dim> domain;

  for(int i=0; i<Dim; i++)
      domain[i] = domain[i] = Index(nr[i] + 1);

  for (int d=0; d < Dim; ++d)
      decomp[d] =  PARALLEL;
  
  // create mesh and layout objects for this problem domain

  Mesh_t *mesh      = new Mesh_t(domain);
  FieldLayout_t *FL = new FieldLayout_t(*mesh, decomp, procs);

  Field<char,Dim> b;
  b.initialize(*mesh, *FL, GuardCellSizes<Dim>(1));

  msg << *mesh << endl;

  msg << *FL << endl;

  int a1 = b.getLayout().getLocalNDIndex()[0].length();
  int a2 = b.getLayout().getLocalNDIndex()[1].length();
  int a3 = b.getLayout().getLocalNDIndex()[2].length();

  msg2all << b.getLayout().getLocalNDIndex() << " size= " << a1*a2*a3 << endl;


  Ippl::Comm->barrier();
  msg << "Particle test PIC3d: End." << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

