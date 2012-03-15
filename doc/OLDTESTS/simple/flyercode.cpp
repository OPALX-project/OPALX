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
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// -*- C++ -*-
// flyercode.cpp , Tim Williams 10/10/1997

// This is the 2D diffusion example from 1997 IPPL flyer (SC97 handout).

#include "Ippl.h"
int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv); // Initialize library, parse command-line arguments
  Inform d2dout("d2d",0);  // Like cout except only virtual processor 0 prints

  // Uniform 2D cartesian mesh,129x129 vertices, default spacing = (1.0,1.0):
  unsigned N = 129; Index Iverts(N); Index Jverts(N);
  typedef UniformCartesian<2> M;
  M mesh(Iverts,Jverts);
 
  // Boundary conditions--zero on all faces:
  Field<double,2>::bcond_container bc;
  for (int f=0; f<4; f++) bc[f] = new ConstantFace<double,2>(f,0.0);

  // Construct a 2D Field of double's representing scalar fluid field U(x,y).
  // Cell-centering on mesh, meaning 128x128 Field elements:
  CenteredFieldLayout<2,M,Cell> layout(mesh,PARALLEL,PARALLEL);  
  Field<double,2,M,Cell> U(layout,GuardCellSizes<2>(1),bc);  

  // Initial conditions, zero except one-cell spike in center of box:
  U = 0.0;  U[N/2][N/2] = 1000.0;
  double dt = 0.1; // Timestep

  // Compute global sum of values from field (diagnostic):
  double sumU = sum(U);
  d2dout << "sum at t = 0 is " << sumU << endl;

  Index I(128), J(128);  // Index objects of whole-system size (all cells)

  // Run the diffusion stencil through 15 timesteps:
  for (int itime=0; itime < 15; itime++)
    U[I][J] += dt*(U[I+1][J] + U[I-1][J] - 4*U[I][J] + U[I][J+1] + U[I][J-1]);

  // Recompute global sum of values from field (diagnostic):
  d2dout << "sum at t = " << 15*dt << " is " << sum(U) << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: flyercode.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:07 $
 * IPPL_VERSION_ID: $Id: flyercode.cpp,v 1.1.1.1 2000/11/30 21:06:07 adelmann Exp $ 
 ***************************************************************************/
