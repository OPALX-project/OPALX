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

#include "Ippl.h"

int main(int argc, char *argv[]) {
  // initialize the IPPL environment
  Ippl ippl(argc,argv);

  // define constants which describe how the Field's are partitioned
  // and sized
  const unsigned Dim = 2;
  int vnodes         = 16;
  int sizeX          = 50, sizeY   = 50;
  int centerX        = 25, centerY = 25;
  int iterations     = 200;
  double fact        = 1.0/9.0;

  // I and J are Index objects, which represent the size of the domain
  Index I(sizeX);
  Index J(sizeY);

  // layout maintains the parallel data distribution information; it
  // breaks the domain up into 'vnodes' blocks
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL, vnodes);

  // create a Field using this parallel layout, with one layer of
  // guard cells around each vnode block
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));

  // initialize the Field with a few spikes
  A = 0.0;
  A[centerX][centerY] = 200.0*iterations;
  A[centerX/2][centerY/3] = 100.0*iterations;
  A[3*centerX/2][centerY] = 400.0*iterations;
  
  // make the Field A available for visualization, if desired
  A.connect("A");
  A.updateConnection();
  A.updateConnection();
  A.interact(argc > 1 ? argv[1] : 0);

  // for several iterations, compute a new value of A using a
  // 9-point stencil
  for(int iter = 0 ; iter < iterations ; iter++) {
    A[I][J] = fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] + 
		    A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] + 
		    A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1]);

    // update the visualization display, if it is being used
    A.updateConnection();
    A.interact();
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: doof2d_vis.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: doof2d_vis.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/

