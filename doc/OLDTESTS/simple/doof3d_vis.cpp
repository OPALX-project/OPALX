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
#include <iostream.h>


int main(int argc, char *argv[]) {

  Ippl ippl(argc,argv);
  Inform msg(argv[0]);

  int vnodes     = 16;
  int sizeX      = 50, sizeY   = 50, sizeZ   = 50;
  int iterations = 100;

  // use command-line args if they are given
  if (argc > 1) {
    if (argc != 6) {
      msg << "ERROR: bad command-line options." << endl;
      msg << "Usage: " << argv[0] << " <nx> <ny> <nz> <vnodes> <iters>";
      msg << endl;
      return 1;
    }

    sizeX = atoi(argv[1]);
    sizeY = atoi(argv[2]);
    sizeZ = atoi(argv[3]);
    vnodes = atoi(argv[4]);
    iterations = atoi(argv[5]);
  }

  int centerX    = sizeX/2, centerY = sizeY/2, centerZ = sizeZ/2;
  int a1 = sizeX/10;
  int a2 = sizeY/4;
  int a3 = sizeZ/3;

  const unsigned Dim=3;
  Index I(sizeX);
  Index J(sizeY);
  Index K(sizeZ);
  FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(2));
  Field<double,Dim> B(layout);

  A = 0.0;
  A[centerX][centerY][centerZ] = 512.0*iterations;
  A[centerX+a1][centerY+a2][0] = 512.0*iterations;
  A[centerX-a1][centerY-a2][centerZ+a3] = 128.0*iterations;
  A[sizeX-1][centerY+a2][centerZ] = 512.0*iterations;
  A[centerX+a1][centerY][centerZ-a3] = 512.0*iterations;

  A.connect("A");
  A.updateConnection();
  A.updateConnection();
  A.interact(argc > 1 ? argv[1] : 0);

  double fact = 1.0/27.0;
  for(int iter = 0 ; iter < iterations ; iter++ ) {
    msg << "Computing new A ..." << endl;
    IpplCounter ca("3D stencil in 4 expressions");
    ca.startCounter();
    assign(B[I][J][K], fact*(A[I  ][J  ][K+1] + A[I  ][J  ][K-1] +
		             A[I  ][J+1][K  ] + A[I  ][J-1][K  ] +
		             A[I+1][J  ][K  ] + A[I-1][J  ][K  ]));
    assign(B[I][J][K], B[I][J][K] + fact*(
               A[I+1][J+1][K  ] + A[I+1][J-1][K  ] +
               A[I][J][K] + 
	       A[I-1][J+1][K  ] + A[I-1][J-1][K  ]));
    assign(B[I][J][K], B[I][J][K] + fact*(
               A[I+1][J+1][K+1] + A[I+1][J  ][K+1] + A[I+1][J-1][K+1] +
	       A[I  ][J+1][K+1] +                    A[I  ][J-1][K+1] +
	       A[I-1][J+1][K+1] + A[I-1][J  ][K+1] + A[I-1][J-1][K+1]));
    assign(B[I][J][K], B[I][J][K] + fact*(
               A[I+1][J+1][K-1] + A[I+1][J  ][K-1] + A[I+1][J-1][K-1] +
	       A[I  ][J+1][K-1] +                    A[I  ][J-1][K-1] +
	       A[I-1][J+1][K-1] + A[I-1][J  ][K-1] + A[I-1][J-1][K-1]));
    ca.stopCounter();
    ca.printIt();

    assign(A,B);
    msg << "iter = " << iter << ", sum = " << sum(A) << endl;
    msg << "Updating connection ..." << endl;
    A.updateConnection();
    A.interact();
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: doof3d_vis.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:55 $
 * IPPL_VERSION_ID: $Id: doof3d_vis.cpp,v 1.1.1.1 2000/11/30 21:05:55 adelmann Exp $ 
 ***************************************************************************/
