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
  int sizeX      = 20, sizeY   = 20, sizeZ   = 20;
  int centerX    = 10, centerY = 10, centerZ = 10;

  int iterations = 10;

  const unsigned Dim=3;
  Index I(sizeX);
  Index J(sizeY);
  Index K(sizeZ);
  FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));
  Field<double,Dim> B(layout);
  DataConnect *conn = DataConnectCreator::create("Ack.config");
  A.connect("Ack");

  A = 0.0;
  A[centerX][centerY][centerZ] = 512.0*iterations;
  A[centerX+3][centerY+7][0] = 512.0*iterations;
  A[centerX-4][centerY-5][centerZ+8] = 128.0*iterations;
  A[sizeX-1][centerY+6][centerZ] = 512.0*iterations;
  A[centerX+2][centerY][centerZ-5] = 512.0*iterations;

  double fact = 1.0/27.0;

  for(int iter = 0 ; iter < iterations ; iter++ ) {
    msg << "Computing new A ..." << endl;
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
    assign(A,B);
    msg << "iter = " << iter << ", sum = " << sum(A) << endl;
    msg << "Updating connection ..." << endl;
    A.updateConnection();
    msg << "iter = " << iter << ", sum = " << sum(A) << endl;
    A.interact();
  }

  delete conn;
  return 0;
}

/***************************************************************************
 * $RCSfile: doof3d_file.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:54 $
 * IPPL_VERSION_ID: $Id: doof3d_file.cpp,v 1.1.1.1 2000/11/30 21:05:54 adelmann Exp $ 
 ***************************************************************************/
