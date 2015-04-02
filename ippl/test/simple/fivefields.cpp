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
  
  Ippl ippl(argc,argv);
  Inform msg(argv[0]);

  if (argc != 5) {
    msg << "Usage: " << argv[0] << " <size> <vnodes> <iterations> <printsum(1 or 0)>" << endl;
    exit(1);
  }
  int size = atoi(argv[1]);
  int vnodes = atoi(argv[2]);
  int iterations = atoi(argv[3]);
  int printsum = atoi(argv[4]);

  const unsigned Dim=3;
  Index I(size);
  Index J(size);
  Index K(size);
  FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(2));
  Field<double,Dim> B(layout,GuardCellSizes<Dim>(2));
  Field<double,Dim> C(layout,GuardCellSizes<Dim>(2));
  Field<double,Dim> D(layout,GuardCellSizes<Dim>(2));
  Field<double,Dim> E(layout,GuardCellSizes<Dim>(1));

  FieldLayout<2> layout2(I,J,PARALLEL,PARALLEL, vnodes);

  Field<double,2> D2(layout2,GuardCellSizes<2>(2));


  A = 0.0;
  B = 0.0;
  C = 0.0;
  D = 0.0;
  E = 0.0;
  D2 = 0.0;

  A[size/2][size/2][size/2] = 512.0*(iterations + 1);
  B[size/2][size/2][size/2] = 512.0*(iterations + 2);
  C[size/2][size/2][size/2] = 512.0*(iterations + 3);
  D[size/2][size/2][size/2] = 512.0*(iterations + 4);

  for(int iter = 0 ; iter < iterations ; iter++ ) {
    if (printsum != 0 || iter % 10 == 0)
      msg << "Computing new values at iteration " << iter << " ..." << endl;

    E[I][J][K]  = A[I][J][K+2] + B[I][J+2][K] + C[I-2][J][K] + D[I][J][K];

    A = E;

    C = A + B;

    B = E;

    if (printsum != 0)
      {
	msg << "  iter = " << iter << ", sum(A) = " << sum(A) << endl;
	msg << "  iter = " << iter << ", sum(C) = " << sum(C) << endl;
	msg << "  iter = " << iter << ", sum(E) = " << sum(E) << endl;
      }
  }

  msg << A[I][J][1] << endl;
  

  return 0;
}

/***************************************************************************
 * $RCSfile: fivefields.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: fivefields.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
