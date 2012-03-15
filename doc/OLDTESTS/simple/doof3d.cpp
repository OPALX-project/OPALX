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
  TAU_PROFILE("main", "int (int, char**)", TAU_DEFAULT);
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
  Field<double,Dim> B(layout);

  //  A[I][J][K] = IpplRandom;
  A[size/2][size/2][size/2] = 512.0*(iterations + 1);
  /*
  A[size/2+3][size/2+7][0] = 512.0*iterations;
  A[size/2-4][size/2-15][size/2+20] = 128.0*iterations;
  A[size/2-1][size/2+10][size/2] = 512.0*iterations;
  A[size/2+12][size/2][size/2-5] = 512.0*iterations;
  */

  double fact = 1.0/27.0;

  for(int iter = 0 ; iter < iterations ; iter++ ) {
    // Do some scalar code first on A and B
    if (printsum != 0) {
      msg << "Scalar code ..." << endl;
      A.prepareForScalarCode();
      msg << "  First element = " << *(A.begin()) << endl;
      A.finishScalarCode(false);
    }

    msg << "Computing new A ..." << endl;
    B[I][J][K]  = fact * (
		    A[I  ][J  ][K+1] + A[I  ][J  ][K-1] +
	            A[I  ][J+1][K  ] + A[I  ][J-1][K  ] +
	            A[I+1][J  ][K  ] + A[I-1][J  ][K  ]
                  );
    B[I][J][K] += fact * (
                    A[I+1][J+1][K  ] + A[I+1][J-1][K  ] +
                    A[I  ][J  ][K  ] + A[I-1][J+1][K  ] +
                    A[I-1][J-1][K  ]
                  );
    B[I][J][K] += fact * (
                    A[I+1][J+1][K+1] + A[I+1][J  ][K+1] +
                    A[I+1][J-1][K+1] + A[I  ][J+1][K+1] +
                    A[I  ][J-1][K+1] + A[I-1][J+1][K+1] +
                    A[I-1][J  ][K+1] + A[I-1][J-1][K+1]
                  );
    B[I][J][K] += fact * (
                    A[I+1][J+1][K-1] + A[I+1][J  ][K-1] +
                    A[I+1][J-1][K-1] + A[I  ][J+1][K-1] +
                    A[I  ][J-1][K-1] + A[I-1][J+1][K-1] +
                    A[I-1][J  ][K-1] + A[I-1][J-1][K-1]
                  );
    A = B;
    if (printsum != 0)
      msg << "iter = " << iter << ", sum = " << sum(A) << endl;
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: doof3d.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:53 $
 * IPPL_VERSION_ID: $Id: doof3d.cpp,v 1.1.1.1 2000/11/30 21:05:53 adelmann Exp $ 
 ***************************************************************************/
