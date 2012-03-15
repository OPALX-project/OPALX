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

// test program to demonstrate use of the SOffset class
#include "Ippl.h"


void report(const char *str, bool result) {
  Inform msg("Results");
  msg << "Test " << str << ": ";
  msg << (result ? "PASSED" : "FAILED") << endl;
}


int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);

  const unsigned Dim=2;

  // testing SOffset creation
  SOffset<Dim> A(1, 2);
  SOffset<Dim> B(1, 1);
  report("SOffset create",
	 A[0] == 1 && A[1] == 2 && B[0] == 1 && B[1] == 1);

  // testing SOffset [] operator
  SOffset<Dim> soLeft;
  SOffset<Dim> soRight;
  for (int d=0; d < Dim; d++) {
    soLeft[d] = 0;
    soRight[d] = (d == 1 ? 1 : 0);
  }
  report ("SOffset bracket create",
	  soLeft[0] == 0 && soLeft[1] == 0 &&
	  soRight[0] == 0 && soRight[1] == 1);

  // testing SOffset addition, subtraction, and copy constructor
  SOffset<Dim> C(A + B);
  SOffset<Dim> D(A - B);
  report("+", C[0] == 2 && C[1] == 3);
  report("-", D[0] == 0 && D[1] == 1);

  // testing SOffset +=, -= operators
  D += C;
  report("+=", D[0] == 2 && D[1] == 4);
  C -= D;
  report("-=", C[0] == 0 && C[1] == -1);

  // testing SOffset comparisons
  report("B  < A", B < A);
  report("B <= A", B <= A);
  report("D  > C", D > C);
  report("D >= C", D >= C);
  report("A != B", A != B);
  report("B == B", B == B);

  // testing containment check
  NDIndex<Dim> N1(Index(0,1), Index(0,1));
  report("A not inside 0..1,0..1 NDIndex", ! A.inside(N1));
  report("B inside 0..1,0..1 NDIndex", B.inside(N1));

  return 0;
}

/***************************************************************************
 * $RCSfile: soffset.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:26 $
 * IPPL_VERSION_ID: $Id: soffset.cpp,v 1.1.1.1 2000/11/30 21:06:26 adelmann Exp $ 
 ***************************************************************************/
