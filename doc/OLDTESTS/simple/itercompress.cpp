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

/****************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by the Regents of the University of
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

//
// NOTE: This test will print FAILED now for some tests, since it
// tests a capability of IPPL that has been removed.  THe capability
// was "automatic compression via iterators" for Barefield.  This was
// removed because it caused poor performance in many situations.
// DO NOT use this test code for regression tests, etc.
//

#include "Ippl.h"

int main(int argc, char* argv[]) {

  Ippl ippl(argc,argv);

  Index I(10), J(10);
  FieldLayout<2> FL(I,J,PARALLEL,PARALLEL,4);
  Field<double,2> A(FL), B(FL);
  if (Ippl::noFieldCompression) {
    if (A.CompressedFraction() == 0.0)
      cout << "PASSED: Constructed without compression." << endl;
    else
      cout << "FAILED: Constructed without compression." << endl;
    if (B.CompressedFraction() == 0.0)
      cout << "PASSED: Constructed without compression." << endl;
    else
      cout << "FAILED: Constructed without compression." << endl;
  }
  else {
    if (A.CompressedFraction() == 1.0)
      cout << "PASSED: Constructed with compression." << endl;
    else
      cout << "FAILED: Constructed with compression." << endl;
    if (B.CompressedFraction() == 1.0)
      cout << "PASSED: Constructed with compression." << endl;
    else
      cout << "FAILED: Constructed with compression." << endl;
  }

  A.Uncompress();
  {
    Field<double,2>::iterator Ait = A.begin(), Bit = B.begin();
    for (int i = 0; i < 5; ++i, ++Ait, ++Bit) {
      *Ait = 3.0;
      *Bit = 3.0;
    }
  }
  A.Compress();

  cout << "A = " << A << endl;
  cout << "B = " << B << endl;

  double sumA = sum(A);
  double sumB = sum(B);

  if (sumA == 15)
    cout << "PASSED: Iterator assignment" << endl;
  else 
    cout << "FAILED: Iterator assignment" << endl;
  if (sumB == 15)
    cout << "PASSED: Automatic uncompression." << endl;
  else
    cout << "FAILED: Automatic uncompression." << endl;

  if (Ippl::noFieldCompression) {
    if (A.CompressedFraction() == 0.0)
      cout << "PASSED: No compression." << endl;
    else
      cout << "FAILED: No compression." << endl;
    if (B.CompressedFraction() == 0.0)
      cout << "PASSED: No compression." << endl;
    else
      cout << "FAILED: No compression." << endl;
  }
  else {
    if (A.CompressedFraction() == 0.75)
      cout << "PASSED: Manual compression." << endl;
    else
      cout << "FAILED: Manual compression." << endl;
    if (B.CompressedFraction() == 0.75)
      cout << "PASSED: Automatic compression." << endl;
    else
      cout << "FAILED: Automatic compression." << endl;
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: itercompress.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:10 $
 * IPPL_VERSION_ID: $Id: itercompress.cpp,v 1.1.1.1 2000/11/30 21:06:10 adelmann Exp $ 
 ***************************************************************************/
