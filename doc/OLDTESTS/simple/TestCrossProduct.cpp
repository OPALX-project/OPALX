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

// -*- C++ -*-
//-----------------------------------------------------------------------------
// The IPPL Framework - Visit http://people.web.psi.ch/adelmann/ for more details
//
// This program was prepared by the Regents of the University of California at
// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "AppTypes/Vektor.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"

#ifdef IPPL_TEST_CROSS_PRODUCT_DEBUG
#include "Utility/FieldDebug.h"
#endif

#include <math.h>

// set dimensionality and problem size 
const unsigned Dim=3;
int size = 2;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  Vektor<double, Dim > v1(1.0,0.0,0.0);
  Vektor<double, Dim > v2(0.0,2.0,0.0);
  Vektor<double, Dim > v3;
  Vektor<double, Dim > v1Crossv2(0.0,0.0,2.0);

  Index I(size), J(size), K(size);
  FieldLayout<Dim> layout(I,J,K);
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim > > 
    V1(layout), V2(layout), V3(layout), V1CrossV2(layout);
  Field<bool,Dim,UniformCartesian<Dim> > checkCross(layout);

  // test Vektor Cross-product operations:
  // Scalar:
  v3 = cross(v1,v2);
  // Field:
  V1 = v1;
  V2 = v2;
  V1CrossV2 = v1Crossv2;
  V3 = cross(V1,V2);
  
#ifdef IPPL_TEST_CROSS_PRODUCT_DEBUG
  testmsg << " v3 " << endl << v3 << endl << endl;
  testmsg << "((((((((((((((V3))))))))))))))" << endl;
  setInform(testmsg);
  setFormat(1,2,0);
  fp3(V3);
#endif // IPPL_TEST_CROSS_PRODUCT_DEBUG

  // Check results:
  if (v3 != v1Crossv2) {
    passed = false;
    testmsg << "(v3 != v1Crossv2)" << endl;
  }
  checkCross = where(ne(V3, V1CrossV2), true, false);
  bool collectCrossField = true;
  collectCrossField = sum(checkCross);
  if (collectCrossField) {
    passed = false;
    testmsg << "(V3 != V1CrossV2)" << endl;
    testmsg << "collectCrossField = " << collectCrossField << endl;
  }

  // Another test:
  // Scalar:
  double s2i = 1.0/sqrt(2.0);
  v1 = Vektor<double,Dim>(0.0, 1.0, 0.0);
  v2 = Vektor<double,Dim>(s2i, s2i, 0.0);
  v1Crossv2 = Vektor<double,Dim>(0.0, 0.0, -s2i);
  v3 = 3.0*cross(v1,v2);
  // Field:
  V1 = v1;
  V2 = v2;
  V1CrossV2 = v1Crossv2;
  V3 = cross(V1,V2);

  // Check results:
  if (v3 != 3.0*v1Crossv2) {
    passed = false;
    testmsg << "Test 2:" << endl;
    testmsg << "(v3 != 3.0*v1Crossv2)" << endl;
    testmsg << " v3 = " << v3 << endl << endl;
    testmsg << " 3.0*v1Crossv2 = " << 3.0*v1Crossv2 << endl << endl;
  }
  checkCross = where(ne(V3, V1CrossV2), true, false);
  collectCrossField = true;
  collectCrossField = sum(checkCross);
  if (collectCrossField) {
    passed = false;
    testmsg << "Test 2:" << endl;
    testmsg << "(V3 != V1CrossV2)" << endl;
    testmsg << "collectCrossField = " << collectCrossField << endl;
  }

  if (!passed) {
    testmsg << "For reference, true=" << true << " ; false=" << false << endl;
  }

#ifdef IPPL_TEST_CROSS_PRODUCT_CHECK_UNDEFINED_DIMS
  testmsg << "Trying to do 1D and 2D..." << endl;
  testmsg << ( (false) ? "PASSED" : "FAILED" ) << endl;
  Vektor<double,1> va1(1.0);
  Vektor<double,1> va2(0.0);
  Vektor<double,1> va3;
  va3 = cross(va1,va2);
  Vektor<double,2> vb1(1.0,0.0);
  Vektor<double,2> vb2(0.0,2.0);
  Vektor<double,2> vb3;
  vb3 = cross(vb1,vb2);
#endif // IPPL_TEST_CROSS_PRODUCT_CHECK_UNDEFINED_DIMS
    
  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: TestCrossProduct.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:01 $
 * IPPL_VERSION_ID: $Id: TestCrossProduct.cpp,v 1.1.1.1 2000/11/30 21:05:01 adelmann Exp $ 
 ***************************************************************************/
