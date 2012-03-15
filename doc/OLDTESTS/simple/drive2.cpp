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

#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0],INFORM_ALL_NODES);
  bool passed = true; // Pass/fail test

  TAU_PROFILE("main","",TAU_DEFAULT);

  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(domain);
  Field<double,Dim> A(layout);
  Field<double,Dim> B(layout);
  Field<double,Dim> C(layout);
  Field<double,Dim> D(layout);
  Field<double,Dim> E(layout);
  Field<double,Dim>::iterator p;
  double localSum;
  double globalSum;
  double eps = 1.0e-07;
  int testNum = 1;
  Index Itest, Jtest; // test Indices

// -------------------------------------------------------
// testing Binary Operations: Field + Field 

  // test 1
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = B + C;
#else
  A << B + C;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  // This gets the global sum to every PE's copy when multiprocessing:
  globalSum = localSum;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: ExprElem + Field 

  // test 2
  A = 0.0;
  B = 1.0;
  C = 2.0;
  D = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = (B + C) + D;
#else
  A << (B + C) + D;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 150.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;


// -------------------------------------------------------
// testing Binary Operations: Field + ExprElem 

  // test 3
  A = 0.0;
  B = 1.0;
  C = 2.0;
  D = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = B + (C + D);
#else
  A << B + (C + D);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 150.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;


// -------------------------------------------------------
// testing Binary Operations: ExprElem + ExprElem 

  // test 4
  A = 0.0;
  B = 1.0;
  C = 2.0;
  D = 3.0;
  E = 4.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = (B + C) + (D + E);
#else
  A << (B + C) + (D + E);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 250.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;


// -------------------------------------------------------
// testing Binary Operations: Field + Scalar 

  // test 5
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = B + 1.0;
#else
  A << B + 1.0;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 50.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Scalar + Field 

  // test 6
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = 1.0 + B;
#else
  A << 1.0 + B;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 50.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: ExprElem + Scalar 

  // test 7
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = (B + C) + 1.0;
#else
  A << (B + C) + 1.0;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 100.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Scalar + ExprElem 

  // test 8
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = 1.0 + (B + C);
#else
  A << 1.0 + (B + C);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 100.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Index + Index 

  // test 9
  A = 0.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = I + J;
#else
  A[I][J] << I + J;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 100.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Index + Scalar 

  // test 10
  A = 0.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = I + 1.0;
#else
  A[I][J] << I + 1.0;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Scalar + Index 

  // test 11
  A = 0.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = 1.0 + I;
#else
  A[I][J] << 1.0 + I;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: ExprElem + Index 

  // test 12
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = (B[I][J] + C[I][J]) + I;
#else
  A[I][J] << (B[I][J] + C[I][J]) + I;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 125.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Index + ExprElem 

  // test 13
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = I + (B[I][J] + C[I][J]);
#else
  A[I][J] << I + (B[I][J] + C[I][J]);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 125.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: IndexingField + Index 

  // test 14
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = B[I][J] + I;
#else
  A[I][J] << B[I][J] + I;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: Index + IndexingField 

  // test 15
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = I + B[I][J];
#else
  A[I][J] << I + B[I][J];
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: ExprElem + IndexingField 

  // test 16
  A = 0.0;
  B = 1.0;
  C = 2.0;
  D = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = (B[I][J] + C[I][J]) + D[I][J];
#else
  A[I][J] << (B[I][J] + C[I][J]) + D[I][J];
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 150.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: IndexingField + ExprElem 

  // test 17
  A = 0.0;
  B = 1.0;
  C = 2.0;
  D = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = B[I][J] + (C[I][J] + D[I][J]);
#else
  A[I][J] << B[I][J] + (C[I][J] + D[I][J]);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 150.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: IndexingField + IndexingField 

  // test 18
  A = 0.0;
  B = 1.0;
  C = 2.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = B[I][J] + C[I][J];
#else
  A[I][J] << B[I][J] + C[I][J];
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 75.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: IndexingField + T 
//
  // test 19
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = B[I][J] + 1.0;
#else
  A[I][J] << B[I][J] + 1.0;
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 50.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

// -------------------------------------------------------
// testing Binary Operations: T + IndexingField 

  // test 20
  A = 0.0;
  B = 1.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = 1.0 + B[I][J];
#else
  A[I][J] << 1.0 + B[I][J];
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - 50.0) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

  // -------------------------------------------------------
  // testing Binary Operations: pow(Field,Field)
  // test 21
  A = 0.0;
  B = 2.0;
  C = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = pow(B , C);
#else
  A << pow(B , C);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - (25*8)) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

  // -------------------------------------------------------
  // testing Binary Operations: pow(Field,scalar)
  // test 22
  A = 0.0;
  B = 2.0;
  C = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = pow(B , 4.0);
#else
  A << pow(B , 4.0);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - (25*16)) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;


  // -------------------------------------------------------
  // testing Binary Operations: pow(scalar,Field)
  // test 23
  A = 0.0;
  B = 2.0;
  C = 3.0;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A = pow(4 , C);
#else
  A << pow(4 , C);
#endif
  localSum = 0.0;
  for (p=A.begin() ; p!=A.end(); ++p) localSum += *p;
  globalSum = localSum;
  reduce(globalSum, globalSum, OpAddAssign());
  if( fabs(globalSum - (25*64)) < eps ) {
    //    testmsg << "test "<< testNum << " passed" << endl;
  } else {
    passed = false;
    testmsg << "test "<< testNum << " failed  <--" << endl;
  }
  testNum++;

  Inform testmsg0(argv[0],0);
  testmsg0 << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;

}
/***************************************************************************
 * $RCSfile: drive2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:00 $
 * IPPL_VERSION_ID: $Id: drive2.cpp,v 1.1.1.1 2000/11/30 21:06:00 adelmann Exp $ 
 ***************************************************************************/
