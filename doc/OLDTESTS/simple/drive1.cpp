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
  Inform testmsg(argv[0]);

  TAU_PROFILE("main","",TAU_DEFAULT);

  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout);
  Field<double,Dim> B(layout);
  Field<double,Dim> C(layout);
  Field<double,Dim>::iterator p;
  double sumval, testVal;
  double a, b, c;
  double eps = 1.0e-07;
  int testNum = 1;
  Index Itest; // a test Index

  // Flags to keep track of pass/fail results:
  bool passed = true;           // True at end means passed all tests
  bool passedSingleTest = true; // For individual tests

  // -----------------------------------------------------
  // testing accumulation: Field += Field

  // test 1
  A = 0.0;
  B = 1.0;
  a = 0.0;
  b = 1.0;
  A += B;
  a += b;
  passedSingleTest = true;

  for (p=A.begin() ; p!=A.end(); ++p)
    {
      if (fabs(*p - a) > eps)
	{
	  passedSingleTest = false; passed = false;
	  p = A.end();
	}
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 2
  A = 0.0;
  B = 1.0;
  a = 0.0;
  b = 1.0;
  A -= B;
  a -= b;
  passedSingleTest = true;

  for (p=A.begin() ; p!=A.end(); ++p) 
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 3
  A = 1.0;
  B = 2.0;
  a = 1.0;
  b = 2.0;
  A *= B;
  a *= b;
  passedSingleTest = true;

  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      testmsg << "test " <<  testNum << " failed <--" << endl;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 4
  A = 1.0;
  B = 5.0;
  a = 1.0;
  b = 5.0;
  A /= B;
  a /= b;
  passedSingleTest = true;

  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: Field += T

  // test 5
  A = 0.0;
  a = 0.0;
  A += 1.0;
  a += 1.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (*p != a)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 6
  A = 0.0;
  a = 0.0;
  A -= 1.0;
  a -= 1.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (*p != a)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 7
  A = 1.0;
  a = 1.0;
  A *= 2.0;
  a *= 2.0;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (*p != a)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 8
  A = 1.0;
  a = 1.0;
  A /= 5.0;
  a /= 5.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: Field += ExprElem

  // test 9
  A = 0.0;
  B = 1.0;
  C = 2.0;
  a = 0.0;
  b = 1.0;
  c = 2.0;
  A += B + C;
  a += b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 10
  A = 0.0;
  B = 1.0;
  C = 2.0;
  a = 0.0;
  b = 1.0;
  c = 2.0;
  A -= B + C;
  a -= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 11
  A = 1.0;
  B = 2.0;
  C = 3.0;
  a = 1.0;
  b = 2.0;
  c = 3.0;
  A *= B + C;
  a *= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 12
  A = 1.0;
  B = 2.0;
  C = 3.0;
  a = 1.0;
  b = 2.0;
  c = 3.0;
  A /= B + C;
  a /= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: IndexingField += IndexingField

  // test 13
  A = 0.0;
  B = 1.0;
  a = 0.0;
  b = 1.0;
  A[I][J] += B[I][J];
  a += b;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 14
  A = 0.0;
  B = 1.0;
  a = 0.0;
  b = 1.0;
  A[I][J] -= B[I][J];
  a -= b;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 15
  A = 1.0;
  B = 2.0;
  a = 1.0;
  b = 2.0;
  A[I][J] *= B[I][J];
  a *= b;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 16
  A = 1.0;
  B = 5.0;
  a = 1.0;
  b = 5.0;
  A[I][J] /= B[I][J];
  a /= b;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: IndexingField += T

  // test 17
  A = 0.0;
  a = 0.0;
  A[I][J] += 1.0;
  a += 1.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 18
  A = 0.0;
  a = 0.0;
  A[I][J] -= 1.0;
  a -= 1.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 19
  A = 1.0;
  a = 1.0;
  A[I][J] *= 2.0;
  a *= 2.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 20
  A = 1.0;
  a = 1.0;
  A[I][J] /= 5.0;
  a /= 5.0;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: IndexingField += Index

  // test 21
  A = 0.0;
  a = 0.0;
  A[I][J] += I;
  sumval = sum(A);
  passedSingleTest = true;
  if( fabs(sumval - 50.0) >= eps ) {
    passedSingleTest = false; passed = false;;
  }
  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 22
  A = 0.0;
  A[I][J] -= I;
  sumval = sum(A);
  passedSingleTest = true;
  if( fabs(sumval + 50.0) >= eps ) {
    passedSingleTest = false; passed = false;;
  }
  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 23
  A = 1.0;
  A[I][J] *= I;
  sumval = sum(A);
  passedSingleTest = true;
  if( fabs(sumval - 50.0) >= eps ) {
    passedSingleTest = false; passed = false;;
  }
  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // test 24
  A = 1.0;
  Itest = I + 1;
  A[I][J] /= Itest;
  testVal = 5.0*(1 + (1.0/2.0) + (1.0/3.0) + (1.0/4.0) + (1.0/5.0));
  sumval = sum(A);
  passedSingleTest = true;
  if( fabs(sumval - testVal) >= eps ) {
    passedSingleTest = false; passed = false;;
  }
  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;

  // -----------------------------------------------------
  // testing accumulation: IndexingField += ExprElem
  // (with IndexingFields on rhs

  // test 25
  A = 0.0;
  B = 1.0;
  C = 2.0;
  a = 0.0;
  b = 1.0;
  c = 2.0;
  A[I][J] += B[I][J] + C[I][J];
  a += b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 26
  A = 0.0;
  B = 1.0;
  C = 2.0;
  a = 0.0;
  b = 1.0;
  c = 2.0;
  A[I][J] -= B[I][J] + C[I][J];
  a -= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 27
  A = 1.0;
  B = 2.0;
  C = 3.0;
  a = 1.0;
  b = 2.0;
  c = 3.0;
  A[I][J] *= B[I][J] + C[I][J];
  a *= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // test 28
  A = 1.0;
  B = 2.0;
  C = 3.0;
  a = 1.0;
  b = 2.0;
  c = 3.0;
  A[I][J] /= B[I][J] + C[I][J];
  a /= b + c;
  passedSingleTest = true;
  for (p=A.begin() ; p!=A.end(); ++p)
    if (fabs(*p - a) > eps)
    {
      passedSingleTest = false; passed = false;;
      p = A.end();
    }

  if (!passedSingleTest) testmsg << "test "<<  testNum << " failed" << endl;
  testNum++;


  // -----------------------------------------------------
  // testing accumulation: IndexingField += ExprElem
  // (with Fields on rhs

//  // test 29
//  A = 0.0;
//  B = 1.0;
//  C = 2.0;
//  A[I][J] += B + C;
//  testmsg << A << endl;
//  sum = 0.0;
//  for (p=A.begin() ; p!=A.end(); ++p) sum += *p;
//  testmsg << sum << endl;
//  if( fabs(sum - 75.0) < eps ) {
//    testmsg << "test "<<  testNum << " passed" << endl;
//  } else {
//    testmsg << "test "<<  testNum << " failed  <--" << endl;
//  }
//  testNum++;
//
//  // test 30
//  A = 0.0;
//  B = 1.0;
//  C = 2.0;
//  A[I][J] -= B + C;
//  sum = 0.0;
//  for (p=A.begin() ; p!=A.end(); ++p) sum += *p;
//  if( fabs(sum + 75.0) < eps ) {
//    testmsg << "test "<<  testNum << " passed" << endl;
//  } else {
//    testmsg << "test "<<  testNum << " failed  <--" << endl;
//  }
//  testNum++;
//
//  // test 31
//  A = 1.0;
//  B = 2.0;
//  C = 3.0;
//  A[I][J] *= B + C;
//  sum = 0.0;
//  for (p=A.begin() ; p!=A.end(); ++p) sum += *p;
//  if( fabs(sum - 125.0) < eps ) {
//    testmsg << "test "<<  testNum << " passed" << endl;
//  } else {
//    testmsg << "test "<<  testNum << " failed  <--" << endl;
//  }
//  testNum++;
//
//  // test 32
//  A = 1.0;
//  B = 2.0;
//  C = 3.0;
//  A[I][J] /= B + C;
//  sum = 0.0;
//  for (p=A.begin() ; p!=A.end(); ++p) sum += *p;
//  if( fabs(sum - 5.0) < eps ) {
//    testmsg << "test "<<  testNum << " passed" << endl;
//  } else {
//    testmsg << "test "<<  testNum << " failed  <--" << endl;
//  }
//  testNum++;

  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: drive1.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:57 $
 * IPPL_VERSION_ID: $Id: drive1.cpp,v 1.1.1.1 2000/11/30 21:05:57 adelmann Exp $ 
 ***************************************************************************/
