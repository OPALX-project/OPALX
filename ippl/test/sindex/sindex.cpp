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
#include "Index/SIndex.h"

/***************************************************************************
  A simple program to test the capabilities of the SIndex and SOffset classes.
 ***************************************************************************/

int main(int argc, char *argv[]) {
  
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  const unsigned Dim=2;
  Index I(25);
  Index J(50);
  Index I2(2);
  Index J2(5);
  NDIndex<Dim> NDX(I,J);
  int IP[Dim];
  IP[0] = -1;
  IP[1] = -1;

  SOffset<Dim> so1(1,1);
  SOffset<Dim> so2(IP);
  testmsg << "Created SOffset so1 = " << so1 << endl;
  testmsg << "Created SOffset so2 = " << so2 << endl;
  testmsg << "Adding IP to so1  = " << so1 + IP << endl;
  testmsg << "Adding so1 to IP = " << IP + so1 << endl;
  testmsg << "Adding so2 to so1 = " << so1 + so2 << endl;
  testmsg << "Adding so1 to NDX " << NDX << " = " << NDX + so1 << endl;
  testmsg << "Adding NDX " << NDX << " to so1 = " << so1 + NDX << endl;
  testmsg << "Multipplying NDX " << NDX << " by so2 = " << NDX * so2 << endl;
  testmsg << "Multipplying so2 by NDX " << NDX << " = " << so2 * NDX  << endl;
  so1 -= IP;
  so2 += IP;
  testmsg << "Accumulated IP from so1 = " << so1 << endl;
  testmsg << "Accumulated IP into so2 = " << so2 << endl;

  FieldLayout<Dim> layout(I, J, PARALLEL, PARALLEL, 2*Ippl::getNodes());
  Field<double,Dim> A(layout);
  Field<bool,Dim> B(layout);

  SIndex<Dim> s1(layout);
  SIndex<Dim> s2 = s1(1,-1);
  SIndex<Dim> x3 = s1(IP);

  testmsg << "Created s1 = " << s1 << endl;
  testmsg << "Created s2 = " << s2 << endl;
  testmsg << "Created x3 = " << x3 << endl;

  s1.addIndex(NDIndex<Dim>(Index(2), Index(3)));
  s2.addIndex(SOffset<Dim>(0,0));
  s2.addIndex(NDIndex<Dim>(Index(20,23), Index(45,46)));

  testmsg << "Added new points, s1 = "; s1.printDebug(testmsg);
  testmsg << "Added new points, s2 = "; s2.printDebug(testmsg);

  SIndex<Dim> s3(s1);
  s3 = NDIndex<Dim>(Index(1,5), Index(2,4));

  testmsg << "Created s3 = " << s3 << endl;

  s3 &= s1;

  testmsg << "Intersection of s3 and s1 = " << s3 << endl;

  // now, test assigment of a Field expression to an SIndex
  A[I][J] = I + J;
  B = lt(A,10);
  s3 = lt(A,10);

  // do a union with a slightly different condition
  s3 |= (gt(A,15) && lt(A,20));
  testmsg << "union of s3 and where 15 < A < 20 = " << s3 << endl;

  // do an intersection with an overlapping condition
  s3 &= (gt(A,8) && lt(A,12));
  testmsg << "intersection of s3 and where 8 < A < 12 = " << s3 << endl;

  // do an indexed assignment to s3
  s3[I2][J2] = (lt(A[I2 + 1][J2+1], 5) && gt(A[I2 + 2][J2 + 2], 0));
  testmsg << "s3[I2][J2] = expr ==> s3 = " << s3 << endl;

  // now use s3 in a Field expression
  testmsg << "Originally, A = I + J ... after A[s3] = A[s3(1,1)]: A = ";
  A[s3] = A[s3(1,1)] + 10;
  testmsg << A << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: sindex.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
 * IPPL_VERSION_ID: $Id: sindex.cpp,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
 ***************************************************************************/
