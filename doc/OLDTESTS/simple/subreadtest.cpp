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

// A simple test code to test writing out data to a DiscField,
// then reading in only a portion of it into a second field with
//   1. The same domain as the first, but different layout
//   2. A different domain as well as a different layout

#include "Ippl.h"

int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);
  Inform msg(argv[0]);

  if (argc != 5) {
    msg << "Usage: " << argv[0] << " <X> <Y> <vnodesout> <vnodesin>" << endl;
    exit(1);
  }

  int sizeX   = atoi(argv[1]);
  int sizeY   = atoi(argv[2]);
  int vnodes  = atoi(argv[3]);
  int vnodes2 = atoi(argv[4]);
  int fields  = 4;

  if (vnodes < 1 || fields < 1 || sizeX < 1 || sizeY < 1) {
    msg << "Error in input values." << endl;
    exit(1);
  }
    
  const unsigned Dim=2;
  typedef int T;
  Index I(1, sizeX), I2(1 + sizeX/3, 1 + 2*sizeX/3), I3(0,2*sizeX + 1);
  Index J(1, sizeY), J2(1 + sizeY/3, 1 + 2*sizeY/3), J3(0,2*sizeY + 1);
  NDIndex<Dim> domain(I, J);
  NDIndex<Dim> subdom(I2, J2);
  NDIndex<Dim> bigdom(I3, J3);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL, vnodes);
  FieldLayout<Dim> layout2(I,J,PARALLEL,PARALLEL, vnodes2);
  FieldLayout<Dim> layout3(I3,J3,PARALLEL,PARALLEL, vnodes2);
  Field<T,Dim> A(layout);
  Field<T,Dim> A2(layout);
  Field<T,Dim> B(layout2);
  Field<T,Dim> B2(layout);
  Field<T,Dim> C(layout3);
  Field<T,Dim> C2(layout3);

  // initialize data
  msg << "Initializing A ..." << endl;
  A[I][J] = I + 10*(J-1);
  //FieldPrint<T,Dim> fp(A);
  //FieldPrint<T,Dim> fpb(B);
  //msg << "Initial A:" << endl;
  //fp.print(domain);
  //msg << "A subset:" << endl;
  //fp.print(subdom);

  // write A to disk N times
  msg << "====================== writing ====================" << endl;
  DiscField<Dim> dcf("ackfiledata", "ackfiledata.config", fields, "data");
  msg << "Writing field " << 0 << " ..." << endl;
  dcf.write(A, 0);
  A2 = A;
  for (int i=1; i < fields; ++i) {
    msg << "Incrementing field ..." << endl;
    A2 += 100;

    msg << "Writing field " << i << " ..." << endl;
    dcf.write(A2, i);
  }

  // read in the first and last field
  msg << "====================== reading ====================" << endl;
  DiscField<Dim> dcf2("ackfiledata", "ackfiledata.config");
  msg << "Reading in field 0 ..." << endl;
  dcf2.read(B, 0);
  //msg << "Field just read: " << endl;
  //fpb.print(domain);
  B2 = B;
  B2 -= A;
  NDIndex<Dim> mloc;
  msg << "Checking field just read: min and max of diff should be zero." << endl;
  msg << "     min(B-A) = " << min(B2, mloc) << endl;
  msg << "  minloc(B-A) = " << mloc << endl;
  msg << "     max(B-A) = " << max(B2, mloc) << endl;
  msg << "  maxloc(B-A) = " << mloc << endl;
  msg << "       sum(A) = " << sum(A) << endl;
  msg << "       sum(B) = " << sum(B) << endl;

  A2 = 0;
  A2[I2][J2] = I2 + 10*(J2-1);
  msg << "Reading in field 0 again, with subdomain " << subdom <<" ..."<<endl;
  B = 0;
  dcf2.read(B, subdom, 0);
  //msg << "Field just read: " << endl;
  //fpb.print(domain);
  B2 = B;
  B2 -= A2;
  msg << "Checking field just read: min and max of diff should be zero." << endl;
  msg << "     min(B-A2) = " << min(B2, mloc) << endl;
  msg << "  minloc(B-A2) = " << mloc << endl;
  msg << "     max(B-A2) = " << max(B2, mloc) << endl;
  msg << "  maxloc(B-A2) = " << mloc << endl;
  msg << "       sum(A2) = " << sum(A2) << endl;
  msg << "        sum(B) = " << sum(B) << endl;

  msg << "Reading in field 0 again into big array, with subdomain " << subdom <<" ..."<<endl;
  C = 0;
  C2 = 0;
  C2[I2][J2] = I2 + 10*(J2-1);
  dcf2.read(C, subdom, 0);
  msg << "Checking field just read: min and max of diff should be zero." << endl;
  msg << "     min(C2-C) = " << min(C2-C, mloc) << endl;
  msg << "  minloc(C2-C) = " << mloc << endl;
  msg << "     max(C2-C) = " << max(C2-C, mloc) << endl;
  msg << "  maxloc(C2-C) = " << mloc << endl;
  msg << "       sum(C2) = " << sum(C2) << endl;
  msg << "        sum(C) = " << sum(C) << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: subreadtest.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:28 $
 * IPPL_VERSION_ID: $Id: subreadtest.cpp,v 1.1.1.1 2000/11/30 21:06:28 adelmann Exp $ 
 ***************************************************************************/
