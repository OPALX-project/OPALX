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


/***************************************************************************
  Test out filling the guard cells for a 2D field, printing out their
  values in the process.
 ***************************************************************************/

template<class T, unsigned Dim>
void printfield(Inform& testmsg, BareField<T,Dim>& A) {
  BareField<int,Dim>::iterator_if lf_i = A.begin_if();
  BareField<int,Dim>::iterator_if lf_e = A.end_if();
  for ( ; lf_i != lf_e; ++lf_i) {
    LField<int,Dim> &lf = *((*lf_i).second);
    const NDIndex<Dim>& lo = lf.getAllocated();
    LField<int,Dim>::iterator data = lf.begin(lo);
    testmsg << "LField on domain " << lf.getOwned() << ":" << endl;
    for (unsigned int j=0; j < lo[1].length(); ++j) {
      for (unsigned int i=0; i < lo[0].length(); ++i) {
	testmsg << "  " << data.offset(i,j);
      }
      testmsg << endl;
    }
  }
}


int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  const unsigned Dim=2;
  const unsigned GCSize=1;
  Index I(5), J(5);
  FieldLayout<Dim> LC(I, J, SERIAL, PARALLEL, 4);
  GuardCellSizes<Dim> gc(GCSize);
  BConds<int,Dim> bc;
  bc[0] = new ZeroFace<int,Dim>(0);
  bc[1] = new ZeroFace<int,Dim>(1);
  bc[2] = new ZeroFace<int,Dim>(2);
  bc[3] = new ZeroFace<int,Dim>(3);
  /*
  bc[0] = new PeriodicFace<int,Dim>(0);
  bc[1] = new PeriodicFace<int,Dim>(1);
  bc[2] = new PeriodicFace<int,Dim>(2);
  bc[3] = new PeriodicFace<int,Dim>(3);
  */

  Field<int,Dim,UniformCartesian<Dim>,Cell> A(LC, gc, bc);
  Field<int,Dim,UniformCartesian<Dim>,Cell> B(LC, gc);
  DebugFieldPrint<int,Dim> dfp(true);

  // test 1: assignment of constant value != zero
  A[I][J] = 3;
  B[I][J] = 3;
  testmsg << "Results of A[I][J] = 3:" << endl;
  printfield(testmsg, A);
  testmsg << "Results of B[I][J] = 3:" << endl;
  printfield(testmsg, B);
  testmsg << "Comparison of A to DebugFieldPrint output:" << endl;
  dfp.print(A, testmsg);

  // test 2: assigment of index values
  testmsg << "\n\n-------------------------------------------" << endl;
  A[I][J] = J + 1;
  B[I][J] = I + 1;
  testmsg << "Results of A[I][J] = J+1:" << endl;
  printfield(testmsg, A);
  testmsg << "Results of B[I][J] = I+1:" << endl;
  printfield(testmsg, B);
  testmsg << "Comparison of A to DebugFieldPrint output:" << endl;
  dfp.print(A, testmsg);

  // test 3: assigment of index values
  testmsg << "\n\n-------------------------------------------" << endl;
  A[I][J] = I + J + 1;
  B[I][J] = I + J + 1;
  testmsg << "Results of A[I][J] = I+J+1:" << endl;
  printfield(testmsg, A);
  testmsg << "Results of B[I][J] = I+J+1:" << endl;
  printfield(testmsg, B);
  testmsg << "Comparison of A to DebugFieldPrint output:" << endl;
  dfp.print(A, testmsg);

  return 0;
}

/***************************************************************************
 * $RCSfile: gcprint.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: gcprint.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
