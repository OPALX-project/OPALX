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
  A simple program to test SubField assignments using two fields, one on
  a cell-centered layout, the other on a vertex-centered (but otherwise
  similar) layout.
 ***************************************************************************/

int main(int argc, char *argv[]) {
  TAU_PROFILE("main()", "int (int, char **)", TAU_PROFILE);
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  const unsigned Dim=2;
  Index I(4), J(4), U(5), V(5);
  FieldLayout<Dim> LC(I, J, PARALLEL, PARALLEL, 4);
  FieldLayout<Dim> LV(U, V, PARALLEL, PARALLEL, 4);

  Field<int,Dim,UniformCartesian<Dim>,Cell> A(LC, GuardCellSizes<Dim>(1));
  Field<int,Dim,UniformCartesian<Dim>,Vert> B(LV, GuardCellSizes<Dim>(1));

  A[I][J] = J + 1;
  B[U][V] = -V - 1;
  testmsg << "A at start: " << A << endl;
  testmsg << "B at start: " << B << endl;

  SIndex<Dim> SI(LC);
  SI = eq(A,2) || eq(A,4);
  testmsg << "SI at start: " << SI << endl;

  A[SI] = A[SI(0,-1)] + B[SI(1,1)];
  testmsg << "A after assignment: " << A << endl;
  testmsg << "B after assignment: " << B << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: vertcell.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:56 $
 * IPPL_VERSION_ID: $Id: vertcell.cpp,v 1.1.1.1 2000/11/30 21:06:56 adelmann Exp $ 
 ***************************************************************************/
