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

double neg(const double& x) { return -x; }
double twice(const double& x) { return 2.0*x; }
double zero(const double&) { return 0.0; }
double one(const double&) { return 1.0; }

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=2;
  Index I(5),J(5);
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> B(layout);

  // Set initial boundary conditions.
  BConds<double,Dim> bc;
  bc[0] = new FunctionFace<double,Dim>(neg,0);
  bc[1] = new FunctionFace<double,Dim>(twice,1);
  bc[2] = new FunctionFace<double,Dim>(zero,2);
  bc[3] = new FunctionFace<double,Dim>(one,3);

  // An array for testing.
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1),bc);

  A = 0.0;
  A[I][J] += I + J*10;
  testmsg << A << endl;
  B[I][J] = A[I-1][J-1];
  testmsg << B << endl;
  B[I][J] = A[I+1][J+1];
  testmsg << B << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: funcface.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: funcface.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
