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

#include <iostream.h>
#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(domain);
  Field<double,Dim> A(layout);
  Field<double,Dim> B(layout);

  A = 0;

  A[I][J] += (I-2)*(I-2) + (J-2)*(J-2);
  testmsg << A << endl;
  A -= 1.0;
  testmsg << A << endl;
  B[I][J] = I;
  testmsg << B << endl;
  A *= B;
  testmsg << A << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: t5.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:35 $
 * IPPL_VERSION_ID: $Id: t5.cpp,v 1.1.1.1 2000/11/30 21:06:35 adelmann Exp $ 
 ***************************************************************************/
