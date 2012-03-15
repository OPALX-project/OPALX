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
  Field<double,Dim> B(layout,GuardCellSizes<Dim>(1));

  Field<double,Dim>::iterator p;
  int i;
  A = -1.0 ;
  for (p=B.begin(), i=0; p!=B.end(); ++p, ++i) *p = i+1;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  B= B*2.0;
#else
  assign( B , B*2.0 );
#endif

  testmsg << A << endl;
  testmsg << B << endl;

#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = (B[I+1][J+1]+B[I+1][J-1]+B[I-1][J+1]+B[I-1][J-1])/8.0;
#else
  assign(A[I][J] , (B[I+1][J+1]+B[I+1][J-1]+B[I-1][J+1]+B[I-1][J-1])/8.0);
#endif
  testmsg << A << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: t3.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:33 $
 * IPPL_VERSION_ID: $Id: t3.cpp,v 1.1.1.1 2000/11/30 21:06:33 adelmann Exp $ 
 ***************************************************************************/
