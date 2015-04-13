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

#include <iostream>
#include <math.h>
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
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout);
  Field<double,Dim> B(layout);
  Field<double,Dim> C(layout);

  A[I][J] << (I-1)*(I-1)+(J-1)*(J-1) + 1;
  NDIndex<Dim> maxloc,minloc;
  double maxval = max(A,maxloc);
  double minval = min(A,minloc);

  double known_max = 19;
  double known_min = 1;
  NDIndex<Dim> known_minloc(Index(1,1),Index(1,1));
  NDIndex<Dim> known_maxloc(Index(4,4),Index(4,4));
  bool test = true;
  test = test && (maxval == known_max);
  test = test && (minval == known_min);
  test = test && (maxloc == known_maxloc);
  test = test && (minloc == known_minloc);
  if ( test )
    testmsg << "PASSED" << endl;
  else
    testmsg << "FAILED" << endl;
}
/***************************************************************************
 * $RCSfile: reduceloc.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: reduceloc.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
