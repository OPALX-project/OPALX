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
#include <stdlib.h>

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  if (argc != 3) {
    testmsg << "Usage: " << argv[0] << " <nx> <ny>" << endl;
    exit(1);
  }

  const unsigned Dim=2;
  Index I(atoi(argv[1]));
  Index J(atoi(argv[2]));
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(domain);
  Field<double,Dim> A(layout);

  A[I][J] << I+1;

  Timer timer;
  double minv, maxv;

  timer.clear(); timer.start();
  minv = min(A);
  timer.stop();
  testmsg << "min(A) = " << minv << ", cputime = " << timer.cpu_time() << endl;

  timer.clear(); timer.start();
  maxv = max(A);
  timer.stop();
  testmsg << "max(A) = " << maxv << ", cputime = " << timer.cpu_time() << endl;

  timer.clear(); timer.start();
  minmax(A, minv, maxv);
  timer.stop();
  testmsg << "minmax(A) = (" << minv << ", " << maxv << ")";
  testmsg << ", cputime = " << timer.cpu_time() << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: reducespeed.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:07:03 $
 * IPPL_VERSION_ID: $Id: reducespeed.cpp,v 1.1.1.1 2000/11/30 21:07:03 adelmann Exp $ 
 ***************************************************************************/
