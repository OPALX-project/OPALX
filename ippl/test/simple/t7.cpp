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
#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=1;
  Index I(10);				// Index on [0..9]
  Index Low(5);				// Index on [0..4]
  Index High(5,9);			// Index on [5..9]
  Index Even(0,9,2);			// Index on [0..9] stride 2
  Index Odd(1,9,2);			// Index on [1..9] stride 2
  Index Edge(0,9,9);                    // Index on [0..9] stride 9; the edge cells
  
  NDIndex<Dim> domain;
  domain[0] = I;
  FieldLayout<Dim> layout(domain);
  Field<double,Dim> A(layout);

  A = 0;

  testmsg << A << endl;

  A[I-5] = 1.0;

  testmsg << A << endl;

  A[I+5] = 2.0;

  testmsg << A << endl;

  A[I*2] -= 1.0;

  testmsg << A << endl;

  A[I*2+1] -= 1.0;

  testmsg << A << endl;
  
  A[I*9] = 3.0;
  
  testmsg << A << endl;

  testmsg << "Should be the same as :" << endl;

  A = 0;

  testmsg << A << endl;

  A[Low] = 1.0;

  testmsg << A << endl;

  A[High] = 2.0;

  testmsg << A << endl;

  A[Even] -= 1.0;

  testmsg << A << endl;

  A[Odd] -= 1.0;

  testmsg << A << endl;
  
  A[Edge] = 3.0;
  
  testmsg << A << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: t7.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: t7.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
