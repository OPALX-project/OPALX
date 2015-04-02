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

  const int N=5;
  Index I(N), J(N);

  FieldLayout<2> layout2(I,J);
  Field<double,2> B(layout2),A(layout2);
  Field<float,2>  C(layout2),D(layout2);
  double d=1;
  float  f=2;
  int i=3;

  A = 1;
  B = 2.0;
  C = 3;
  D = 4.0;

  A = f;
  B = d;
  C = d;
  D = f;

  B += 2.0*A;			// d = d*d
  B += 2.0*C;			// d = d*f
  B += 2.0F*A;			// d = f*d
  B += 2.0F*C;			// d = f*f
  B += 2*A;			// d = i*d
  B += 2*C;			// d = i*f
  C += 2.0*A;			// f = d*d
  C += 2.0*C;			// f = d*f
  C += 2.0F*A;			// f = f*d
  C += 2.0F*C;			// f = f*f
  C += 2*A;			// f = f*d
  C += 2*C;			// f = f*f

  B += d*A;			// d = d*d
  B += d*C;			// d = d*f
  B += f*A;			// d = f*d
  B += f*C;			// d = f*f
  B += i*A;			// d = i*d
  B += i*C;			// d = i*f
  C += d*A;			// f = d*d
  C += d*C;			// f = d*f
  C += f*A;			// f = f*d
  C += f*C;			// f = f*f
  C += i*A;			// f = i*d
  C += i*C;			// f = i*f

  testmsg << "Results:" << endl;
  testmsg << "A = " << A << endl;
  testmsg << "B = " << B << endl;
  testmsg << "C = " << C << endl;
}
/***************************************************************************
 * $RCSfile: float.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: float.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
