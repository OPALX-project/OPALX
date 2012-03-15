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

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#include "Ippl.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform msg(argv[0]);

  const int Dim = 3;
  typedef double Typ;

  Index I(5), J(3), K(4);
  int i, j, k;

  // create a layout and field with 2 vnodes/processor
  FieldLayout<Dim> layout(I, J, K,
			  PARALLEL, PARALLEL, PARALLEL,
			  2*Ippl::getNodes());
  Field<Typ,Dim> A(layout);

  // assign values to 2D field
  A[I][J][K] = I + J*K;

  // for each slice of A, print out sum
  msg << "----------------------------------" << endl;
  for (i=0; i < I.length(); ++i)
    {
      Typ s = sum(A[i][J][K]);
      Typ comps = 0;
      for (j=0; j < J.length(); ++j)
	for (k=0; k < K.length(); ++k)
	  comps += (i + j*k);
      msg << "Sum of A[" << i << "][J][K] = " << s;
      msg << "  (should be " << comps << ") ... ";
      msg << (s == comps ? "PASSED" : "FAILED") << endl;
    }

  // for each slice of A, print out sum
  msg << "----------------------------------" << endl;
  for (j=0; j < J.length(); ++j)
    {
      Typ s = sum(A[I][j][K]);
      Typ comps = 0;
      for (i=0; i < I.length(); ++i)
	for (k=0; k < K.length(); ++k)
	  comps += (i + j*k);
      msg << "Sum of A[I][" << j << "][K] = " << s;
      msg << "  (should be " << comps << ") ... ";
      msg << (s == comps ? "PASSED" : "FAILED") << endl;
    }

  // for each slice of A, print out sum
  msg << "----------------------------------" << endl;
  for (k=0; k < K.length(); ++k)
    {
      Typ s = sum(A[I][J][k]);
      Typ comps = 0;
      for (i=0; i < I.length(); ++i)
	for (j=0; j < J.length(); ++j)
	  comps += (i + j*k);
      msg << "Sum of A[I][J][" << k << "] = " << s;
      msg << "  (should be " << comps << ") ... ";
      msg << (s == comps ? "PASSED" : "FAILED") << endl;
    }

  // test 1D slices
  msg << "----------------------------------" << endl;
  for (i=0; i < I.length(); ++i)
    {
      for (j=0; j < J.length(); ++j)
	{
	  Typ s = sum(A[i][j][K]);
	  Typ comps = 0;
	  for (k=0; k < K.length(); ++k)
	    comps += (i + j*k);
	  msg << "Sum of A[" << i << "][" << j << "][K] = " << s;
	  msg << "  (should be " << comps << ") ... ";
	  msg << (s == comps ? "PASSED" : "FAILED") << endl;
	}
    }

  msg << "----------------------------------" << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: sumslice.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:30 $
 * IPPL_VERSION_ID: $Id: sumslice.cpp,v 1.1.1.1 2000/11/30 21:06:30 adelmann Exp $ 
 ***************************************************************************/
