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
  Ippl ippl(argc,argv);

  

  int sizeX, sizeY, iterations;

  cout << "input sizeX" << endl;
  cin >> sizeX;
  cout << "input sizeY" << endl;
  cin >> sizeY;
  cout << "input the iterations" << endl;
  cin >> iterations;
   
  const unsigned Dim=2;
  Index I(sizeX);
  Index J(sizeY);
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));

  A = 0.0;
  A[sizeX/2][sizeY/2] = 1.0;

  double fact = 1.0/9.0;
  for(int iter = 0 ; iter < iterations ; iter++ ) 
  {
    A[I][J] = fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] + 
		    A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] + 
		    A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1]);
  }
  cout << "Done ; A[sizeX/2][sizeX/2] = " << A[sizeX/2][sizeX/2] << endl;
  return 0;
}
/***************************************************************************
 * $RCSfile: doof2d_simple.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: doof2d_simple.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
