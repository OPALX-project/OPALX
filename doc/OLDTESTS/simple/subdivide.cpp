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
  Ippl ippl(argc,argv);

  const int Dim=2;
  int N=5, i, j;
  Index II(N+2),JJ(N+2);
  Index I(0,N),J(0,N);
  FieldLayout<Dim> layout(II,JJ,PARALLEL,PARALLEL,4);
  Field<double,Dim> A(layout),B(layout);

  A=-1.0;
  assign(B[II][JJ], II+JJ*10);
  A[I][J] = B[J+1][I+1];

  for (j=0; j<=N+1; ++j) {
    for(i=0; i<=N+1; ++i) 
      cout << " " << B[i][j];
    cout << endl;
  }
  cout << endl;

  for (Field<double,Dim>::iterator p=B.begin(); p!=B.end(); ++p)
    cout << " " << *p;
  cout << endl;

  for (j=0; j<=N+1; ++j) {
    for(i=0; i<=N+1; ++i) 
      cout << " " << A[i][j].get();
    cout << endl;
  }
  cout << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: subdivide.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:27 $
 * IPPL_VERSION_ID: $Id: subdivide.cpp,v 1.1.1.1 2000/11/30 21:06:27 adelmann Exp $ 
 ***************************************************************************/
