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
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "Ippl.h"

//////////////////////////////////////////////////////////////////////

void compute(double *a, double *b, int n)
{
  int i,j;
  int n2 = n+1;
  double *aa = a + 2 + n;
  for (j=0; j<n; j++)
    for (i=0; i<n; i++)
      b[i+n*j]=0.25*(aa[i-1+n2*j]+aa[i+1+n2*j]+aa[i+n2+n2*j]+aa[i-n2+n2*j]);
}

//////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[] )
{
  Ippl ippl(argc,argv);

  int N=atoi(argv[1]);
  int n_loops=atoi(argv[2]);

  Index I(N),J(N);
  const unsigned Dim=2;
  FieldLayout<Dim> layout(I,J);
  BConds<double,Dim> bc;
  bc[0] = new PeriodicFace<double,Dim>(0);
  bc[1] = new PeriodicFace<double,Dim>(1);
  bc[2] = new PeriodicFace<double,Dim>(2);
  bc[3] = new PeriodicFace<double,Dim>(3);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1)), B(layout);
  A = 0.0;
  B = 0.0;

  double *a = new double[(N+1)*(N+1)];
  double *b = new double[N*N];
  int i;
  for (i=0; i<(N+1)*(N+1); i++) a[i] = 0.0;
  for (i=0; i<N*N; i++) b[i] = 0.0;

  Timer tippl;
  Timer tc;

  tippl.clear();
  tc.clear();

  tippl.start();
  for (i=0; i<n_loops; i++)
    assign( B[I][J] , 0.25*(A[I+1][J]+A[I-1][J]+A[I][J+1]+A[I][J-1]) );
  tippl.stop();
  tc.start();
  for (i=0; i<n_loops; i++)
    compute(a,b,N);
  tc.stop();

  double flops = 4*n_loops*N*N;
  
  cout << "ippl mflops=" << 1e-6*flops/tippl.cpu_time() << endl;
  cout << "c     mflops=" << 1e-6*flops/tc.cpu_time() << endl;

  cout << "ippl usec=" << 1e6*tippl.cpu_time()/n_loops << endl;
  cout << "c     usec=" << 1e6*tc.cpu_time()/n_loops << endl;
}
/***************************************************************************
 * $RCSfile: speed1.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
 * IPPL_VERSION_ID: $Id: speed1.cpp,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
 ***************************************************************************/
