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

#ifdef IPPL_USE_STANDARD_HEADERS
#include <iostream>
#include <fstream>
using namespace std;
#else
#include <iostream.h>
#include <fstream.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

ostream *logfile;

//////////////////////////////////////////////////////////////////////

double mode(int i1, int i2, Field<double,2>& A,
	    Field<double,2>& X, Field<double,2>& Y)
{
  NDIndex<2> domain = A.getDomain();
  double w1 = 2.0/domain[0].length();
  double w2 = 2.0/domain[1].length();
  Field<double,2U> B( A.getLayout() );
#ifdef IPPL_USE_MEMBER_TEMPLATES
  B = A*sin(i1*X)*sin(i2*Y)*(w1*w2) ;
#else
  assign(B , A*sin(i1*X)*sin(i2*Y)*(w1*w2));
#endif
  return sum(B);
}

//////////////////////////////////////////////////////////////////////

void modes(Field<double,2>& A, Field<double,2>& X, Field<double,2>& Y)
{
  cout<< " "
    << mode(1,1,A,X,Y) << " "
    << mode(2,1,A,X,Y) << " "
    << mode(1,2,A,X,Y) << " " 
    << mode(2,2,A,X,Y) << endl;
}

void print1(Field<double,2>& A, Inform& testmsg)
{
  Field<double,2>::iterator p = A.begin();
  int i,j;
  NDIndex<2> idx = A.getDomain();
  for (j=idx[1].max(); j>=idx[1].min(); --j) {
    for (i=idx[0].min(); i<=idx[0].max(); ++i) {
      printf(" %10.6f",*p);
      ++p;
    }
    printf("\n");
  }
  printf("\n");
}

void print(Field<double,2>& A, Inform& testmsg)
{
  Field<double,2>::iterator p = A.begin();
  
  NDIndex<2> idx = A.getLocalDomain();
  for (int j=idx[1].max(); j>=idx[1].min(); --j) {
    for (int i=idx[0].min(); i<=idx[0].max(); ++i) {
      testmsg << *p << " ";
      ++p;
    }
  }
  testmsg << endl;
}

int main(int argc, char *argv[] )
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  char name[100];
  sprintf(name,"out.%d",Ippl::Comm->myNode());
  logfile = new ofstream(name);
  *logfile << "Starting to run!" << endl;
  Ippl::Comm->barrier();

  const int N=5;
  Index I(N),J(N);
  const unsigned Dim=2;

  // Set initial boundary conditions.
  BConds<double,Dim> bc;
  bc[0] = new NegReflectFace<double,Dim>(0);
  bc[1] = new NegReflectFace<double,Dim>(1);
  bc[2] = new NegReflectFace<double,Dim>(2);
  bc[3] = new NegReflectFace<double,Dim>(3);

  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1),bc);
  Field<double,Dim> X(layout),Y(layout);
  Field<double,Dim> B(layout);

  const double Pi = 4.0*atan(1.0);
#ifdef IPPL_USE_MEMBER_TEMPLATES
  X[I][J] = Pi*I/(double(N)) ;
  Y[I][J] = Pi*J/(double(N)) ;
  A = sin(X)*sin(Y) ;
#else
  assign(X[I][J] , Pi*I/(double(N)));
  assign(Y[I][J] , Pi*J/(double(N)));
  assign(A , sin(X)*sin(Y));
#endif
  A = 0.0;
  A[N/2][N/2] = 1.0;

  int istep = 0;
  testmsg << A << endl;
  for (istep=1; istep<=4; istep++) {
#ifdef IPPL_USE_MEMBER_TEMPLATES
    B[I][J] = 0.25*(A[I+1][J] + A[I-1][J] + A[I][J-1] + A[I][J+1]) ;
#else
    assign(B[I][J] , 0.25*(A[I+1][J] + A[I-1][J] + A[I][J-1] + A[I][J+1]) );
#endif
    testmsg << B << endl;
  }

  print(A,testmsg);
  print(B,testmsg);
}

/***************************************************************************
 * $RCSfile: p3.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:16 $
 * IPPL_VERSION_ID: $Id: p3.cpp,v 1.1.1.1 2000/11/30 21:06:16 adelmann Exp $ 
 ***************************************************************************/
