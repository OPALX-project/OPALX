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

#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "IndexedField.h"
#include "Field.h"
#include "Expressions.h"
#include "GuardCellSizes.h"
#include "Reduction.h"
#include "Communicate.h"

typedef Field<double,1> F;

/**********************************************************************/

int c1(Index& I, F& A, F& B, double b0)
{
  assign( A[I] , B[I]*b0 );
  return I.length();
}

int c2(Index& I, F& A, F& B, double b0, double b1)
{
  assign( A[I] , B[I]*b0 + B[I+1]*b1 );
  return I.length()*3;
}

int c3(Index& I, F& A, F& B, double b0, double b1, double b2)
{
  assign( A[I] , B[I]*b0 + B[I+1]*b1 + B[I+2]*b2 );
  return I.length()*5;
}

int c4(Index& I, F& A, F& B, double b0, double b1, double b2, double b3)
{
  assign( A[I] , B[I]*b0 + B[I+1]*b1 + B[I+2]*b2 + B[I+3]*b3 );
  return I.length()*7;
}

int c5(Index& I, F& A, F& B, double b0, double b1, double b2, double b3, double b4)
{
  assign( A[I] , B[I]*b0 + B[I+1]*b1 + B[I+2]*b2 + B[I+3]*b3 + B[I+4]*b4 );
  return I.length()*9;
}

int c6(Index& I, F& A,F& B,double b0,double b1,double b2,double b3,double b4,double b5)
{
  assign(A[I],B[I]*b0+B[I+1]*b1+B[I+2]*b2+B[I+3]*b3+B[I+4]*b4+B[I+5]*b5);
  return I.length()*11;
}

/**********************************************************************/

int main(int argc, char *argv[])
{
  IpplInit(argc,argv);

  int total_its = atoi(argv[1]);
  int min_len = atoi(argv[2]);
  int max_len = atoi(argv[3]);
  int len=min_len;
  Index domain(max_len+10);
  FieldLayout<1> layout(domain);
  F A(layout),B(layout);

  A[domain] = domain;
  B[domain] = domain;

  printf("%8s %10s %10s %10s %10s %10s %10s\n","length",
	 "mflops 1","mflops 2","mflops 3", "mflops 4","mflops 5","mflops 6");

  while (len<max_len) {
    clock_t tic0,tic1,tic2,tic3,tic4,tic5,tic6;
    double t1,t2,t3,t4,t5,t6;
    double mf1,mf2,mf3,mf4,mf5,mf6;
    int op1,op2,op3,op4,op5,op6;
    int its = total_its / len;
    int it,next_len;
    Index I(len);

    tic0 = clock();
    for (it=0; it<its; ++it)
      op1 = c1(I,A,B,1.0);
    tic1 = clock();
    for (it=0; it<its; ++it)
      op2 = c2(I,A,B,1.0,2.0);
    tic2 = clock();
    for (it=0; it<its; ++it)
      op3 = c3(I,A,B,1.0,2.0,3.0);
    tic3 = clock();
    for (it=0; it<its; ++it)
      op4 = c4(I,A,B,1.0,2.0,3.0,4.0);
    tic4 = clock();
    for (it=0; it<its; ++it)
      op5 = c5(I,A,B,1.0,2.0,3.0,4.0,5.0);
    tic5 = clock();
    for (it=0; it<its; ++it)
      op6 = c6(I,A,B,1.0,2.0,3.0,4.0,5.0,6.0);
    tic6 = clock();

    t1 = 1e6*(tic1-tic0)/(its*(double)CLOCKS_PER_SEC);
    t2 = 1e6*(tic2-tic1)/(its*(double)CLOCKS_PER_SEC);
    t3 = 1e6*(tic3-tic2)/(its*(double)CLOCKS_PER_SEC);
    t4 = 1e6*(tic4-tic3)/(its*(double)CLOCKS_PER_SEC);
    t5 = 1e6*(tic5-tic4)/(its*(double)CLOCKS_PER_SEC);
    t6 = 1e6*(tic6-tic5)/(its*(double)CLOCKS_PER_SEC);
    mf1 = op1/t1;
    mf2 = op2/t2;
    mf3 = op3/t3;
    mf4 = op4/t4;
    mf5 = op5/t5;
    mf6 = op6/t6;

//    printf("%8i %10g %10g %10g %10g %10g %10g\n",len,
//	   mf1,mf2,mf3,mf4,mf5,mf6);

    printf("%8i %10g %10g %10g %10g %10g %10g\n",len,
	   t1,t2,t3,t4,t5,t6);

    next_len = len*1.1;
    if ( next_len == len )
      len += 1;
    else
      len = next_len;
  }
}
/***************************************************************************
 * $RCSfile: stencil1.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:07:05 $
 * IPPL_VERSION_ID: $Id: stencil1.cpp,v 1.1.1.1 2000/11/30 21:07:05 adelmann Exp $ 
 ***************************************************************************/
