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

// include files 

#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"
#include "Field/BCond.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  const int N=5;
  Index I(N), J(N);
  BConds<double,2> bc;
  if (Ippl::getNodes() == 1) {
    bc[0] = new PeriodicFace<double,2>(0);
    bc[1] = new PeriodicFace<double,2>(1);
    bc[2] = new PeriodicFace<double,2>(2);
    bc[3] = new PeriodicFace<double,2>(3);
  }
  else {
    bc[0] = new ParallelPeriodicFace<double,2>(0);
    bc[1] = new ParallelPeriodicFace<double,2>(1);
    bc[2] = new ParallelPeriodicFace<double,2>(2);
    bc[3] = new ParallelPeriodicFace<double,2>(3);
  }

  FieldLayout<1> layout1(I);
  FieldLayout<2> layout2(I,J);
  Field<double,2> B(layout2);
  Field<double,1> C(layout1);
  Field<double,2> T2(layout2);
  Field<double,1> T1(layout1);

  int Guards = 0;
  int i,j;

  {
    Field<double,2> A(layout2,GuardCellSizes<2>(Guards),bc);

    //----------------------------------------

    assign(A[I][J] , I+J*10);
    T2 = A;
    for (i=0;i<N;++i)
      for (j=0;j<N;++j) {
	T2[i][j] -= i+j*10.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " initializing A" << endl;

    //----------------------------------------

    B = -1.0;
    B[I][J] = A[J][I];
    T2 = B;
    for (i=0;i<N;++i)
      for (j=0;j<N;++j) {
	T2[j][i] -= i+j*10.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " transposing A" << endl;

    //----------------------------------------

    B = -1.0;
    B[I][J] = A[I+1][J+1];
    T2 = B;
    for (i=0;i<N;++i) 
      for (j=0;j<N;++j)  {
	if ( (i<N-1)&&(j<N-1) )
	  T2[i][j] -= (i+1)+(j+1)*10.0;
	else 
	  T2[i][j] += 1.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " shifting A" << endl;

    //----------------------------------------

    B = -1.0;
    B[4][J] = A[J][4];
    T2 = B;
    for (i=0;i<N;++i) {
      for (j=0;j<N;++j) {
	if ( i==4 )
	  T2[i][j] -= j+40.0;
	else
	  T2[i][j] += 1.0;
      }
    }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " copying a slice" << endl;

    //----------------------------------------

    C[I] = 0.1*I;
    B = -1.0;
    B[I][4] = C[I] ;
    T2 = B;
    for (i=0;i<N;++i) {
      for (j=0;j<N;++j) {
	if ( j==4 )
	  T2[i][j] -= i*0.1;
	else
	  T2[i][j] += 1.0;
      }
    }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " inserting a slice" << endl;

    //----------------------------------------

    C[I] = A[2][I] ;
    for (i=0; i<N; ++i)
      C[i] -= 2.0 + i*10.0;
    C *= C;
    if ( sum(C) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " extracting a slice" << endl;

  }

  // Now the same tests with 1 layer of guard cells and periodic bc.
  // The answers for shifting are slightly different.
  Guards = 1;

  {
    Field<double,2> A(layout2,GuardCellSizes<2>(Guards),bc);

    //----------------------------------------

    A[I][J] = I+J*10 ;
    T2 = A;
    for (i=0;i<N;++i)
      for (j=0;j<N;++j) {
	T2[i][j] -= i+j*10.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " initializing guarded A" << endl;

    //----------------------------------------

    B = -1.0;
    B[I][J] = A[J][I];
    T2 = B;
    for (i=0;i<N;++i)
      for (j=0;j<N;++j) {
	T2[j][i] -= i+j*10.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " transposing A" << endl;

    //----------------------------------------

    B = -1.0;
    B[I][J] = A[I+1][J+1];
    T2 = B;
    for (i=0;i<N;++i) 
      for (j=0;j<N;++j)  {
	int ii = (i+1)%N;
	int jj = (j+1)%N;
	T2[i][j] -= ii+jj*10.0;
      }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " shifting A" << endl;

    //----------------------------------------

    B = -1.0;
    B[4][J] = A[J][4];
    T2 = B;
    for (i=0;i<N;++i) {
      for (j=0;j<N;++j) {
	if ( i==4 )
	  T2[i][j] -= j+40.0;
	else
	  T2[i][j] += 1.0;
      }
    }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " copying a slice" << endl;

    //----------------------------------------

    C[I] = 0.1*I;
    B = -1.0;
    B[I][4] = C[I] ;
    T2 = B;
    for (i=0;i<N;++i) {
      for (j=0;j<N;++j) {
	if ( j==4 )
	  T2[i][j] -= i*0.1;
	else
	  T2[i][j] += 1.0;
      }
    }
    T2 *= T2;
    if ( sum(T2) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " inserting a slice" << endl;

    //----------------------------------------

    C[I] = A[2][I] ;
    for (i=0; i<N; ++i)
      C[i] -= 2.0 + i*10.0;
    C *= C;
    if ( sum(C) != 0 )
      testmsg << "Failed";
    else
      testmsg << "PASSED";
    testmsg << " extracting a slice" << endl;

  }
}
/***************************************************************************
 * $RCSfile: patches.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:17 $
 * IPPL_VERSION_ID: $Id: patches.cpp,v 1.1.1.1 2000/11/30 21:06:17 adelmann Exp $ 
 ***************************************************************************/
