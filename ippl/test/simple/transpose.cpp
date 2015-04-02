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

#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"

const unsigned Dim=3;
Inform testmsg;

  
void check( Field<int,Dim>& f, int s1, int s2, int s3, int test)
{
  Index I = f.getIndex(0);
  Index J = f.getIndex(1);
  Index K = f.getIndex(2);
  f[I][J][K] -= s1*I + s2*J + s3*K;
  int sum_f = sum(f);
  if ( sum_f==0 )
    testmsg << "PASSED test " << test << "" << endl;
  else
    testmsg << "FAILED test " << test << " sum = " << sum_f << endl;
}

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  testmsg.setPrintNode(-1);

  const int N=10;

  Index I(N),J(N),K(N);
  FieldLayout<Dim> layout_ppp(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);
  FieldLayout<Dim> layout_spp(I,J,K,SERIAL,PARALLEL,PARALLEL,8);
  FieldLayout<Dim> layout_psp(I,J,K,PARALLEL,SERIAL,PARALLEL,8);
  FieldLayout<Dim> layout_pps(I,J,K,PARALLEL,PARALLEL,SERIAL,8);

  Field<int,Dim> A(layout_ppp);
  Field<int,Dim> B(layout_spp);
  Field<int,Dim> C(layout_psp);
  Field<int,Dim> D(layout_pps);
  int ii = 1;
  int jj = 10;
  int kk = 100;
  
  assign( A[I][J][K] , ii*I + jj*J + kk*K );
  B = A;
  C = A;
  D = A;
  check(B,ii,jj,kk,1);
  check(C,ii,jj,kk,2);
  check(D,ii,jj,kk,3);

  B[I][J][K] = A[I][J][K];
  C[I][J][K] = A[I][J][K];
  D[I][J][K] = A[I][J][K];
  check(B,ii,jj,kk,4);
  check(C,ii,jj,kk,5);
  check(D,ii,jj,kk,6);

  B[I][K][J] = A[I][J][K];
  C[I][K][J] = A[I][J][K];
  D[I][K][J] = A[I][J][K];
  check(B,ii,kk,jj,7);
  check(C,ii,kk,jj,8);
  check(D,ii,kk,jj,9);

  B[J][I][K] = A[I][J][K];
  C[J][I][K] = A[I][J][K];
  D[J][I][K] = A[I][J][K];
  check(B,jj,ii,kk,10);
  check(C,jj,ii,kk,11);
  check(D,jj,ii,kk,12);

  B[J][K][I] = A[I][J][K];
  C[J][K][I] = A[I][J][K];
  D[J][K][I] = A[I][J][K];
  check(B,jj,kk,ii,13);
  check(C,jj,kk,ii,14);
  check(D,jj,kk,ii,15);

  B[K][I][J] = A[I][J][K];
  C[K][I][J] = A[I][J][K];
  D[K][I][J] = A[I][J][K];
  check(B,kk,ii,jj,16);
  check(C,kk,ii,jj,17);
  check(D,kk,ii,jj,18);

  B[K][J][I] = A[I][J][K];
  C[K][J][I] = A[I][J][K];
  D[K][J][I] = A[I][J][K];
  check(B,kk,jj,ii,19);
  check(C,kk,jj,ii,20);
  check(D,kk,jj,ii,21);
}
/***************************************************************************
 * $RCSfile: transpose.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: transpose.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
