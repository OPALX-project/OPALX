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
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"

// const unsigned D1 = 1;
const unsigned D2 = 2;
const unsigned D3 = 3;
const unsigned D4 = 4;
// const unsigned D5 = 5;
// const unsigned D6 = 6;
Inform testmsg;

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  testmsg.setPrintNode(-1);

  const int N=10;

  Index I(N),J(N),K(N),L(N);
  NDIndex<D2> IJ(I,J);
  NDIndex<D3> IJK(IJ,K);
  NDIndex<D4> IJKL(IJK,L);
  e_dim_tag dist[D4] = { PARALLEL,PARALLEL,SERIAL,SERIAL };
  FieldLayout<D2> layout2(IJ,dist);
  FieldLayout<D3> layout3(IJK,dist);
  FieldLayout<D4> layout4(IJKL,dist);

  Field<int,D2> A1(layout2),A2(layout2);
  Field<int,D3> B1(layout3),B2(layout3);
  Field<int,D4> C1(layout4),C2(layout4);
  int ii = 1;
  int jj = 10;
  int kk = 100;
  int ll = 1000;
  int s;
  
  A1[I][J] = ii*I + jj*J;
  A2[J][I] = A1[I][J] ;
  A2[I][J] -= ii*J + jj*I ;
  A2 *= A2;
  s = sum(A2);
  testmsg << "2D General transpose: " << (s?"FAILED":"PASSED") << endl;

  B1[I][J][K] = ii*I + jj*J + kk*K;
  B2[I][J][K] = B1[I][K][J] ;
  B2[I][J][K] -= ii*I + jj*K + kk*J;
  B2 *= B2;
  s = sum(B2);
  testmsg << "3D General transpose: " << (s?"FAILED":"PASSED") << endl;

  C1[I][J][K][L] = ii*I + jj*J + kk*K + ll*L ;
  C2[I][J][L][K] = C1[I][J][K][L] ;
  C2[I][J][K][L] -= ii*I + jj*J + kk*L + ll*K ;
  C2 *= C2;
  s = sum(C2);
  testmsg << "4D serial transpose: "<< (s?"FAILED":"PASSED") << endl;
      
  int a,b;
  s = 0;
  for (a=0; a<N; ++a)
    {
      A1[I][J] = B1[a][I][J];
      A1[I][J] -= a*ii + I*jj + J*kk;
      A1 *= A1;
      s += sum(A1);
    }
  testmsg << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << endl;

  s = 0;
  for (a=0; a<N; ++a)
    {
      A1[I][J] = B1[I][a][J];
      A1[I][J] -= I*ii + a*jj + J*kk;
      A1 *= A1;
      s += sum(A1);
    }
  testmsg << "General 2 from 3 " << (s ? "FAILED " : "PASSED ") << a << endl;

  s = 0;
  for (a=0; a<N; ++a)
    {
      A1[I][J] = B1[I][J][a];
      A1[I][J] -= I*ii + J*jj + a*kk;
      A1 *= A1;
      s += sum(A1);
    }
  testmsg << "Serial 2 from 3 " << (s ? "FAILED " : "PASSED ") << endl;

  A1[I][J] = ii*I + jj*J;
  for (a=0; a<N; ++a)
    B1[a][I][J] = A1[I][J] ;
  B1[K][I][J] -= I*ii + J*jj;
  B1 *= B1;
  s = sum(B1);
  testmsg << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

  for (a=0; a<N; ++a)
    B1[I][a][J] = A1[I][J] ;
  B1[I][K][J] -= I*ii + J*jj;
  B1 *= B1;
  s = sum(B1);
  testmsg << "General 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

  for (a=0; a<N; ++a)
    B1[I][J][a] = A1[I][J] ;
  B1[I][J][K] -= I*ii + J*jj;
  B1 *= B1;
  s = sum(B1);
  testmsg << "Serial 3 from 2 " << (s ? "FAILED " : "PASSED ") << endl;

  for (a=0; a<N; ++a)
    for (b=0; b<N; ++b)
      C1[I][J][a][b] = A1[I][J] ;
  C1[I][J][K][L] -= I*ii + J*jj;
  C1 *= C1;
  s = sum(C1);
  testmsg << "Serial 4 from 2 " << (s ? "FAILED " : "PASSED ") << endl;
}
/***************************************************************************
 * $RCSfile: slices.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: slices.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
