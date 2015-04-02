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
// ----------------------------------------------------------------------------
// The IPPL Framework - Visit http://people.web.psi.ch/adelmann/ for more details
// 
// This program was prepared by the Regents of the University of California at
// TestEureka.cpp
// Test Eureka boundary condition

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/CartesianCentering.h"
#include "AppTypes/Vektor.h"

#include <stdio.h>

// define some helper functions for printing out a Field

template <class M, class C>
void print(Field<double,2,M,C>& f)
{
  NDIndex<2> domain = f.getDomain();
  for (int j=domain[1].max(); j>=domain[1].min(); --j)
    {
      for (int i=domain[0].min(); i<=domain[0].max(); ++i)
	printf("%6.0f ",f[i][j].get());
      printf("\n");
    }
  return;
}

template <class M, class C>
void print(Field<Vektor<double,2>,2,M,C>& f)
{
  NDIndex<2> domain = f.getDomain();
  for (int j=domain[1].max(); j>=domain[1].min(); --j)
    {
      for (int i=domain[0].min(); i<=domain[0].max(); ++i)
	{
	  Vektor<double,2> x = f[i][j].get();
	  printf("(%2.0f,%2.0f)",x[0],x[1]);
	}
      printf("\n");
    }
  return;
}


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform out(argv[0]);
  const int N = 6;
  Index I(1,N), J(1,N);
  Index I0(-1,N+2), J0(-1,N+2);
  Index I1(2,N-1),J1(2,N-1);

  FieldLayout<2> layout(I,J);
  FieldLayout<2> layout0(I0,J0);
  GuardCellSizes<2> gc(2);

  // Test all cell centering.
  BConds<double,2> bc1;
  bc1[1] = new EurekaFace<double,2>(0);
  bc1[2] = new EurekaFace<double,2>(1);
  bc1[3] = new EurekaFace<double,2>(2);
  bc1[4] = new EurekaFace<double,2>(3);

  Field<double,2> A1(layout,gc,bc1), A0(layout0);

  A1[I][J] = 10.0*I + 100.0*J;
  //cout << "A1:" << endl;
  //print(A1);

  // Assign the full domain to A0.
  A0[I0][J0] = A1[I0][J0];

  // See if we got the right answer.
  // s1 should be zero.
  A0[I1][J1] -= 10.0*I1 + 100.0*J1;

  double s1 = sum(A0*A0);
  /*
  // Do pass/fail test below
  if ( s1 == 0 )
    out << "PASSED" << endl;
  else
    out << "FAILED" << endl;
  */

  // Test mixed centering.
  typedef Vektor<double,2> T;
  typedef UniformCartesian<2> M;
  typedef CartesianCentering<CCCEnums<2,2,0>::vectorFace,2,2> C;
  BConds<T,2,M,C> bc2;
  bc2[0] = new EurekaFace<T,2,M,C>(0,0);
  bc2[1] = new EurekaFace<T,2,M,C>(1,0);
  bc2[2] = new EurekaFace<T,2,M,C>(2,0);
  bc2[3] = new EurekaFace<T,2,M,C>(3,0);
  bc2[4] = new EurekaFace<T,2,M,C>(0,1);
  bc2[5] = new EurekaFace<T,2,M,C>(1,1);
  bc2[6] = new EurekaFace<T,2,M,C>(2,1);
  bc2[7] = new EurekaFace<T,2,M,C>(3,1);

  Field<T,2,M,C> B1(layout,gc,bc2), B0(layout0);

  // Fill with nontrivial data.
  B1[I][J] = T(1,1)*(I + 10.0*J);

  // Pull it out into a field that shows the guard layers.
  B0[I0][J0] = B1[I0][J0];

  // See if we got the right answer.
  B0[I1][J1] -= T(1,1)*(I1+10.0*J1);
  B0[1][J1] -= T(0,1)*(1+10.0*J1);
  B0[N][J1] -= T(0,1)*(N+10.0*J1);
  B0[I1][1] -= T(1,0)*(I1+10.0);
  B0[I1][N] -= T(1,0)*(I1+10.0*N);

  Vektor<T,2> s2 = sum(B0*B0);
  if ( s1 == 0 && s2 == Vektor<T,2>(0,0) )
    out << "PASSED" << endl;
  else
    out << "FAILED" << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: TestEureka.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: TestEureka.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
