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
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  const unsigned Dim=2;
  const int N=5;
  const int Guards=1;
  Index I(N);
  Index J(N);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL,2,2,true);
  Field<int,Dim> A(layout),AA(layout);

  BConds<int,Dim> bc;
  if (Ippl::getNodes() == 1) {
    bc[0] = new PeriodicFace<int,Dim>(0);
    bc[1] = new PeriodicFace<int,Dim>(1);
    bc[2] = new PeriodicFace<int,Dim>(2);
    bc[3] = new PeriodicFace<int,Dim>(3);
  }
  else {
    bc[0] = new ParallelPeriodicFace<int,Dim>(0);
    bc[1] = new ParallelPeriodicFace<int,Dim>(1);
    bc[2] = new ParallelPeriodicFace<int,Dim>(2);
    bc[3] = new ParallelPeriodicFace<int,Dim>(3);
  }
  Field<int,Dim> B0(layout, GuardCellSizes<Dim>(Guards), bc);
  Field<int,Dim> B1(layout);

  A = 0.0;
  B0 = 0.0;
  B1 = 0.0;
  int errors = 0;
  for (int i0=0; i0<N; ++i0)
    for (int j0=0; j0<N; ++j0) {
      B0 = 0;
      B0[i0][j0] = 1;

      for (int i1=-Guards; i1<N+Guards; ++i1)
	for (int j1=-Guards; j1<N+Guards; ++j1) {
	  int i2 = (i1+N)%N;
	  int j2 = (j1+N)%N;
	  double v = ( (i2==i0) && (j2==j0) ? 1 : 0 );
	  double compv = B0[i1][j1].get();
	  if ( compv != v ) {
	    ++errors;
	    testmsg << "Test FAILED: v=" << v << " != B0[" << i1 << "][";
	    testmsg << i2 << "]=" << compv << endl;
	    testmsg << "   (i0,j0)=(" << i0 << "," << j0 << ")" << endl;
	    testmsg << "   (i1,j1)=(" << i1 << "," << j1 << ")" << endl;
	    testmsg << "   (i2,j2)=(" << i2 << "," << j2 << ")" << endl;
	  } else {
#ifndef __MWERKS__ // Output overflows WinSIOUS buffer
	    testmsg << "Test PASSED" << v << endl;
#endif
	  }
	}
    }
  if ( errors != 0 ) {
    testmsg << "test FAILED for " << errors << " cases." << endl;
  } else {
    testmsg << "all tests PASSED" << endl;
  }
}
/***************************************************************************
 * $RCSfile: drive3.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:03 $
 * IPPL_VERSION_ID: $Id: drive3.cpp,v 1.1.1.1 2000/11/30 21:06:03 adelmann Exp $ 
 ***************************************************************************/
