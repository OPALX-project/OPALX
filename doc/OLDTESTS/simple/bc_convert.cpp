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

//
// NOTE: This test no longer compiles because the "convert_type()" method
// for BConds objects has been disabled.
//

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"


int main(int argc, char *argv[])
{
  Ippl(argc,argv);
  Inform testmsg(argv[0]);
  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  NDIndex<Dim> domain;
  domain[0] = I;
  domain[1] = J;
  FieldLayout<Dim> layout(domain);
  typedef UniformCartesian<Dim> Mesh;

  // Set initial boundary conditions.
  BConds<double,Dim,Mesh> Abc;
  Abc[0] = new PeriodicFace<double,Dim>(2);
  Abc[1] = new PeriodicFace<double,Dim>(3);            
  Abc[2] = new PosReflectFace<double,Dim>(0);
  Abc[3] = new ZeroFace<double,Dim>(1);
  // An array for the base.
  Field<double,Dim,Mesh> A(layout,GuardCellSizes<Dim>(2),Abc);
  
  // An array that gets its configuration from A.
  BConds<int,Dim,Mesh> *Bbc =
    (BConds<int,Dim,Mesh>*)A.getBConds().convert_type(int());
  Field<int,Dim,Mesh> B(A.getLayout(),A.getGuardCellSizes(),*Bbc);

  BConds<double,Dim,Mesh>::iterator pd;
  BConds<int,Dim,Mesh>::iterator pi;
  for (pd=Abc.begin(); pd!=Abc.end(); ++pd)
    testmsg << "A:" << (*pd).second[0] << endl;
  for (pi=Bbc->begin(); pi!=Bbc->end(); ++pi)
    testmsg << "B:" << (*pi).second[0] << endl;

  delete Bbc;

#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = I*10 + J ;
  B=A;

  A[I][J]=
	 A[I+1][J+1]+100.*A[I-1][J-1]+10000.*A[I-1][J+1]+1000000.*A[I+1][J-1];
  B[I][J]=
	 B[I+1][J+1]+100*B[I-1][J-1]+10000*B[I-1][J+1]+1000000*B[I+1][J-1];
#else
  A[I][J] << I*10 + J ;
  B << A;

  A[I][J]<<
	 A[I+1][J+1]+100.*A[I-1][J-1]+10000.*A[I-1][J+1]+1000000.*A[I+1][J-1];
  B[I][J]<<
	 B[I+1][J+1]+100*B[I-1][J-1]+10000*B[I-1][J+1]+1000000*B[I+1][J-1];
#endif

  testmsg << "A:" << A << endl;
  testmsg << "B:" << B << endl;
  A -= B;
  A *= A;
  double s = sum(A);
  testmsg << (s==0 ? "PASSED" : "FAILED") << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: bc_convert.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:32 $
 * IPPL_VERSION_ID: $Id: bc_convert.cpp,v 1.1.1.1 2000/11/30 21:05:32 adelmann Exp $ 
 ***************************************************************************/
