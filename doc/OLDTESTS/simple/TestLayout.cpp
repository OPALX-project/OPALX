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

// this test files checks that the output from three conforming Fields
// with different layouts produces the same data

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  int vnodes = 4;
  int size = 8;
  const unsigned Dim=2;
  Index I(size), J(size);
  FieldLayout<Dim> layout1(I,J,PARALLEL,PARALLEL,vnodes);
  Field<double,Dim> A(layout1);
  Field<double,Dim> Acomp(layout1);
  FieldLayout<Dim> layout2(I,J,SERIAL,PARALLEL,vnodes);
  Field<double,Dim> B(layout2);
  FieldLayout<Dim> layout3(I,J,PARALLEL,SERIAL,vnodes);
  Field<double,Dim> C(layout3);

  // Create a user-specified layout as well
  NDIndex<Dim> *domlist = new NDIndex<Dim>[vnodes];
  domlist[0] = NDIndex<Dim>(Index(0,7), Index(0,2));
  domlist[1] = NDIndex<Dim>(Index(0,1), Index(3,5));
  domlist[2] = NDIndex<Dim>(Index(2,7), Index(3,5));
  domlist[3] = NDIndex<Dim>(Index(0,7), Index(6,7));
  int *nodelist = new int[vnodes];
  for (int n=0; n < vnodes; ++n)
    nodelist[n] = n % Ippl::getNodes();
  FieldLayout<Dim> layout4(NDIndex<Dim>(I,J), domlist, domlist + vnodes,
			   nodelist, nodelist + vnodes);
  Field<double,Dim> D(layout4);

  // Also create a CenteredFieldLayout as a test
  UniformCartesian<Dim> mesh(NDIndex<Dim>(I,J));
  CenteredFieldLayout<Dim, UniformCartesian<Dim>, Vert> cfl(mesh,
	domlist, domlist + vnodes, nodelist, nodelist + vnodes);

  // clean up lists of domains and nodes
  delete [] domlist;
  delete [] nodelist;

  testmsg << "Layout 1:" << endl;
  testmsg << "---------" << endl;
  testmsg << layout1 << endl << endl;

  testmsg << "Layout 2:" << endl;
  testmsg << "---------" << endl;
  testmsg << layout2 << endl << endl;

  testmsg << "Layout 3:" << endl;
  testmsg << "---------" << endl;
  testmsg << layout3 << endl << endl;

  testmsg << "Layout 4:" << endl;
  testmsg << "---------" << endl;
  testmsg << layout4 << endl << endl;

  testmsg << "Layout cfl:" << endl;
  testmsg << "-----------" << endl;
  testmsg << cfl << endl << endl;

  A[I][J] = I + 100*J;
  testmsg << "A = " << A << endl;
  //A.write("atest");

  B[I][J] = I + 100*J;
  testmsg << "B = " << B << endl;
  //B.write("btest");

  C[I][J] = I + 100*J;
  testmsg << "C = " << C << endl;
  //C.write("ctest");

  D[I][J] = I + 100*J;
  testmsg << "D = " << D << endl;
  //D.write("dtest");

  Acomp = B;
  double bsum = sum( (A - Acomp) * (A - Acomp) );
  Acomp = C;
  double csum = sum( (A - Acomp) * (A - Acomp) );
  Acomp = D;
  double dsum = sum( (A - Acomp) * (A - Acomp) );
  testmsg << "B OK ? " << (bsum == 0) << " (MSD = " << bsum << ")" << endl;
  testmsg << "C OK ? " << (csum == 0) << " (MSD = " << csum << ")" << endl;
  testmsg << "D OK ? " << (dsum == 0) << " (MSD = " << dsum << ")" << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: TestLayout.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:13 $
 * IPPL_VERSION_ID: $Id: TestLayout.cpp,v 1.1.1.1 2000/11/30 21:05:13 adelmann Exp $ 
 ***************************************************************************/
