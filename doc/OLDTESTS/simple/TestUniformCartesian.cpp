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

// TestUniformCartesian.cpp
// Various tests of UniformCartesian class

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  const unsigned Dim3 = 3;
  const unsigned nx = 6, ny = 6 , nz = 6;
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  GuardCellSizes<Dim3> gc(1);
  Index I(nx);
  Index J(ny);
  Index K(nz);
  FieldLayout<Dim3> layoutVert(I,J,K);
  Index Ic(nx-1);
  Index Jc(ny-1);
  Index Kc(nz-1);
  FieldLayout<Dim3> layoutCell(Ic,Jc,Kc);

  // Test Average() functions:
  // Scalar Field, Cell-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> C(layoutCell,gc);
  C = 1.0;
  // Scalar weight Field, Cell-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> wC(layoutCell,gc);
  wC = 2.0;
  // Scalar Field, Vert-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> V(layoutVert,gc);
  V = 1.0;
  // Scalar weight Field, Vert-Centered
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> wV(layoutVert,gc);
  wV = 2.0;
  // Field's to hold weighted averages:
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> avgToC(layoutCell,gc);
  Field<double,Dim3,UniformCartesian<Dim3>,Vert> avgToV(layoutVert,gc);
  // Weighted average from Cell to Vert:
  assign(avgToV, Average(C, wC, avgToV));
  // Weighted average from Vert to Cell:
  assign(avgToC, Average(V, wV, avgToC));
  // Check results:
  if (sum(avgToV) != nx*ny*nz) {
    testmsg << "avgToC values wrong" << endl;
    passed = false;
  }
  if (sum(avgToC) != (nx-1)*(ny-1)*(nz-1)) {
    testmsg << "avgToC values wrong" << endl;
    passed = false;
  }

  // Try out Divergence Vektor/Vert -> Scalar/Cell:
  // Create mesh object:
  typedef UniformCartesian<Dim3,double> M;
  double delX[Dim3] = {1.0, 2.0, 3.0};
  Vektor<double,Dim3> origin;
  origin(0) = 1.0; origin(1) = 2.0; origin(2) = 3.0; 
  NDIndex<Dim3> verts;
  verts[0] = Index(nx); verts[1] = Index(nz); verts[2] = Index(nz); 

  //  M mesh(verts, delX, origin);
  // test constructing mesh first and adding origin, spacing later
  M mesh(verts);
  mesh.set_origin(origin);
  mesh.set_meshSpacing(delX);

  // Now use mesh object just created to construct Field's:
  Field<Vektor<double,Dim3>,Dim3,M,Vert> vectVert(mesh,layoutVert,gc);
  Field<double,Dim3,M,Cell> scalarCell(mesh,layoutCell,gc);
  // Assign values into the vert-centered Field<Vektor>:
  assign(vectVert, mesh.getVertexPositionField(vectVert));
  scalarCell = Div(vectVert, scalarCell);
  // The value should be 3.0 for all elements; test this:
  if (abs(sum(scalarCell)/((nx-1)*(ny-1)*(nz-1)) - 3.0) > 1.0e-6) {
    testmsg << "Div(vert position field) != const=3.0)" << endl;
    passed = false;
  }

  // Mechanism for constructing a Field from a mesh (two-step; must
  // first construct CenteredFieldLayout object):
  CenteredFieldLayout<Dim3,UniformCartesian<Dim3>,Cell> 
    cl(mesh, PARALLEL, SERIAL, SERIAL, 4);
  // Now construct a Field using it, (and the mesh object)
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> A(mesh,cl);
  assign(A, 5.0);

  // Some accessor function tests:
  double theVolume, theVolume2, theVolume3;
  NDIndex<Dim3> ndi;
  ndi[0] = Index(2,2,1);
  ndi[1] = Index(2,2,1);
  ndi[2] = Index(2,2,1);
  theVolume = mesh.getCellVolume(ndi);
  ndi[0] = Index(0,2,1);
  ndi[1] = Index(0,2,1);
  ndi[2] = Index(0,2,1);
  theVolume2 = mesh.getCellRangeVolume(ndi);
  if (theVolume2 != (3.0*6.0*9.0)) {
    testmsg << "volume of cells [0:2][0:2][0:2] != 162" << endl;
    testmsg << "volume of cells [0:2][0:2][0:2] = " << theVolume2 << endl;
    passed = false;
  }
  ndi[0] = Index(0,3,1);
  ndi[1] = Index(0,3,1);
  ndi[2] = Index(0,3,1);
  theVolume3 = mesh.getVertRangeVolume(ndi);
  if (theVolume3 != theVolume2) {
    testmsg << "volume within vertices [0:3][0:3][0:3] != "
	    << "vol of cells [0:2][0:2][0:2]" << endl;
    testmsg << "volume w/in vertices [0:3][0:3][0:3] = " << theVolume3 << endl;
    testmsg << "volume of cells [0:2][0:2][0:2] = " << theVolume2 << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------
  Field<double,Dim3,UniformCartesian<Dim3>,Cell> theVolumes(cl);
  mesh.getCellVolumeField(theVolumes);
  if ((sum(theVolumes)/((nx-1)*(ny-1)*(nz-1))) != theVolume) {
    testmsg << "(sum(theVolumes)/((nx-1)*(ny-1)*(nz-1))) != cell vol" << endl;
    testmsg << "(sum(theVolumes)/((nx-1)*(ny-1)*(nz-1))) = " 
	    << (sum(theVolumes)/((nx-1)*(ny-1)*(nz-1))) << endl;
    testmsg << "cell vol = " << theVolume << endl;
  }
  //---------------------------------------------------------------------------
  Vektor<double,Dim3> v;
  v(0) = 1.5; v(1) = 4.5; v(2) = 9.5;
  ndi = mesh.getNearestVertex(v);
  if ( ((ndi[0].first() != 1) || ndi[0].length() != 1) || 
       ((ndi[1].first() != 1) || ndi[1].length() != 1) || 
       ((ndi[2].first() != 2) || ndi[2].length() != 1) ) {
    testmsg << "NEAREST VERTEX TO (1.5,4.5,9.5) = " << ndi << " != (1,1,2)"
	    << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------
  Vektor<double,Dim3> v1;
  v1 = mesh.getVertexPosition(ndi);
  v(0) = 2.0; v(1) = 4.0; v(2) = 9.0; // Correct value
  if (v1 != v) {
    testmsg << "VERT POSITION OF" << ndi << " = " << v1 << " != " << v << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------
  CenteredFieldLayout<Dim3,UniformCartesian<Dim3>,Vert> 
    clVert(mesh);
  Field<Vektor<double,Dim3>,Dim3,UniformCartesian<Dim3>,Vert> 
    thePositions(clVert);
  mesh.getVertexPositionField(thePositions);
  //---------------------------------------------------------------------------
  v = mesh.getDeltaVertex(ndi);
  Vektor<double,Dim3> vcorrect;
  vcorrect(0) = 1.0; vcorrect(1) = 2.0; vcorrect(2) = 3.0;
  if (v != vcorrect) {
    testmsg << "DELTA-VERTEX OF" << ndi << " = " << v 
	    << " != " << vcorrect << endl;
    passed = false;
  }

  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: TestUniformCartesian.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:24 $
 * IPPL_VERSION_ID: $Id: TestUniformCartesian.cpp,v 1.1.1.1 2000/11/30 21:05:24 adelmann Exp $ 
 ***************************************************************************/
