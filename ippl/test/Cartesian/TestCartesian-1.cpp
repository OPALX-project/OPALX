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

// TestCartesian.cpp
// Various tests of Cartesian class

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Field/GuardCellSizes.h"
#include "Meshes/Cartesian.h"
#include "AppTypes/Vektor.h"
#include "Utility/FieldDebug.h"


// Some of the tests are specific to 3D, meaning you'll get a compile error if
// you try to run 2D or 1D. By using this C-preprocessor macro to specify the
// dimensionality, the problematic tests are excluded if you set DIM to
// something other than 3:

#define DIM 3

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  Inform pout(NULL,0);
  setInform(pout);
  setFormat(2,3);
  bool passed = true; // Pass/fail test

  const unsigned D = DIM; // Hardwire dimensionality
  const unsigned nv = 6;  // Hardwire number of vertices in every direction
  unsigned vnodes = 4;    // Hardwire 4 vnodes

  double roundOffError = 1.0e-10;

  // Sizes:
  unsigned nverts[D], ncells[D];
  unsigned totverts=1, totcells=1;
  int d;
  for (d=0; d<D; d++) {
    ncells[d] = nv - 1;
    nverts[d] = nv;
    totcells *= ncells[d];
    totverts *= nverts[d];
  }
  NDIndex<D> verts, cells;
  for (d=0; d<D; d++) {
    verts[d] = Index(nverts[d]);
    cells[d] = Index(ncells[d]);
  }

  //---------------------------------------------------------------------------
  // Construct some CenteredFieldLayout's and Field's to be used below:

  // Create cartesian mesh object:
  typedef Cartesian<D,double> M;
  
  double* delX[D];

  for (d=0; d<D; d++)
      delX[d] = new double[nverts[d]];

  Vektor<double,D> origin;
  for (d=0; d<D; d++) 
      origin(d) = d + 1.0;

  // Assign nonuniform mesh-spacing values to each component (linear ramps):
  for (d=0; d<D; d++) {
      double multipplier = (d + 1)*1.0;
      for (int vert=0; vert < nverts[d]; vert++) {
	  (delX[d])[vert] = multipplier*(1 + vert);
      }
  }
  
  // Mesh boundary conditions:
  MeshBC_E mbc[2*D];
  for (unsigned b=0; b < (2*D); b++) 
      mbc[b] = Reflective;

  // Test constructing mesh, and then setting spacing, origin, BC's
  M mesh(verts);
  mesh.set_origin(origin);
  mesh.set_meshSpacing(delX);
  mesh.set_MeshBC(mbc);

  // Clean up mesh spacing arrays
  for (d=0; d<D; d++)
      delete [] delX[d];

  // ada have to cross check Div() fails without this
  mesh.storeSpacingFields(); 

  BareField<Vektor<double,D>,D>& vertSpacings = *(mesh.VertSpacings);
  BareField<Vektor<double,D>,D>& cellSpacings = *(mesh.CellSpacings);

  // Construct CenteredFieldLayout's using this for Vert and Cell centering:
  e_dim_tag edt[D];
  for (d=0; d<D; d++) 
      edt[d] = PARALLEL;
  CenteredFieldLayout<D,M,Cell> cl(mesh, edt, vnodes);
  CenteredFieldLayout<D,M,Vert> vl(mesh, edt, vnodes);

  // Use 1 guard layer in all Field's:
  GuardCellSizes<D> gc(1);



  // Vectors:
  BConds<Vektor<double,D>,D,M,Vert> vvbc;
  BConds<Vektor<double,D>,D,M,Cell> vcbc;

  // Scalars:
  BConds<double,D,M,Cell> scbc;

  // Symmetric tensors:
  BConds<SymTenzor<double,D>,D,M,Cell> stcbc;

  // Tensors:
  BConds<Tenzor<double,D>,D,M,Cell> tcbc;

  // Use linear negative reflecting conditions:
  for (int face=0; face<2*D; face++) {
    vvbc[face]  = new NegReflectFace<Vektor<double,D>,D,M,Vert>(face);
    vcbc[face]  = new NegReflectFace<Vektor<double,D>,D,M,Cell>(face);
    scbc[face]  = new NegReflectFace<double,D,M,Cell>(face);
    stcbc[face] = new NegReflectFace<SymTenzor<double,D>,D,M,Cell>(face);
    tcbc[face] =  new NegReflectFace<Tenzor<double,D>,D,M,Cell>(face);
  }

  // Now use all this to construct some Field's:
  Field<Vektor<double,D>,D,M,Vert> vectorVert(mesh, vl, gc, vvbc);
  Field<Vektor<double,D>,D,M,Cell> vectorCell(mesh, cl, gc, vcbc);
  Field<SymTenzor<double,D>,D,M,Cell> symtCell(mesh, cl, gc, stcbc);
  Field<Tenzor<double,D>,D,M,Cell> tensorCell(mesh, cl, gc, tcbc);
  Field<double,D,M,Cell> scalarCell(mesh, cl, gc, scbc);

  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  // Try out Divergence Vektor/Vert -> Scalar/Cell:
  // Assign values into the vert-centered Field<Vektor>:
  assign(vectorVert, mesh.getVertexPositionField(vectorVert));
  scalarCell = Div(vectorVert, scalarCell);
  // The value should be 3.0 for all elements; test this:
  if (abs(sum(scalarCell)/totcells - 1.0*D) > roundOffError) {
    testmsg << "Div(vert position field) != const=" << D << ")" << endl;
    testmsg << "sum(scalarCell)/totcells = " 
	    << sum(scalarCell)/totcells << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Try out Gradient Scalar/Cell -> Vektor/Vert: 

  // Use mesh object and vectorVert and scalarCell Field's constructed above.
  vectorCell = mesh.getCellPositionField(vectorCell);
  vectorCell -= mesh.get_origin();
  // Assign positive-sloping linear ramp values into the cell-centered
  // Field<scalar>:
  scalarCell = 0.0;
  for (d=0; d<D; d++) scalarCell[cells] += vectorCell[cells](d);
  // Now take the gradient:
  vectorVert = Grad(scalarCell, vectorVert);
  // The value should be (1.0,1.0,1.0) for all elements one at least one
  // removed from the last-physical-layer elements. Last-physical-layer
  // elements will be different because the BC available in IPPL don't really
  // do the kind of linear extrapolation appropriate for the needs here:
  Vektor<double,D> unit; for (d=0; d<D; d++) unit[d] = 1.0;
  Vektor<double,D> sumVectorVert;
  // Use temporary, smaller BareField as a reduced-by-two vector Field to hold
  // only the boundary-exclusive elements (needed because of limitations of
  // IPPL reductions ops):
  NDIndex<D> bev;
  for (d=0; d<D; d++) bev[d] = Index(1,nverts[d]-2,1);
  FieldLayout<D> templayout(bev);
  BareField<Vektor<double,D>,D> temp(templayout);
  temp[bev] = vectorVert[bev];
  sumVectorVert = sum(temp);
  unsigned totred=1; for (d=0; d<D; d++) totred *= nverts[d] - 2;
  sumVectorVert /= totred;
  Vektor<double,D> diffVectorVert;
  diffVectorVert = sumVectorVert - unit;
  double magDiffVectorVert = 0.0;
  for (d=0; d<D; d++) magDiffVectorVert += diffVectorVert(d)*diffVectorVert(d);
  magDiffVectorVert = sqrt(magDiffVectorVert);
  if (abs(magDiffVectorVert) > roundOffError) {
    testmsg << "Grad(cell position field) != const=(1.0,1.0,....))" << endl;
    testmsg << "sum(vectorVert)/totverts = " 
 	    << sumVectorVert << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Try out Gradient Scalar/Cell -> Vektor/Cell: 

  // Use mesh object and vectorVert and scalarCell Field's constructed above.
  vectorCell = mesh.getCellPositionField(vectorCell);
  vectorCell -= mesh.get_origin();
  // Assign positive-sloping linear ramp values into the cell-centered
  // Field<scalar>:
  scalarCell = 0.0;
  for (d=0; d<D; d++) scalarCell[cells] += vectorCell[cells](d);
  // Now take the gradient:
  vectorCell = Grad(scalarCell, vectorCell);
  // The value should be (1.0,1.0,1.0) for all elements one at least one
  // removed from the last-physical-layer elements. Last-physical-layer
  // elements will be different because the BC available in IPPL don't really
  // do the kind of linear extrapolation appropriate for the needs here:
  for (d=0; d<D; d++) unit[d] = 1.0;
  Vektor<double,D> sumVectorCell;
  // Use temporary, smaller BareField as a reduced-by-two vector Field to hold
  // only the boundary-exclusive elements (needed because of limitations of
  // IPPL reductions ops):
  NDIndex<D> bec;
  for (d=0; d<D; d++) bec[d] = Index(1,ncells[d]-2,1);
  FieldLayout<D> templayout2(bec);
  BareField<Vektor<double,D>,D> temp2(templayout);
  temp2[bec] = vectorCell[bec];
  sumVectorCell = sum(temp2);
  unsigned totredc=1; for (d=0; d<D; d++) totredc *= ncells[d] - 2;
  sumVectorCell /= totredc;
  Vektor<double,D> diffVectorCell;
  diffVectorCell = sumVectorCell - unit;
  double magDiffVectorCell = 0.0;
  for (d=0; d<D; d++) magDiffVectorCell += diffVectorCell(d)*diffVectorCell(d);
  magDiffVectorCell = sqrt(magDiffVectorCell);
  if (abs(magDiffVectorCell) > roundOffError) {
    testmsg << "Grad(cell position field) != const=(1.0,1.0,....))" << endl;
    testmsg << "sum(vectorCell)/totcells = " 
 	    << sumVectorCell << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  // Try out Divergence SymTenzor/Cell -> Vektor/Vert:

  // Use CenteredFieldLayout's from above object to construct SymTenzor Field:
  // Assign values into the cell-centered Field<SymTenzor>; use values from
  // cell-centered scalar Field scalarCell set up above:
  SymTenzor<double,D> unitSymTenzor = 1.0;
  symtCell = unitSymTenzor*scalarCell;
  // Now take the divergence:
  vectorVert = Div(symtCell, vectorVert);
  // The value should be (D,D,D,....) for all elements; test this:
  // Use temporary, smaller BareField as a reduced-by-two symtensor Field to
  // hold only the boundary-exclusive elements (needed because of limitations
  // of IPPL reductions ops):
  temp[bev] = vectorVert[bev];
  sumVectorVert = sum(temp);
  sumVectorVert /= totred;
  Vektor<double,D> deesVector; for (d=0; d<D; d++) deesVector(d) = 1.0*D;
  diffVectorVert = sumVectorVert - deesVector;
  magDiffVectorVert = 0.0;
  for (d=0; d<D; d++) magDiffVectorVert += diffVectorVert(d)*diffVectorVert(d);
  magDiffVectorVert = sqrt(magDiffVectorVert);
  if (abs(magDiffVectorVert) > roundOffError) {
    testmsg << "Div(cell position symtensor field) != const=(D,D,....))" 
	    << endl;
    testmsg << "sum(vectorVert)/totverts = " << sumVectorVert << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // Try out Gradient Vektor/Vert -> Tenzor/Cell: 

  // Set up input values in Vektor/Vert field:
  vectorVert = mesh.getVertexPositionField(vectorVert);
  // Now take the gradient:
  tensorCell = Grad(vectorVert, tensorCell);
  // Since this is the gradient of the position vector (x*x_hat + y* y_hat +
  // z*z_hat), the result should be the identity tensor (NRL Plasma Formulary
  // Vector Identities section):
  Tenzor<double,D> identityTensor = 0.0;
  for (d=0; d<D; d++) identityTensor(d,d) = 1.0;
  Tenzor<double,D> sumTensorCell = sum(tensorCell);
  sumTensorCell /= totcells;
  Tenzor<double,D> diffTensorCell;
  diffTensorCell = sumTensorCell - identityTensor;
  double magDiffTensorCell = 0.0;
  for (d=0; d<D; d++) {
    for (int d2=0; d2<D; d2++) {
      magDiffTensorCell += diffTensorCell(d,d2)*diffTensorCell(d,d2);
    }
  }
  if (abs(magDiffTensorCell) > roundOffError) {
    testmsg << "magDiffTensorCell = " << magDiffTensorCell << endl;
    testmsg << "Grad(vert position vector field) != identity tensor)" << endl;
    testmsg << "diffTensorCell = " << diffTensorCell << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------

  /* THIS TEST DOES NOT COMPILE

  // Test Average() functions:
  // Scalar Field, Cell-Centered
  Field<double,D,Cartesian<D>,Cell> C(mesh, cl, gc);
  C = 1.0;
  // Scalar weight Field, Cell-Centered
  Field<double,D,Cartesian<D>,Cell> wC(mesh, cl, gc);
  wC = 2.0;
  // Scalar Field, Vert-Centered
  Field<double,D,Cartesian<D>,Vert> V(mesh, vl, gc);
  V = 1.0;
  // Scalar weight Field, Vert-Centered
  Field<double,D,Cartesian<D>,Vert> wV(mesh, vl, gc);
  wV = 2.0;
  // Field's to hold weighted averages:
  Field<double,D,Cartesian<D>,Cell> avgToC(mesh, cl, gc);
  Field<double,D,Cartesian<D>,Vert> avgToV(mesh, vl, gc);

  assign(avgToV, Average(C, wC, avgToV));
  assign(avgToC, Average(V, wV, avgToC));

  // Weighted average from Cell to Vert:
  // ada  does not work assign(avgToV, Average(C, wC, avgToV));
  // Weighted average from Vert to Cell:
  // ada dones not work assign(avgToC, Average(V, wV, avgToC));
  // Check results:
  if (sum(avgToV) != totverts) {
    testmsg << "avgToV values wrong" << endl;
    testmsg << "sum(avgToV) = " << sum(avgToV) << " ; totverts = " << totverts
	    << endl;
    // passed = false;
  }
  if (sum(avgToC) != totcells) {
    testmsg << "avgToC values wrong" << endl;
    testmsg << "sum(avgToC) = " << sum(avgToC) << " ; totcells = " << totcells
	    << endl;
    // passed = false;
  }

  */

  //---------------------------------------------------------------------------

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if DIM == 3 // This block-o-tests only compiles in 3D

  //---------------------------------------------------------------------------
  // Some accessor function tests:
  double theVolume, theVolume2, theVolume3;
  NDIndex<D> ndi;
  ndi[0] = Index(2,2,1);
  ndi[1] = Index(2,2,1);
  ndi[2] = Index(2,2,1);
  theVolume = mesh.getCellVolume(ndi);
  ndi[0] = Index(0,2,1);
  ndi[1] = Index(0,2,1);
  ndi[2] = Index(0,2,1);
  theVolume2 = mesh.getCellRangeVolume(ndi);
  if (theVolume2 != (6.0*12.0*18.0)) {
    testmsg << "volume of cells [0:2][0:2][0:2] != 1296" << endl;
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
  Field<double,D,Cartesian<D>,Cell> theVolumes(mesh, cl);
  mesh.getCellVolumeField(theVolumes);
  if ((sum(theVolumes)/totcells) != theVolume) {
    testmsg << "(sum(theVolumes)/totcells) != cell vol" << endl;
    testmsg << "(sum(theVolumes)/totcells) = " 
	    << (sum(theVolumes)/totcells) << endl;
    testmsg << "cell vol = " << theVolume << endl;
  }
  //---------------------------------------------------------------------------
  Vektor<double,D> v;
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
  Vektor<double,D> v1;
  v1 = mesh.getVertexPosition(ndi);
  v(0) = 2.0; v(1) = 4.0; v(2) = 12.0; // Correct value
  if (v1 != v) {
    testmsg << "VERT POSITION OF" << ndi << " = " << v1 << " != " << v << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------
  CenteredFieldLayout<D,Cartesian<D>,Vert> 
    clVert(mesh);
  Field<Vektor<double,D>,D,Cartesian<D>,Vert> 
    thePositions(clVert);
  mesh.getVertexPositionField(thePositions);
  //---------------------------------------------------------------------------
  v = mesh.getDeltaVertex(ndi);
  Vektor<double,D> vcorrect;
  vcorrect(0) = 2.0; vcorrect(1) = 4.0; vcorrect(2) = 9.0;
  if (v != vcorrect) {
    testmsg << "DELTA-VERTEX OF" << ndi << " = " << v 
	    << " != " << vcorrect << endl;
    passed = false;
  }
  //---------------------------------------------------------------------------

#endif // DIM == 3

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: TestCartesian.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: TestCartesian.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
