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

#include <iostream>
#include "Ippl.h"
 
template < unsigned Dim > 
void initVectorVert(Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>&);
template < unsigned Dim > 
void initScalarVert(Field<double,Dim,UniformCartesian<Dim>,Vert>&);

const unsigned Dim = 3;
int size = 8;

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);

  NDIndex<Dim> cellDomain, vertDomain;
  for (int i = 0 ; i < Dim ; i++ ) {
    cellDomain[i] = Index(size);
    vertDomain[i] = Index(size+1);
  }

  // Some boundary condition objects:
  BConds<double,Dim,UniformCartesian<Dim>,Vert> bcScalarVert;
  BConds<double,Dim,UniformCartesian<Dim>,Cell> bcScalarCell;
  for (int d=0; d<2*Dim; d++) {
    bcScalarVert[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Vert>(d);
    bcScalarCell[d] = new ZeroFace<double,Dim,UniformCartesian<Dim>,Cell>(d);
  }
  BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> bcVectorVert;
  BConds<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> bcVectorCell;
  for (int d=0; d<2*Dim; d++) {
    bcVectorCell[d] = 
      new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>(d);
    bcVectorVert[d] = 
      new ZeroFace<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert>(d);
  }

  // perform the vertex->cell divergenace operation 
  FieldLayout<Dim> layoutCell(cellDomain), layoutVert(vertDomain);
  Field<double,Dim,UniformCartesian<Dim>,Cell> 
    ScaC1(layoutCell, bcScalarCell, GuardCellSizes<Dim>(1));
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> 
    VecV1(layoutVert, bcVectorVert, GuardCellSizes<Dim>(1));
  
  initVectorVert(VecV1);

  VecV1.writeb("vtest1");
  VecV1.fillGuardCells();
  VecV1.writeb("vtest2");

  ScaC1 = Div(VecV1,ScaC1);
 
  ScaC1.writeb("ScaC1");

  // now perform the vertex->cell gradient operation 
  Field<double,Dim,UniformCartesian<Dim>,Vert> 
    ScaV1(layoutVert, bcScalarVert, GuardCellSizes<Dim>(1));
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> 
    VecC1(layoutCell, bcVectorCell, GuardCellSizes<Dim>(1));
  
  initScalarVert(ScaV1);

  ScaV1.writeb("stest1");
  ScaV1.fillGuardCells();
  ScaV1.writeb("stest2");

  VecC1 = Grad(ScaV1,VecC1);

  VecC1.writeb("VecC1");
  return 0;
}
//----------------------------------------------------------------------
void initVectorVert(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Vert>& V1) {
  double X[1] = { 1.0 } ;
  V1.get_mesh().set_meshSpacing(X);
  FieldLayout<1U>& layoutVert = V1.getLayout();
  const GuardCellSizes<1U>& guardCells = V1.getGuardCellSizes();
  Field<double,1U,UniformCartesian<1U>,Vert> 
    Vscalar1(layoutVert, guardCells);
  Field <Vektor<double,1U>,1U,UniformCartesian<1U>,Vert >::iterator pv;
  Field <double,1U,UniformCartesian<1U>,Vert >::iterator pf1;
  Index I = V1.getIndex(0);
  Vscalar1[I]= I;
  V1.Uncompress();
  for (pv=V1.begin(), pf1=Vscalar1.begin();
       pv!=V1.end() ; ++pv, ++pf1) {
    (*pv)[0] = *pf1;
  }
}
//----------------------------------------------------------------------
void initVectorVert(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Vert>& V1) {
  double X[2] = { 1.0, 1.0 } ;
  V1.get_mesh().set_meshSpacing(X);
  FieldLayout<2U>& layoutVert = V1.getLayout();
  const GuardCellSizes<2U>& guardCells = V1.getGuardCellSizes();
  Field<double,2U,UniformCartesian<2U>,Vert> 
    Vscalar1(layoutVert, guardCells),
    Vscalar2(layoutVert, guardCells);
  Field <Vektor<double,2U>,2U,UniformCartesian<2U>,Vert >::iterator pv;
  Field <double,2U,UniformCartesian<2U>,Vert >::iterator pf1, pf2;
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  Vscalar1[I][J]= I;
  Vscalar2[I][J]= J;
  V1.Uncompress();
  for (pv=V1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ; 
       pv!=V1.end() ; ++pv, ++pf1, ++pf2) {
    (*pv)[0] = *pf1;
    (*pv)[1] = *pf2;
  }
}
//----------------------------------------------------------------------
void initVectorVert(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Vert>& V1) {
  double X[3] = { 1.0, 1.0, 1.0 } ;
  V1.get_mesh().set_meshSpacing(X);
  FieldLayout<3U>& layoutVert = V1.getLayout();
  const GuardCellSizes<3U>& guardCells = V1.getGuardCellSizes();
  Field<double,3U,UniformCartesian<3U>,Vert> 
    Vscalar1(layoutVert, guardCells),
    Vscalar2(layoutVert, guardCells),
    Vscalar3(layoutVert, guardCells);
  Field <Vektor<double,3U>,3U,UniformCartesian<3U>,Vert >::iterator pv;
  Field <double,3U,UniformCartesian<3U>,Vert >::iterator pf1, pf2, pf3;
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  Index K = V1.getIndex(2);
  Vscalar1[I][J][K]= I;
  Vscalar2[I][J][K]= J;
  Vscalar3[I][J][K]= K;
  V1.Uncompress();
  for (pv=V1.begin(), pf1=Vscalar1.begin(), 
	 pf2=Vscalar2.begin(), pf3=Vscalar3.begin() ;  
       pv!=V1.end() ; ++pv, ++pf1, ++pf2, ++pf3) {
    (*pv)[0] = *pf1;
    (*pv)[1] = *pf2;
    (*pv)[2] = *pf3;
  }
}
//----------------------------------------------------------------------
void initScalarVert(Field<double,1U,UniformCartesian<1U>,Vert>& V1) {
  Index I = V1.getIndex(0);
  V1[I]= I ;
}
//----------------------------------------------------------------------
void initScalarVert(Field<double,2U,UniformCartesian<2U>,Vert>& V1) {
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  V1[I][J]= I ;
}
//----------------------------------------------------------------------
void initScalarVert(Field<double,3U,UniformCartesian<3U>,Vert>& V1) {
  Index I = V1.getIndex(0);
  Index J = V1.getIndex(1);
  Index K = V1.getIndex(2);
  V1[I][J][K]= I ;
}
//----------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: vc.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: vc.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
