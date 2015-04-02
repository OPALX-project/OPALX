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
void 
initVektorCell(Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>&);
template < unsigned Dim > 
void 
initTenzorCell(Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim>,Cell>&);

const unsigned Dim = 2;
int size = 4;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  
  NDIndex<Dim> cellDomain, vertDomain;
  for (int i = 0 ; i < Dim ; i++ ) {
    cellDomain[i] = Index(size);
    vertDomain[i] = Index(size+1);
  }

  FieldLayout<Dim> layoutCell(cellDomain), layoutVert(vertDomain);

  // perform the cell->vertex divergence operation 
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> 
    VecV1(layoutCell);
  Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> 
    TzrC1(layoutCell, GuardCellSizes<Dim>(1));
  
  initTenzorCell(TzrC1);
  
  TzrC1.writeb("ttest1");
  TzrC1.fillGuardCells();
  TzrC1.writeb("ttest2");
  
  Div(TzrC1,VecV1);

  VecV1.writeb("VecV1");

  // now perform the cell->vertex gradient operation 
  Field<Vektor<double,Dim>,Dim,UniformCartesian<Dim>,Cell> 
    VecC1(layoutCell, GuardCellSizes<Dim>(1));
  Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim>,Vert> 
    TzrV1(layoutVert);
  
  initVektorCell(VecC1);

  VecC1.writeb("vtest1");
  VecC1.fillGuardCells();
  VecC1.writeb("vtest2");

  Grad(VecC1,TzrV1);

  TzrV1.writeb("TrzV1");

  return 0;
}
//----------------------------------------------------------------------
void 
initTenzorCell(Field<Tenzor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1) {
  double X[1] = { 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<1U>& layoutCell = C1.getLayout();
  const GuardCellSizes<1U>& guardCells = C1.getGuardCellSizes();
  Field<double,1U,UniformCartesian<1U>,Cell> 
    Vscalar1(layoutCell, guardCells);
  Field <Tenzor<double,1U>,1U,UniformCartesian<1U>,Cell >::iterator pv;
  Field <double,1U,UniformCartesian<1U>,Cell >::iterator pf1;
  Index I = C1.getIndex(0);
  Vscalar1[I]= I;
  C1.Uncompress();
  for (pv=C1.begin(), pf1=Vscalar1.begin() ;
       pv!=C1.end() ; ++pv, ++pf1) {
    (*pv)(0,0) = *pf1;
  }
}
//----------------------------------------------------------------------
void 
initTenzorCell(Field<Tenzor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1) {
  double X[2] = { 1.0, 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<2U>& layoutCell = C1.getLayout();
  const GuardCellSizes<2U>& guardCells = C1.getGuardCellSizes();
  Field<double,2U,UniformCartesian<2U>,Cell> 
    Vscalar1(layoutCell, guardCells),
    Vscalar2(layoutCell, guardCells);
  Field <Tenzor<double,2U>,2U,UniformCartesian<2U>,Cell >::iterator pv;
  Field <double,2U,UniformCartesian<2U>,Cell >::iterator pf1;
  Index I = C1.getIndex(0);
  Index J = C1.getIndex(1);
  Vscalar1[I][J]= I;
  C1.Uncompress();
  for (pv=C1.begin(), pf1=Vscalar1.begin() ; pv!=C1.end() ; ++pv, ++pf1) {
    (*pv)(0,0) = *pf1;
    (*pv)(0,1) = 0.0;
    (*pv)(1,0) = 0.0;
    (*pv)(1,1) = *pf1;
  }
}
//----------------------------------------------------------------------
void 
initTenzorCell(Field<Tenzor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1) {
  double X[3] = { 1.0, 1.0, 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<3U>& layoutCell = C1.getLayout();
  const GuardCellSizes<3U>& guardCells = C1.getGuardCellSizes();
  Field<double,3U,UniformCartesian<3U>,Cell> 
    Vscalar1(layoutCell, guardCells),
    Vscalar2(layoutCell, guardCells),
    Vscalar3(layoutCell, guardCells);
  Field <Tenzor<double,3U>,3U,UniformCartesian<3U>,Cell >::iterator pv;
  Field <double,3U,UniformCartesian<3U>,Cell >::iterator pf1, pf2, pf3;
  Index I = C1.getIndex(0);
  Index J = C1.getIndex(1);
  Index K = C1.getIndex(2);
  Vscalar1[I][J][K]= I;
  C1.Uncompress();
  for (pv=C1.begin(), pf1=Vscalar1.begin() ; pv!=C1.end() ; ++pv, ++pf1) {
    (*pv)(0,0) = *pf1;
    (*pv)(0,1) = 0.0;
    (*pv)(0,2) = 0.0;
    (*pv)(1,0) = 0.0;
    (*pv)(1,1) = *pf1;
    (*pv)(1,2) = 0.0;
    (*pv)(2,0) = 0.0;
    (*pv)(2,1) = 0.0;
    (*pv)(2,2) = *pf1;
  }
}
//----------------------------------------------------------------------
void 
initVektorCell(Field<Vektor<double,1U>,1U,UniformCartesian<1U>,Cell>& C1) {
  double X[1] = { 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<1U>& layoutCell = C1.getLayout();
  const GuardCellSizes<1U>& guardCells = C1.getGuardCellSizes();
  Field<double,1U,UniformCartesian<1U>,Cell> 
    Vscalar1(layoutCell, guardCells);
  Field <Vektor<double,1U>,1U,UniformCartesian<1U>,Cell >::iterator pv;
  Field <double,1U,UniformCartesian<1U>,Cell >::iterator pf1;
  Index I = C1.getIndex(0);
  Vscalar1[I]= I;
  for (pv=C1.begin(), pf1=Vscalar1.begin() ;
       pv!=C1.end() ; ++pv, ++pf1) {
    (*pv)[0] = *pf1;
  }
}
//----------------------------------------------------------------------
void 
initVektorCell(Field<Vektor<double,2U>,2U,UniformCartesian<2U>,Cell>& C1) {
  double X[2] = { 1.0, 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<2U>& layoutCell = C1.getLayout();
  const GuardCellSizes<2U>& guardCells = C1.getGuardCellSizes();
  Field<double,2U,UniformCartesian<2U>,Cell> 
    Vscalar1(layoutCell, guardCells),
    Vscalar2(layoutCell, guardCells);
  Field <Vektor<double,2U>,2U,UniformCartesian<2U>,Cell >::iterator pv;
  Field <double,2U,UniformCartesian<2U>,Cell >::iterator pf1, pf2;
  Index I = C1.getIndex(0);
  Index J = C1.getIndex(1);
  Vscalar1[I][J]= I;
  Vscalar2[I][J]= J;
  C1.Uncompress();
  Vscalar1.Uncompress();
  Vscalar2.Uncompress();
  for (pv=C1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ; 
       pv!=C1.end() ; ++pv, ++pf1, ++pf2) {
    (*pv)[0] = *pf1;
    (*pv)[1] = *pf2;
  }
}
//----------------------------------------------------------------------
void 
initVektorCell(Field<Vektor<double,3U>,3U,UniformCartesian<3U>,Cell>& C1) {
  double X[3] = { 1.0, 1.0, 1.0 } ;
  C1.get_mesh().set_meshSpacing(X);
  FieldLayout<3U>& layoutCell = C1.getLayout();
  const GuardCellSizes<3U>& guardCells = C1.getGuardCellSizes();
  Field<double,3U,UniformCartesian<3U>,Cell> 
    Vscalar1(layoutCell, guardCells),
    Vscalar2(layoutCell, guardCells),
    Vscalar3(layoutCell, guardCells);
  Field <Vektor<double,3U>,3U,UniformCartesian<3U>,Cell >::iterator pv;
  Field <double,3U,UniformCartesian<3U>,Cell >::iterator pf1, pf2, pf3;
  Index I = C1.getIndex(0);
  Index J = C1.getIndex(1);
  Index K = C1.getIndex(2);
  Vscalar1[I][J][K]= I;
  Vscalar2[I][J][K]= J;
  Vscalar3[I][J][K]= K;
  C1.Uncompress();
  Vscalar1.Uncompress();
  Vscalar2.Uncompress();
  Vscalar3.Uncompress();
  for (pv=C1.begin(), pf1=Vscalar1.begin(), pf2=Vscalar2.begin() ,
       pf3=Vscalar3.begin() ; 
       pv!=C1.end() ; ++pv, ++pf1, ++pf2, ++pf3) {
    (*pv)[0] = *pf1;
    (*pv)[1] = *pf2;
    (*pv)[2] = *pf3;
  }
}
//----------------------------------------------------------------------
/***************************************************************************
 * $RCSfile: tv.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: tv.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
