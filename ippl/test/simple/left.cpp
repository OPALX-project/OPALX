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
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform pout(NULL,0);
  setInform(pout);
  setFormat(8,2);
  bool b1=true; bool b2=false;
  pout << "bool output key: true = " << b1 << " ; false = " << b2 << endl;

  // Sizes
  const unsigned D = 1;
  int ncells = 7, nverts = 8;
  int vnodes = 2;

  // Mesh
  NDIndex<1U> verts; 
  for (int d=0; d<D; d++) verts[d] = Index(nverts);
  typedef UniformCartesian<1U> M1;
  M1 mesh(verts);

  // Layouts
  CenteredFieldLayout<1U,M1,Cell> layoutCell(mesh,PARALLEL,vnodes);
  CenteredFieldLayout<1U,M1,Vert> layoutVert(mesh,PARALLEL,vnodes);

  // Guard cells
  GuardCellSizes<1U> gc(1);

  // Boundary conditions
  BConds<double,1U,M1,Cell> cbc;
  BConds<double,1U,M1,Vert> vbc;
  for (int face=0; face<2*1U; face++) {
    cbc[face] = new NegReflectFace<double,1U,M1,Cell>(face);
    vbc[face] = new NegReflectFace<double,1U,M1,Vert>(face);
  }

  // Fields
  Field<double,1U,M1,Cell> A(layoutCell,gc,cbc);
  Field<double,1U,M1,Vert> B(layoutVert,gc,vbc);
  Field<double,1U,M1,Cell> A0(layoutCell,gc,cbc);
  Field<double,1U,M1,Vert> B0(layoutVert,gc,vbc);

  // Initial values (duplicate in A0,B0 for FieldDebug output w/o changing
  // dirty_m of actual A and B Fields):
  Index I(ncells);
  Index Iv(nverts);
  A[I] = I;
  A0[I] = I;
  B[Iv] = Iv;
  B0[Iv] = Iv;

  // Output State of dirty flags prior to "stencil" assignment:
  pout << "!!!!!!!!!!!!! BEFORE !!!!!!!!!!!!!" << endl;
  pout << "B.isDirty() = " << B.isDirty() << " ; " 
       << "A.isDirty() = " << A.isDirty() << endl;

  // Output (copies of) initial values):
  pout << "[[[[[[[ A ]]]]]]]" << endl;
  fp1(A0);
  pout << "[[[[[[[ B ]]]]]]]" << endl;
  fp1(B0);

  // "Stencil" assignment:
  B[I + 1] = B[I + 1] + A[I + 1];

  // Output State of dirty flags after to "stencil" assignment:
  pout << "!!!!!!!!!!!!! AFTER !!!!!!!!!!!!!" << endl;
  pout << "B.isDirty() = " << B.isDirty() << " ; " 
       << "A.isDirty() = " << A.isDirty() << endl;

  // Output resulting value of B:
  pout << "[[[[[[[ B ]]]]]]]" << endl;
  fp1(B);

  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

