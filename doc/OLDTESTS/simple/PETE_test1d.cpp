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

// include files
#include "Utility/IpplInfo.h"
#include "Index/Index.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Index/Index.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/GuardCellSizes.h"
#include "Meshes/UniformCartesian.h"


int main(int argc, char *argv[]) {

  Ippl ippl(argc,argv); // Standard Ippl startup

  // set some constants
  const unsigned Dim=1;
  const int sizeX      = 50;
  const int centerX    = 25;

  // build the needed Ippl objects
  Index I(sizeX);
  Index I1(sizeX+1);
  UniformCartesian<Dim> myMesh(I1);
  CenteredFieldLayout<Dim,UniformCartesian<Dim>,Cell>
    cell_layout(myMesh,PARALLEL);
  Field<double,Dim> A(myMesh,cell_layout,GuardCellSizes<Dim>(1));
  Field<double,Dim> B(myMesh,cell_layout);
  Field<bool,Dim> mask(myMesh,cell_layout);

  // initialize Fields
  A = 0.0;
  A[centerX] = 100;
  assign( mask[I], lt(I,(double) centerX) );
  const double fact1 = 1.0/6.0;
  const double fact2 = 1.0/8.0;

  assign( B[I], where( mask[I],
                       fact1 * (A[I] + A[I] + A[I] + A[I] + A[I+1] + A[I-1]),
                       fact2 * (A[I+1] + A[I-1] + A[I+1] + A[I-1] +
		                A[I+1] + A[I-1] + A[I+1] + A[I-1]) ) );

  return 0;
}

/***************************************************************************
 * $RCSfile: PETE_test1d.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:04:13 $
 * IPPL_VERSION_ID: $Id: PETE_test1d.cpp,v 1.1.1.1 2000/11/30 21:04:13 adelmann Exp $
 ***************************************************************************/
