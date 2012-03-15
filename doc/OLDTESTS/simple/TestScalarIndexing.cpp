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

// TestScalarIndexing.cpp , Tim Williams 10/23/1996
// This tests indexing Field objects with scalar index values.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"

// set dimensionality and problem size
const unsigned Dim=3;
const unsigned nx = 4, ny = 4, nz = 4;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Index I(nx), J(ny), K(nz);
  FieldLayout<Dim> layout(I,J,K);

  // Instantiate and initialize scalar, vector, tensor fields:
  Field<double,Dim> scalarFld(layout);
  double scalar = 1.0;
  scalarFld << scalar;
  Field<Vektor<double,Dim>,Dim> vectorFld(layout);
  Vektor<double, Dim> vector(1.0,2.0,3.0);
  vectorFld << vector;
  Field<Tenzor<double,Dim>,Dim,UniformCartesian<Dim> > tensorFld(layout);
  Tenzor<double, Dim> tensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
  tensorFld << tensor;
  Field<SymTenzor<double,Dim>,Dim,UniformCartesian<Dim> > symTensorFld(layout);
  SymTenzor<double, Dim> symTensor(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  symTensorFld << symTensor;

  // Now try the scalar indexing:
  double scalar1 = 0.0;
  Vektor<double, Dim> vector1(0.0, 0.0, 0.0);
  Tenzor<double, Dim> tensor1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  SymTenzor<double, Dim> symTensor1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  scalar1 = scalarFld[1][1][1].get();
  vector1 = vectorFld[1][1][1].get();
  tensor1 = tensorFld[1][1][1].get();
  symTensor1 = symTensorFld[1][1][1].get();

  bool works = true;
  if (scalar1 != scalar) works = false;
  if (vector1 != vector) works = false;
  if (tensor1 != tensor) works = false;
  if (symTensor1 != symTensor) works = false;

  if (works) {
    testmsg << "PASSED!" << endl;
  }
  else {
    testmsg << "FAILED!" << endl;
  }
  return 0;
}

/***************************************************************************
 * $RCSfile: TestScalarIndexing.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:20 $
 * IPPL_VERSION_ID: $Id: TestScalarIndexing.cpp,v 1.1.1.1 2000/11/30 21:05:20 adelmann Exp $ 
 ***************************************************************************/
