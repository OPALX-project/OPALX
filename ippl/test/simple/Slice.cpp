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
#include "Utility/Inform.h"
#include "Utility/FieldView.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <iostream>
using namespace std;
#else
#include <iostream>
#endif

#include <math.h>


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  int in;

  cout << "input the number of vnodes" << endl;
  int vnodes;
  cin >> vnodes; 
  cout << "input the number of size" << endl;
  int size;
  cin >> size; 
  const unsigned Dim=3;
  Index I(size), J(size), K(size);
  FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
  Field<double,Dim> X(layout), Y(layout), Z(layout);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));
  cout << " Layout is: "<< endl ;
  testmsg << layout << endl;

  double pi = 2.0*(4.0*atan(1.0))/size;

#ifdef IPPL_USE_MEMBER_TEMPLATES
  X[I][J][K]  = pi * I;
  Y[I][J][K]  = pi * J;
  Z[I][J][K]  = pi * K;
  A = sin(X)*cos(Y)*Z;
#else
  X[I][J][K] << pi * I;
  Y[I][J][K] << pi * J;
  Z[I][J][K] << pi * K;
  A << sin(X)*cos(Y)*Z;
#endif

  A.write("atest");

  FieldView<double,Dim> plotX(0,A);
  int ipplot;
  for(ipplot = 0 ; ipplot<size; ipplot++) {
    plotX.view(ipplot);
  }
  cout << "enter an integer to continue" << endl;
  cin >> in;
  
  FieldView<double,Dim> plotY(1,A);
  for(ipplot = 0 ; ipplot<size; ipplot++) {
    plotY.view(ipplot);
  }
  cout << "enter an integer to continue" << endl;
  cin >> in;
  
  FieldView<double,Dim> plotZ(2,A);
  for(ipplot = 0 ; ipplot<size; ipplot++) {
    plotZ.view(ipplot);
  }
  cout << "enter an integer to continue" << endl;
  cin >> in;

  return 0;
}

/***************************************************************************
 * $RCSfile: Slice.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: Slice.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
