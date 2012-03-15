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

#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Utility/PAssert.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  const int Dim=2;
  int n=5;
  Index I(n),J(n);
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout);
  int i,j;

  for (i=0; i<n; ++i) {
    for (j=0; j<n; ++j) { 
      A[i][j] = i+10*j;
    }
  }
  for (j=n-1; j>=0; --j) {
    for (i=0; i<n; ++i) {
      PInsist(A[i][j].get() == i+10*j, "Failure in test code single.cpp!!");
    }
  }
  testmsg << "PASSED!" << endl;
}
/***************************************************************************
 * $RCSfile: single.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:24 $
 * IPPL_VERSION_ID: $Id: single.cpp,v 1.1.1.1 2000/11/30 21:06:24 adelmann Exp $ 
 ***************************************************************************/
