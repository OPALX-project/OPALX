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
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/LField.h"

// #define PRINT_DETAILS

#ifdef PRINT_DETAILS
// define helper print function
void print_field(Inform& out, Field<int,2>& f)
{
  for (Field<int,2>::iterator_if lf_i=f.begin_if(); lf_i!=f.end_if(); ++lf_i)
    {
      LField<int,2> &lf = *(*lf_i).second;
      for (LField<int,2>::iterator p = lf.begin(); p!=lf.end(); ++p)
	out << " " << *p ;
      out << endl;
    }
  out << endl;
}
#endif


int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform out(argv[0]);
#ifdef PRINT_DETAILS
  out.setPrintNode();
#endif

  const unsigned Dim=2;
  const int N=10;
  Index I(N);
  Index J(N);
  FieldLayout<Dim> layout1(I,J,PARALLEL,SERIAL,4);
  FieldLayout<Dim> layout2(I,J,SERIAL,PARALLEL,8);
  Field<int,Dim> A1(layout1),A2(layout2);

  A1 = 0;
  A1[I][J] += I + 10*J;
  A2 = A1;
#ifdef PRINT_DETAILS
  print_field(out,A1);
  print_field(out,A2);
#endif

  A2[I][J] -= I + 10*J;
  A2 *= A2;
  int s = sum(A2);
  if ( s==0 )
    out << "PASSED" << endl;
  else
    out << "FAILED" << endl;
  
#ifdef PRINT_DETAILS
  print_field(out,A2);
#endif

  return 0;
}

/***************************************************************************
 * $RCSfile: balance.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:26 $
 * IPPL_VERSION_ID: $Id: balance.cpp,v 1.1.1.1 2000/11/30 21:05:26 adelmann Exp $ 
 ***************************************************************************/





