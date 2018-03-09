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

#include "Ippl.h"
#include <vector>
#define PRINT_DETAILS

//////////////////////////////////////////////////////////////////////

const unsigned Dim=2;
const int N=10;

int main(int argc, char *argv[])
{
  Inform out;
  Ippl ippl(argc, argv);

  Index I(N), J(N);
  FieldLayout<Dim> layout1(I,J,PARALLEL,SERIAL);
  Field<int,Dim> A1(layout1);
  Field<int,Dim> A2(layout1);

  assign(A1[I][J], I + 10*J);
  assign(A2[I][J], J + 10*I);

  // an alternative way to generate the list of local vnode domains,
  // directly from a FieldLayout
  vector< NDIndex<Dim> > newDomains2;

  FieldLayout<Dim>       layout2(I,J,SERIAL,PARALLEL);

  FieldLayout<Dim>::iterator_iv lvnode, endvnode=layout2.end_iv();
  for (lvnode = layout2.begin_iv(); lvnode != endvnode; ++lvnode) 
    newDomains2.push_back(NDIndex<Dim>((*(*lvnode).second).getDomain()));

  layout1.Repartition(newDomains2.begin(),newDomains2.end());

  out << "My portion of Field A1 = " << A1 << endl;
  out << "My portion of Field A2 = " << A2 << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: repartition.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: repartition.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/





