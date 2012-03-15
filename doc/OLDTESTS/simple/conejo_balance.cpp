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
#include "Utility/FieldDebug.h"
#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "FieldLayout/ConejoBalancer.h"
#include "Field/Field.h"
#include "Message/Communicate.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <iostream>
using namespace std;
#else
#include <iostream.h>
#endif


int main(int argc, char *argv[])
{
  IpplInfo ippl(argc,argv);
  int procs = ippl.getNodes();
  int vnodes = procs*4;

  const int N=100;
  Index I(N), J(N);
  Inform out;
  out.setPrintNode();

  FieldLayout<2> layout(I,J,PARALLEL,PARALLEL,vnodes);
  out << layout << endl ;

  Field<double,2> weight1(layout);
  Field<double,2> weight2(layout);
  weight1[I][J] = I*J/(10000.0);
  weight2[I][J] = (N-1-I)*(N-1-J)/(10000.0);

  double x;

  x = 0;
  for (Field<double,2>::iterator fp=weight1.begin(); fp!=weight1.end(); ++fp)
    x += *fp;
  out <<"Starting weight: " << x << endl;
  x = 0;
  for (Field<double,2>::iterator fp=weight2.begin(); fp!=weight2.end(); ++fp)
    x += *fp;
  out <<"Starting weight: " << x << endl;

  ConejoBalancer balancer;
  out <<"Add Material" << endl;
  balancer.addMaterial(weight1,false);
  balancer.addMaterial(weight2,false);
  out <<"Redistribute" << endl;
  balancer.redistribute(layout);
  out << layout << endl ;

  x=0;
  for (Field<double,2>::iterator fp=weight1.begin(); fp!=weight1.end(); ++fp)
    x += *fp;
  out <<"Ending weight: " << x << endl;
  x=0;
  for (Field<double,2>::iterator fp=weight2.begin(); fp!=weight2.end(); ++fp)
    x += *fp;
  out <<"Ending weight: " << x << endl;

  Ippl::Comm->barrier();
  return 0;
}

/***************************************************************************
 * $RCSfile: conejo_balance.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:40 $
 * IPPL_VERSION_ID: $Id: conejo_balance.cpp,v 1.1.1.1 2000/11/30 21:05:40 adelmann Exp $ 
 ***************************************************************************/
