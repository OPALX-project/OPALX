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


#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>

#include "Ippl.h"

ostream *logfile;

int main(int argc, char *argv[] )
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);

  char name[100];
  sprintf(name,"out.%d",Ippl::Comm->myNode());
  logfile = new ofstream(name);
  *logfile << "Starting to run!" << endl;
  Ippl::Comm->barrier();

  Index I(5),J(5);
  const unsigned Dim=2;

  // Set initial boundary conditions.
  BConds<double,Dim> bc;
  bc[0] = new NegReflectFace<double,Dim>(0);
  bc[1] = new NegReflectFace<double,Dim>(1);
  bc[2] = new NegReflectFace<double,Dim>(2);
  bc[3] = new NegReflectFace<double,Dim>(3);

  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1),bc), B(layout);

  Field<double,Dim>::iterator p;

#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I][J] = (I+1)+10*(J+1);
#else
  assign(A[I][J], (I+1)+10*(J+1));
#endif

  *logfile << "Calculated A" << endl;
  Ippl::Comm->barrier();
  *logfile << "A = ";
  for (p=A.begin(); p!=A.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  B[I][J] = A[I+1][J+1] ;

  *logfile << "Calculated B" << endl;
  Ippl::Comm->barrier();
  *logfile << "B = ";
  for (p=B.begin(); p!=B.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  *logfile << "sum(A)  = " << sum(A) << endl;
  *logfile << "prod(A) = " << prod(A) << endl;
  *logfile << "max(A)  = " << max(A) << endl;
  *logfile << "min(A)  = " << min(A) << endl;
  Ippl::Comm->barrier();
}
/***************************************************************************
 * $RCSfile: p2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:15 $
 * IPPL_VERSION_ID: $Id: p2.cpp,v 1.1.1.1 2000/11/30 21:06:15 adelmann Exp $ 
 ***************************************************************************/
