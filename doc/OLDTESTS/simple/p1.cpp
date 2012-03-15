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

  Index I(10);
  const unsigned Dim=1;

  // Set initial boundary conditions.
  BConds<double,Dim> bc;
  bc[0] = new PosReflectFace<double,Dim>(1);
  bc[1] = new ZeroFace<double,Dim>(0);

  FieldLayout<Dim> layout(I);
  *logfile << "Built layout" << endl;
  Ippl::Comm->barrier();
  Field<double,Dim> A(layout,GuardCellSizes<Dim>(1),bc), B(layout);
  *logfile << "Built fields" << endl;
  Ippl::Comm->barrier();

  A[I] = I;
  *logfile << "Initialized A" << endl;
  Ippl::Comm->barrier();
  Field<double,Dim>::iterator p;
  *logfile << "A = ";
  for (p=A.begin(); p!=A.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  B = 0.0;
  *logfile << "Initialized B" << endl;
  Ippl::Comm->barrier();
  *logfile << "B = ";
  for (p=B.begin(); p!=B.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  A = 0.0;
  *logfile << "Initialized A" << endl;
  Ippl::Comm->barrier();
  *logfile << "A = ";
  for (p=A.begin(); p!=A.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  B[I] = I;
  *logfile << "Initialized B" << endl;
  Ippl::Comm->barrier();
  *logfile << "B = ";
  for (p=B.begin(); p!=B.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();

  A[I] = I;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  A[I] = A[I] + 0.1*A[I-1] + 10.0*A[I+1] ;
#else
  assign(A[I] , A[I] + 0.1*A[I-1] + 10.0*A[I+1]);
#endif
  *logfile << "Calculated A" << endl;
  Ippl::Comm->barrier();
  *logfile << "A = ";
  for (p=A.begin(); p!=A.end(); ++p) 
    *logfile << " " << *p;
  *logfile << endl;
  Ippl::Comm->barrier();
}
/***************************************************************************
 * $RCSfile: p1.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:14 $
 * IPPL_VERSION_ID: $Id: p1.cpp,v 1.1.1.1 2000/11/30 21:06:14 adelmann Exp $ 
 ***************************************************************************/
