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

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  const unsigned Dim=2;
  Index I(5);
  Index J(5);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL);
  Field<double,Dim> A(layout);

  Vektor<double,3> r[10];
  Vektor<double,3> rtmp[10];

  for (int i=0;i<10;i++)
      for(int d=0; d<3; d++)
          rtmp[i](d)=i*d;
  
  for (int i=0;i<10;i++) {      
      for(int d=0; d<3; d++)
          testmsg << rtmp[i](d) << "\t";
      testmsg << endl;
  }

  reduce(rtmp,rtmp+10,r,OpAddAssign());
  
  for (int i=0;i<10;i++) {      
      for(int d=0; d<3; d++)
          testmsg << r[i](d) << "\t";
      testmsg << endl;
  }


  A[I][J] = I+2*J+1;

  testmsg << "min(A) = " << min(A) << endl;
  testmsg << "max(A) = " << max(A) << endl;
  testmsg << "sum(A) = " << sum(A) << endl;
  testmsg << "prod(A)= " << prod(A) << endl;

  double minv, maxv;
  minmax(A, minv, maxv);
  testmsg << "minmax(A) = (" << minv << ", " << maxv << ")" << endl;

  testmsg << "any(A,1)= " << any(A,1.0) << endl;
  testmsg << "any(A,10)= " << any(A,10.0) << endl;
  testmsg << "any(A,1,OpGT())= " << any(A,1.0,OpGT()) << endl;
  testmsg << "any(A,10,OpLT())= " << any(A,10.0,OpLT()) << endl;

  NDIndex<Dim> Loc;
  double m;
  m = min(A,Loc);
  testmsg << "min(A) = " << m << " at " << Loc << endl;
  m = max(A,Loc);
  testmsg << "max(A) = " << m << " at " << Loc << endl;

  Field<Vektor<double,Dim>,Dim> B(layout);
  B.Uncompress();
  Field<Vektor<double,Dim>,Dim>::iterator Bit, Bend = B.end();
  for (Bit = B.begin(); Bit != Bend; ++Bit)
    *Bit = Vektor<double,Dim>(IpplRandom(),IpplRandom());
  Vektor<double,Dim> Bmin = min(B);
  Vektor<double,Dim> Bmax = max(B);

  testmsg << "Vektor Field B = " << B << endl << endl;
  testmsg << "min(B) = " << Bmin << endl;
  testmsg << "max(B) = " << Bmax << endl;

  return 0;
}
/***************************************************************************
 * $RCSfile: reduce.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: reduce.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
