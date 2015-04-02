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
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  const unsigned Dim=2;
  Index I(5),J(5);
  FieldLayout<Dim> layout(I,J);
  typedef UniformCartesian<Dim,double> Mesh;
#ifdef __MWERKS__
  // Is this a serious bug in Metrowerks CodeWarrior 4? --tjw
  Field<double,Dim,Mesh,Mesh::DefaultCentering> S1(layout),S2(layout),
    S3(layout);
  Field<Vektor<double,Dim>,Dim,Mesh,Mesh::DefaultCentering> V1(layout);
#else
  Field<double,Dim,Mesh> S1(layout),S2(layout),S3(layout);
  Field<Vektor<double,Dim>,Dim,Mesh> V1(layout);
#endif // __MWERKS__

  S1[I][J] = I+10*J ;
  S2[I][J] = -I-10*J ;
#ifdef __MWERKS__
  // Note: operator=() doesn't work for either MWERKS or SGI here
  assign(V1[I][J](0), S1[I][J]);
  assign(V1[I][J](1), S2[I][J]*10.0) ;
#else
  V1[I][J](0) << S1[I][J] ;
  V1[I][J](1) << S2[I][J]*10.0 ;
#endif // __MWERKS__
  S1[I][J] = V1[I][J](0) - (I+10*J);
  S2[I][J] = V1[I][J](1)/10.0 + (I+10*J);
  S1 *= S1;
  S2 *= S2;
  double s1 = sum(S1);
  double s2 = sum(S2);
  testmsg << ( (s1<1e-10) ? "PASSED" : "FAILED" ) << endl;
  testmsg << ( (s2<1e-10) ? "PASSED" : "FAILED" ) << endl;

#ifdef __MWERKS__
  // Note: operator=() doesn't work for either MWERKS or SGI here
  assign(V1[I][J](0), I);
  assign(V1[I][J](1), 27.5);
#else
  V1[I][J](0)  << I ;
  V1[I][J](1)  << 27.5 ;
#endif // __MWERKS__
  S1[I][J] = V1[I][J](0) - I;
  S2[I][J] = V1[I][J](1) - 27.5;
  S1 *= S1;
  S2 *= S2;
  s1 = sum(S1);
  s2 = sum(S2);
  testmsg << ( (s1<1e-10) ? "PASSED" : "FAILED" ) << endl;
  testmsg << ( (s2<1e-10) ? "PASSED" : "FAILED" ) << endl;

  S1[I][J] = I+10*J ;
  S2[I][J] = -I-10*J ;
  V1 = Vektor<double,2>(0,0);
  V1[I][J](0) += S1[I][J] ;
  V1[I][J](1) += S2[I][J]*10.0 ;
  S1[I][J] = V1[I][J](0) - (I+10*J);
  S2[I][J] = V1[I][J](1)/10.0 + (I+10*J);
  S1 *= S1;
  S2 *= S2;
  s1 = sum(S1);
  s2 = sum(S2);
  testmsg << ( (s1<1e-10) ? "PASSED" : "FAILED" ) << endl;
  testmsg << ( (s2<1e-10) ? "PASSED" : "FAILED" ) << endl;
 
  V1 = Vektor<double,2>(0,0);
  V1[I][J](0) += I ;
  V1[I][J](1) += 27.5 ;
  S1[I][J] = V1[I][J](0) - I;
  S2[I][J] = V1[I][J](1) + (-27.5);
  S1 *= S1;
  S2 *= S2;
  s1 = sum(S1);
  s2 = sum(S2);
  testmsg << ( (s1<1e-10) ? "PASSED" : "FAILED" ) << endl;
  testmsg << ( (s2<1e-10) ? "PASSED" : "FAILED" ) << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: component.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: component.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
