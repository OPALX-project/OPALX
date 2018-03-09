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
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"

// set dimensionality and problem size
const unsigned Dim = 2;
int size = 20;


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Index I(size), J(size);
  FieldLayout<Dim> layout(I,J);
  Field<double,Dim> 
    F1(layout), F2(layout), F3(layout), F4(layout), F5(layout); 
  Field< Vektor<double,Dim>, Dim, UniformCartesian<Dim> > 
    V1(layout), V2(layout), V3(layout);
  Field< Tenzor<double,Dim>, Dim, UniformCartesian<Dim> >  
    T1(layout), T2(layout), T3(layout);
  Field< SymTenzor<double,Dim>, Dim, UniformCartesian<Dim> >  
    S1(layout), S2(layout), S3(layout);

  Vektor<double,Dim> Vinit1(1.0,2.0);
  Vektor<double,Dim> Vinit2(10.0,20.0);
  Vektor<double,Dim> Vinit3(100.0,200.0);

  Tenzor<double,Dim> Tinit1(1.0,2.0,2.0,3.0);
  Tenzor<double,Dim> Tinit2(10.0,20.0,20.0,30.0);
  Tenzor<double,Dim> Tinit3(100.0,200.0,200.0,300.0);

  SymTenzor<double,Dim> Sinit1(1.0,2.0,3.0);
  SymTenzor<double,Dim> Sinit2(10.0,20.0,30.0);
  SymTenzor<double,Dim> Sinit3(100.0,200.0,300.0);

  F1 = 1.0;
  F2 = 2.0;
  F3 = 3.0;
  V1 = Vinit1 ;
  V2 = Vinit2 ;
  V3 = Vinit3 ;
  T1 = Tinit1 ;
  T2 = Tinit2 ;
  T3 = Tinit3 ;
  S1 = Sinit1 ;
  S2 = Sinit2 ;
  S3 = Sinit3 ;

  // test Vektor Field operations

#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = V2 + V1;
#else
  V3 << V2 + V1;
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = V2 - V1;
#else
  V3 << V2 - V1;
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = F1 * V1;
#else
  V3 << F1 * V1;
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = V1 * F1;
#else
  V3 << V1 * F1;
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  F3 = dot(V1,V2);
#else
  F3 << dot(V1,V2);
#endif
  testmsg << " F3 " << endl << F3 << endl;

  // test Tenzor Field operations

#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = T2 + T1;
#else
  T3 << T2 + T1;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = T2 - T1;
#else
  T3 << T2 - T1;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = dot(T1,T2);
#else
  T3 << dot(T1,T2);
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  F3 = dotdot(T1,T2);
#else
  F3 << dotdot(T1,T2);
#endif
  testmsg << " F3 " << endl << F3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = dot(V1,T2);
#else
  V3 << dot(V1,T2);
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = dot(T2,V1);
#else
  V3 << dot(T2,V1);
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = F1 * T2;
#else
  T3 << F1 * T2;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = T2 * F1;
#else
  T3 << T2 * F1;
#endif
  testmsg << " T3 " << endl << T3 << endl;

  // test SymTenzor Field operations

#ifdef IPPL_USE_MEMBER_TEMPLATES
  S3 = S1 + S2;
#else
  S3 << S1 + S2;
#endif
  testmsg << " S3 " << endl << S3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  S3 = S1 - S2;
#else
  S3 << S1 - S2;
#endif
  testmsg << " S3 " << endl << S3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  S3 = dot(S1,S2);
#else
  S3 << dot(S1,S2);
#endif
  testmsg << " S3 " << endl << S3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  F3 = dotdot(S1,S2);
#else
  F3 << dotdot(S1,S2);
#endif
  testmsg << " F3 " << endl << F3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = dot(V1,S2);
#else
  V3 << dot(V1,S2);
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  V3 = dot(S2,V1);
#else
  V3 << dot(S2,V1);
#endif
  testmsg << " V3 " << endl << V3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = T1 + S2;
#else
  T3 << T1 + S2;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = T1 - S2;
#else
  T3 << T1 - S2;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = dot(T1,S2);
#else
  T3 << dot(T1,S2);
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  F3 = dotdot(T1,S2);
#else
  F3 << dotdot(T1,S2);
#endif
  testmsg << " F3 " << endl << F3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = S1 + T2;
#else
  T3 << S1 + T2;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = S2 - T1;
#else
  T3 << S2 - T1;
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  T3 = dot(S1,T2);
#else
  T3 << dot(S1,T2);
#endif
  testmsg << " T3 " << endl << T3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  F3 = dotdot(S1,T2);
#else
  F3 << dotdot(S1,T2);
#endif
  testmsg << " F3 " << endl << F3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  S3 = F1 * S2;
#else
  S3 << F1 * S2;
#endif
  testmsg << " S3 " << endl << S3 << endl;
#ifdef IPPL_USE_MEMBER_TEMPLATES
  S3 = S2 * F1;
#else
  S3 << S2 * F1;
#endif
  testmsg << " S3 " << endl << S3 << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: FVecOps.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: FVecOps.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
