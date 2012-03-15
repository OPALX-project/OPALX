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

/***************************************************************************
  A simple program to test SubField assignments using sparse index operations
  (with SIndex, SOffset, and NDIndex objects as the subsetting mechanism).
 ***************************************************************************/

int main(int argc, char *argv[]) {
  TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  const unsigned Dim=2;
  Index I(4);
  Index I2(2);
  Index J(4);
  Index K(2);
  FieldLayout<Dim> layout(I, J, PARALLEL, PARALLEL, 4*Ippl::getNodes());
  FieldLayout<3> layout3(I, J, K, PARALLEL, PARALLEL, PARALLEL,
			  2*Ippl::getNodes());
  BConds<double,Dim> dbc;
  for (int f=0; f < 2*Dim; f++) dbc[f] = new ZeroFace<double,Dim>(f);
  Field<double,Dim> A(layout, dbc, GuardCellSizes<Dim>(1));
  A = 0.0;
  BConds<bool,Dim> bbc;
  for (int f=0; f < 2*Dim; f++) bbc[f] = new ZeroFace<bool,Dim>(f);
  Field<bool,Dim>   B(layout, bbc, GuardCellSizes<Dim>(1));
  B = true;
  Field<double,Dim> C(layout, dbc, GuardCellSizes<Dim>(1));
  C = 0.0;

  BConds<float,3> fbc;
  for (int f=0; f < 2*3; f++) fbc[f] = new ZeroFace<float,3>(f);
  Field<float,3>    A3(layout3, fbc, GuardCellSizes<3>(2));
  A3 = 0.0;
  Field<float,3>    B3(layout3, fbc, GuardCellSizes<3>(2));
  B3 = 0.0;
  SIndex<Dim> s1(layout);
  s1.addIndex(NDIndex<Dim>(Index(1,1), Index(1,2)));

  SIndex<Dim> s2 = s1(1,-1);
  s2.addIndex(SOffset<Dim>(0,0));

  testmsg << "\n************ testing 2D SubParticleAttrib<> ************" << endl;
  A[I][J] = (I+1) + (J*10);
  testmsg << "s1 = " << s1 << endl;
  testmsg << " A = " << A << endl;
  ParticleAttrib<double> PA;
  PA[s1] = A[s1];
  testmsg << "Result of PA[s1] = A[s1] : PA = " << endl;
  testmsg << PA[s1] << endl;

  testmsg << "s2 = " << s2 << endl;
  PA[s2] = A[s2] + A[s2];
  testmsg << "Result of PA[s2] = A[s2] + A[s2] : PA = " << endl;
  testmsg << PA[s2] << endl;

  PA[s2] *= (A[s2] + PA[s2]);
  testmsg << "Result of PA[s2] *= A[s2] + PA[s2] : PA = " << endl;
  testmsg << PA[s2] << endl;

  A = 1.0;
  SIndex<Dim> sX = s1;
  sX.addIndex(SOffset<Dim>(0,0));
  testmsg << "Created sX = " << sX << endl;

  A[sX] += PA[sX];
  testmsg << "Result of A[sX] = 1 + PA[sX] : A = " << A << endl;

  A = 0.0;
  C = 1.0;
  A[sX] = C[sX];
  testmsg << "Result of A[sX] = C[sX], after C = 1.0: A = " << A << endl;

  ParticleAttrib<Vektor<double,2> > PA2;
  PA2(0)[sX] = A[sX] + PA[sX];
  PA2(1)[sX] = -(A[sX] + PA[sX]);
  testmsg << "Result of PA2[sX] = +- (A[sX] + PA[sX]):" << endl;
  testmsg << "  A = " << A << endl;
  testmsg << " PA = " << PA[sX] << endl;
  testmsg << "PA2 = " << PA2[sX] << endl;

  testmsg << "\n************ testing 2D SubBareField<NDIndex> ************" << endl;

  Index IX(0,2);
  Index JX(2,2);
  NDIndex<Dim> sub1(IX, JX);
  SubBareField<double,Dim,NDIndex<Dim> > SA1(A, sub1);
  testmsg << "Created SubBareField SA1 from A and NDIndex " << sub1 << ":" <<endl;
  testmsg << " SA1 = " << SA1 << endl;

  SA1 = 3.0;
  testmsg << "SA1 set to 3:" << endl;
  testmsg << "   A = " << A << endl;

  NDIndex<Dim> sub2(IX + 1, JX);
  SubBareField<double,Dim,NDIndex<Dim> > SA2(A, sub2);
  testmsg << "Created SubBareField SA2 from A and NDIndex " << sub2 << ":" <<endl;
  testmsg << " SA2 = " << SA2 << endl;

  SA1 = SA2 + 4.0;
  testmsg << "SA1 on " << sub1 << " set to (SA2 + 4):" << endl;
  testmsg << "   A = " << A << endl;

  testmsg << "\n************ testing 2D SubBareField<SIndex> ************" << endl;

  SubBareField<double,Dim,SIndex<Dim> > SB1(A, s1);
  testmsg << "Created SubBareField SB1 from A and SIndex " << s1 << ":" << endl;
  testmsg << " SB1 = " << SB1 << endl;

  SB1 = -2.0;
  testmsg << "SB1 set to -2:" << endl;
  testmsg << "   A = " << A << endl;

  SubBareField<double,Dim,SIndex<Dim> > SB2(A, s1(1,-1));
  testmsg << "Created SubBareField SB2 from A and SIndex " << SB2.getDomain();
  testmsg << ":" << endl;
  testmsg << " SB2 = " << SB2 << endl;

  SB1 = SB2 + 4.0;
  testmsg << "SB1 set to (SB2 + 4):" << endl;
  testmsg << "   A = " << A << endl;

  testmsg << "\n************ testing 2D SubBareField<SOffset> ************" << endl;

  SubBareField<double,Dim,SOffset<Dim> > SC1(A, SOffset<Dim>(2,1));
  testmsg << "Created SubBareField SC1 from A and SOffset "  << SC1.getDomain();
  testmsg << ":" << endl;
  testmsg << " SC1 = " << SC1 << endl;

  SC1 = 10.0;
  testmsg << "SC1 set to 10:" << endl;
  testmsg << "   A = " << A << endl;

  SubBareField<double,Dim,SOffset<Dim> > SC2(A, SOffset<Dim>(2,2));
  SubBareField<double,Dim,SOffset<Dim> > SC3(A, SOffset<Dim>(3,2));
  SC1 = SC1 + (SC2 + SC3) * 10;
  testmsg << "SC1 set to 10 + [(2,2) + (3,2)]*10:" << endl;
  testmsg << "   A = " << A << endl;

  testmsg << "\n************ testing 2D Field[SIndex] ************" << endl;

  testmsg << "Current SIndex s1 = " << s1 << endl;
  testmsg << "Results of A[s1] = -1000:" << endl;
  A[s1] = -1000.0;
  testmsg << "   A = " << A << endl;
  B[s1] = ne(A[s1], A[s1]);
  testmsg << "Results of B[s1] = (A[s1] != A[s1]):" << endl;
  testmsg << "   B = " << B << endl;
  C[s1] = 1.0;
  testmsg << "Results of C[s1] = 1:" << endl;
  testmsg << "   C = " << C << endl;
  C[s1] = A[s1(-1,-1)] - 2000;
  testmsg << "Results of C[s1] = (A[s1(-1,-1)] - 2000) :" << endl;
  testmsg << "   C = " << C << endl;

  testmsg << "\n************ testing 2D Field[SIndex,compressed] ****" << endl;

  Field<double,Dim> A2(layout);
  A2 = 1.0;
  testmsg << "   A2 = " << A2 << endl;
  testmsg << "   A2.compressedFraction = " << A2.CompressedFraction() << endl;
  A2[I2][I2] = -10.0;
  testmsg << "Initial settings:" << endl;
  testmsg << "   A2 = " << A2 << endl;
  testmsg << "   A2.compressedFraction = " << A2.CompressedFraction() << endl;
  s1 = ne(A2, 0.0);
  testmsg << "   s1 = " << s1 << endl;
  testmsg << "   A2 = " << A2 << endl;
  testmsg << "   A2.compressedFraction = " << A2.CompressedFraction() << endl;
  A2[s1] = A2[s1] + 100.0;
  testmsg << "Results of A2[s1] = A2[s1] + 100.0 :" << endl;
  testmsg << "   A2 = " << A2 << endl;
  testmsg << "   A2.compressedFraction = " << A2.CompressedFraction() << endl;
  testmsg << "   s1 = " << s1 << endl;

  testmsg << "\n************ testing 3D Field[SIndex] ************" << endl;

  A3[I][J][K] = I + J + K;
  B3 = 0.0;

  SIndex<3> sindx3(layout3);
  sindx3 = eq(A3, 2.0);

  SOffset<3> offset;
  offset[0] = 1;
  offset[1] = 0;
  offset[2] = (-1);

  testmsg << "Current 3D SIndex sindx3 = " << sindx3 << endl;
  testmsg << "Current 3D A3 = " << A3 << endl;
  testmsg << "A3 after A3[sindx3] = -1:" << endl;
  A3[sindx3] = -1.0;
  testmsg << "   A3 = " << A3 << endl;
  B3[sindx3] = 1.0;
  testmsg << "Results of B3[sindx3] = 1:" << endl;
  testmsg << "   B3 = " << B3 << endl;
  B3[sindx3] = A3[sindx3(offset)] - 2000.0;
  testmsg << "Results of B3[sindx3] = A3[sindx3(1,0,-1)] - 2000.0:" << endl;
  testmsg << "   B3 = " << B3 << endl;

  return 0;
}

/***************************************************************************
 * $RCSfile: subfield.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
 * IPPL_VERSION_ID: $Id: subfield.cpp,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
 ***************************************************************************/
