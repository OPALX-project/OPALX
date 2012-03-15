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

// TestFieldNDIndexing.cpp , Tim Williams 11/06/1996
// Tests the use of an NDIndex object in single brackets instead of N Index
// objects in their own brackets for indexing a Field. Includes testing of
// index arithmetic on the NDIndex objects.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/GuardCellSizes.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool success = true;
  const unsigned Dim=2;
  Index I(5),J(5);
  FieldLayout<Dim> layout(I,J);
  GuardCellSizes<Dim> gc(2);
  typedef UniformCartesian<Dim> Mesh;
#ifdef __MWERKS__
  // Is this a serious bug in Metrowerks CodeWarrior 4? --tjw
  Field<double,Dim,Mesh,typename Mesh::DefaultCentering> S1(layout,gc);
  Field< Vektor<double,Dim>, Dim, Mesh, typename Mesh::DefaultCentering >
    V1(layout,gc);
#else
  Field<double,Dim,Mesh> S1(layout,gc);
  Field<Vektor<double,Dim>,Dim,Mesh> V1(layout,gc);
#endif // __MWERKS__

  NDIndex<Dim> allCells, interiorCells;
  allCells[0] = I;
  allCells[1] = J;
  Index interiorI(1,3,1), interiorJ(1,3,1);
  interiorCells[0] = interiorI;
  interiorCells[1] = interiorJ;

  assign(S1,0.0);
  Vektor<double, Dim> vector0(0.0, 0.0), vector1(1.0, 1.0), vector9(9.0, 9.0);
  assign(V1,vector0);

  assign(S1[interiorCells],1.0);
  assign(V1[interiorCells],vector1);

  double summit;
  summit = sum(S1);
  //  testmsg << ( (summit==9.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==9.0;
  // Kludge because sum(Field<Vektor<...>...>) doesn't work:
  Vektor<double, Dim> v1sum;
  assign(S1[allCells], V1[allCells](0));
  v1sum(0) = sum(S1);
  assign(S1[allCells], V1[allCells](1));
  v1sum(1) = sum(S1);
  //  testmsg << ( (v1sum==vector9) ? "PASSED" : "FAILED" ) << endl;
  success = success && v1sum==vector9;

  //------------------------------------------
  // NDIndex arithmetic:

  // [NDIndex + int*] and [NDIndex - int*]:
  Field<double,Dim,Mesh> S2(layout,gc);
  assign(S1,1.0);
  assign(S2,1.0);
  int ones[Dim] = {1, 1};
  assign(S1[interiorCells], 
	 S1[interiorCells + ones] + S1[interiorCells - ones]);
  assign(S1, S1 - S2);
  summit = sum(S1);
  //  testmsg << ( (summit==9.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==9.0;
  // [int* + NDIndex] (and [NDIndex - int*]):
  assign(S1,1.0);
  assign(S1[interiorCells], 
	 S1[ones + interiorCells] + S1[interiorCells - ones]);
  assign(S1, S1 - S2);
  summit = sum(S1);
  //  testmsg << ( (summit==9.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==9.0;
  // [NDIndex * int*]:
  NDIndex<Dim> XTwoCells;
  Index XTwoI(2,2,1);
  XTwoCells[0] = XTwoI;
  XTwoCells[1] = J;
  assign(S1,1.0);
  assign(S1[XTwoCells], 2.0);
  int twoone[2] = {2, 1};
  assign(S1[XTwoCells], S1[XTwoCells*twoone]);
  summit = sum(S1);
  //  testmsg << ( (summit==25.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==25.0;
  // [int* * NDIndex]:
  assign(S1,1.0);
  assign(S1[XTwoCells], 2.0);
  assign(S1[XTwoCells], S1[twoone*XTwoCells]);
  summit = sum(S1);
  //  testmsg << ( (summit==25.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==25.0;
  // [-NDIndex]:
  Index IN(-3,3,1),JN(-3,3,1);
  FieldLayout<Dim> layoutN(IN,JN);
  Field<double,Dim,Mesh> S1N(layoutN,gc);
  NDIndex<Dim> TwoByTwoCells;
  Index YTwoByTwoJ(1,2,1);
  Index XTwoByTwoI(1,2,1);
  TwoByTwoCells[0] = XTwoByTwoI;
  TwoByTwoCells[1] = YTwoByTwoJ;
  assign(S1N,1.0);
  assign(S1N[TwoByTwoCells], 2.0);
  assign(S1N[TwoByTwoCells], S1N[-TwoByTwoCells]);
  summit = sum(S1N);
  //  testmsg << ( (summit==49.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==49.0;
  // [NDIndex / int*]:
  Index I8(8),J8(8);
  FieldLayout<Dim> layout8(I8,J8);
  Field<double,Dim,Mesh> S8(layout8,gc);
  NDIndex<Dim> FourSixFourSixCells;
  Index XFourSixFourSixI(4,6,2);
  Index YFourSixFourSixJ(4,6,2);
  FourSixFourSixCells[0] = XFourSixFourSixI;
  FourSixFourSixCells[1] = YFourSixFourSixJ;
  assign(S8,1.0);
  assign(S8[FourSixFourSixCells], 2.0);
  int twotwo[2] = {2, 2};
  assign(S8[FourSixFourSixCells], S8[FourSixFourSixCells/twotwo]);
  summit = sum(S8);
  //  testmsg << ( (summit==64.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==64.0;
  // [int* - NDIndex]:
  NDIndex<Dim> OneTwoOneTwoCells;
  Index XOneTwoOneTwoI(1,2,1);
  Index YOneTwoOneTwoJ(1,2,1);
  OneTwoOneTwoCells[0] = XOneTwoOneTwoI;
  OneTwoOneTwoCells[1] = YOneTwoOneTwoJ;
  assign(S8,1.0);
  assign(S8[OneTwoOneTwoCells], 2.0);
  int eighteight[2] = {8, 8};
  assign(S8[OneTwoOneTwoCells], S8[eighteight - OneTwoOneTwoCells]);
  summit = sum(S8);
  //  testmsg << ( (summit==64.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==64.0;
  // Check the special "alias" form of the NDIndex constructor.
  // This example duplicates the one just above, but uses an
  // intermediary NDIndex object "SixSevenSixSevenCells" which is constructed
  // using the special "alias" form of the constructor.
  assign(S8,1.0);
  assign(S8[OneTwoOneTwoCells], 2.0);
  NDIndex<Dim> SixSevenSixSevenCells(eighteight - OneTwoOneTwoCells);
  assign(S8[OneTwoOneTwoCells], S8[SixSevenSixSevenCells]);
  summit = sum(S8);
  //  testmsg << ( (summit==64.0) ? "PASSED" : "FAILED" ) << endl;
  success = success && summit==64.0;
  //------------------------------------------
  // print overall test result
  testmsg << (success ? "PASSED" : "FAILED") << endl; 
  return 0;
}

/***************************************************************************
 * $RCSfile: TestFieldNDIndexing.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:11 $
 * IPPL_VERSION_ID: $Id: TestFieldNDIndexing.cpp,v 1.1.1.1 2000/11/30 21:05:11 adelmann Exp $ 
 ***************************************************************************/
