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
// ----------------------------------------------------------------------------
// The IPPL Framework - Visit http://people.web.psi.ch/adelmann/ for more details
// 
// This program was prepared by the Regents of the University of California at
// TestBC2.cpp
// Various tests of BCond* classes (see also TestBC.cpp)

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"
#include "Utility/FieldDebug.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <fstream>
using namespace std;
#else
#include <fstream.h>
#endif

// forward declarations
void hardCodedOutput(char* filename); // Prototype of function defined below.
bool thediff(char* filename1, char* filename2);


int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  // For writing file output to compare against hardcoded correct file output:
  Inform* fdip = new Inform(NULL,"text.test.TestBC2",Inform::OVERWRITE,0);
  Inform& fdi = *fdip;
  setInform(fdi);
  setFormat(9,1,0);

  const unsigned Dim = 2;
  int nCells = 4;
  NDIndex<Dim> allVerts;
  NDIndex<Dim> allCells;
  for (int d=0; d<Dim; d++) allCells[d] = Index(nCells);
  for (int d=0; d<Dim; d++) allVerts[d] = Index(nCells + 1);
  e_dim_tag serialParallel[Dim];
  for (int d=0; d<Dim; d++) serialParallel[d] = PARALLEL;
  unsigned vnodes = 4;    // Hardwire 4 vnodes
  FieldLayout<Dim> layoutCell(allCells, serialParallel, vnodes);
  FieldLayout<Dim> layoutVert(allVerts, serialParallel, vnodes);
  typedef UniformCartesian<Dim> M;

  // Set Cell-centered boundary conditions.
  BConds<double,Dim,M,Cell> cbc;
  for (int face = 0; face < 2*Dim; face++) {
    cbc[face] = new ZeroGuardsAndZeroFace<double,Dim,M,Cell>(face);
  }
  // Cell-centered test Field's:
  Field<double,Dim,M,Cell> cA(layoutCell,GuardCellSizes<Dim>(2),cbc);

  // Set Vert-centered boundary conditions.
  BConds<double,Dim,M,Vert> vbc;
  for (int face = 0; face < 2*Dim; face++) {
    vbc[face] = new ZeroGuardsAndZeroFace<double,Dim,M,Vert>(face);
  }
  // Vert-centered test Field's:
  Field<double,Dim,M,Vert> vA(layoutVert,GuardCellSizes<Dim>(2),vbc);

  // Assign reference values:
  cA = 1.0;
  vA = 2.0;

  // Print reference values, including global guard layers, to see that the
  // proper elements are zeroed out (first, print BConds objects):
  fdi << "++++++++BConds object cbc begin++++++++" << endl;
  fdi << cbc;
  fdi << "++++++++BConds object cbc end++++++++++" << endl << endl;
  fdi << "++++++++++cA+++++++++++" << endl ;
  ggfp2(cA);
  fdi << "++++++++BConds object vbc begin++++++++" << endl;
  fdi << vbc;
  fdi << "++++++++BConds object vbc end++++++++++" << endl << endl;
  fdi << "++++++++++vA+++++++++++" << endl ;
  ggfp2(vA);

  // --------------------------------------------------------------------------
  // Now test the LinearExtrapolateFace BC on the Vert-Centered Field<Vektor>:

  // Construct the Field and BConds:
  BConds<Vektor<double,Dim>,Dim,M,Vert> vvbc;
  for (int face = 0; face < 2*Dim; face++) {
    vvbc[face] = 
      new LinearExtrapolateFace<Vektor<double,Dim>,Dim,M,Vert>(face);
  }
  fdi << "++++++++BConds object vvbc begin++++++++" << endl;
  fdi << vvbc;
  fdi << "++++++++BConds object vvbc end++++++++++" << endl << endl;
  Field<Vektor<double,Dim>,Dim,M,Vert> 
    vvA(layoutVert,GuardCellSizes<Dim>(2),vvbc);

  // Load a reasonable set of spatially-dependent values, like mesh node
  // coordinates:
  for (int d=0; d<Dim; d++) {
    assign(vvA[allVerts](d), allVerts[d]);
  }

  // Print reference values, including global guard layers, to check results:
  setFormat(3,1,0);
  fdi << "++++++++++vvA+++++++++++" << endl ;
  ggfp2(vvA);

  // Now do the same test, but use componentwise BC:

  // Construct the Field and BConds:
  BConds<Vektor<double,Dim>,Dim,M,Vert> vvbcComponent;
  int counter = 0;
  for (int d=0; d<Dim; d++) {
    int loFace = 2*d;
    int hiFace = 2*d + 1;
    for (int i=0; i<Dim; i++) {
      vvbcComponent[counter++] = 
	new ComponentLinearExtrapolateFace
	<Vektor<double,Dim>,Dim,M,Vert>(loFace, i);
      vvbcComponent[counter++] = 
	new ComponentLinearExtrapolateFace
	<Vektor<double,Dim>,Dim,M,Vert>(hiFace, i);
    }
  }
  fdi << "++++++++BConds object vvbcComponent begin++++++++" << endl;
  fdi << vvbcComponent;
  fdi << "++++++++BConds object vvbcComponent end++++++++++" << endl << endl;
  Field<Vektor<double,Dim>,Dim,M,Vert> 
    vvB(layoutVert,GuardCellSizes<Dim>(2),vvbcComponent);

  // Load a reasonable set of spatially-dependent values, like mesh node
  // coordinates:
  for (int d=0; d<Dim; d++) {
    assign(vvB[allVerts](d), allVerts[d]);
  }

  // Print reference values, including global guard layers, to check results:
  setFormat(3,1,0);
  fdi << "++++++++++vvB+++++++++++" << endl ;
  ggfp2(vvB);


  // Test LinearExtrapolateFace BC on the cell-centered Field<double>.

  // Construct the Field and BConds:
  BConds<double,Dim,M,Cell> cbc2;
  for (int face = 0; face < 2*Dim; face++) {
    cbc2[face] = new LinearExtrapolateFace<double,Dim,M,Cell>(face);
  }
  fdi << "++++++++BConds object cbc2 begin++++++++" << endl;
  fdi << cbc2;
  fdi << "++++++++BConds object cbc2 end++++++++++" << endl << endl;
  Field<double,Dim,M,Cell> cB(layoutCell,GuardCellSizes<Dim>(2),cbc2);

  // Assign reference values:
  cB = 3.0;

  // Print reference values, including global guard layers, to check results:
  setFormat(9,1,0);
  fdi << "++++++++++cB+++++++++++" << endl ;
  ggfp2(cB);

  // Uncompress one of the vnodes and look again::
  NDIndex<Dim> oneVnode;
  for (int d=0; d<Dim; d++) { oneVnode[d] = Index(0,1,1); }
  for (int d=0; d<Dim; d++) {
    cB[oneVnode] = cB[oneVnode] + oneVnode[d];
  }

  // Print reference values, including global guard layers, to check results:
  fdi << "++++++++++cBafter+++++++++++" << endl;
  ggfp2(cB);

  // --------------------------------------------------------------------------

  // Write out "by hand" into another file what the previous field-printing
  // functions should have produced; this will be compared with what they
  // actually did produce:
  hardCodedOutput("text.correct.TestBC2");

  // Compare the two files by mocking up the Unix "diff" command:
  delete fdip;
  passed = thediff("text.test.TestBC2",
		   "text.correct.TestBC2");

  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

//-----------------------------------------------------------------------------
// Mock up the Unix "diff" utility to compare two files:
//-----------------------------------------------------------------------------
bool thediff(char* filename1, char* filename2)
{
  bool same = true;
  char ch1, ch2;
  ifstream file1(filename1);
  ifstream file2(filename2);
  while (file1.get(ch1)) {          // Read file 1 char-by-char until eof
    if (file2.get(ch2)) {           // Read equivalent char from file 2
      if (ch1 != ch2) same = false; // If they're different,files are different
    }
    else {
      same = false;                 // If file 2 ends before file 1, different
    }
  }
  return same;
}

//-----------------------------------------------------------------------------
void hardCodedOutput(char* filename)
{
  ofstream of(filename);
  of << "++++++++BConds object cbc begin++++++++" << endl;
  of << "BConds:(" << endl;
  of << "ZeroGuardsAndZeroFace, Face=0 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=1 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=2 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=3" << endl;
  of << ")" << endl;
  of << "" << endl;
  of << "++++++++BConds object cbc end++++++++++" << endl;
  of << "" << endl;
  of << "++++++++++cA+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:5:1, -2:5:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "0.0e+00 0.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "0.0e+00 0.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "0.0e+00 0.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "0.0e+00 0.0e+00 1.0e+00 1.0e+00 1.0e+00 1.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 "
     << endl;
  of << "" << endl;
  of << "++++++++BConds object vbc begin++++++++" << endl;
  of << "BConds:(" << endl;
  of << "ZeroGuardsAndZeroFace, Face=0 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=1 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=2 , " << endl;
  of << "ZeroGuardsAndZeroFace, Face=3" << endl;
  of << ")" << endl;
  of << "" << endl;
  of << "++++++++BConds object vbc end++++++++++" << endl;
  of << "" << endl;
  of << "++++++++++vA+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:6:1, -2:6:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 2.0e+00 2.0e+00 2.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 2.0e+00 2.0e+00 2.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 2.0e+00 2.0e+00 2.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 6" << endl;
  of << "0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00 0.0e+00"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "++++++++BConds object vvbc begin++++++++" << endl;
  of << "BConds:(" << endl;
  of << "LinearExtrapolateFace, Face=0 , " << endl;
  of << "LinearExtrapolateFace, Face=1 , " << endl;
  of << "LinearExtrapolateFace, Face=2 , " << endl;
  of << "LinearExtrapolateFace, Face=3" << endl;
  of << ")" << endl;
  of << "" << endl;
  of << "++++++++BConds object vvbc end++++++++++" << endl;
  of << "" << endl;
  of << "++++++++++vvA+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:6:1, -2:6:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "( -2.0e+00 , -2.0e+00 ) ( -1.0e+00 , -2.0e+00 ) ( 0.0e+00 , -2.0e+00 )"
     << endl;
  of << "( 1.0e+00 , -2.0e+00 ) ( 2.0e+00 , -2.0e+00 ) ( 3.0e+00 , -2.0e+00 )"
     << endl;
  of << "( 4.0e+00 , -2.0e+00 ) ( 5.0e+00 , -2.0e+00 ) ( 6.0e+00 , -2.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "( -2.0e+00 , -1.0e+00 ) ( -1.0e+00 , -1.0e+00 ) ( 0.0e+00 , -1.0e+00 )"
     << endl;
  of << "( 1.0e+00 , -1.0e+00 ) ( 2.0e+00 , -1.0e+00 ) ( 3.0e+00 , -1.0e+00 )"
     << endl;
  of << "( 4.0e+00 , -1.0e+00 ) ( 5.0e+00 , -1.0e+00 ) ( 6.0e+00 , -1.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "( -2.0e+00 , 0.0e+00 ) ( -1.0e+00 , 0.0e+00 ) ( 0.0e+00 , 0.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 0.0e+00 ) ( 2.0e+00 , 0.0e+00 ) ( 3.0e+00 , 0.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 0.0e+00 ) ( 5.0e+00 , 0.0e+00 ) ( 6.0e+00 , 0.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "( -2.0e+00 , 1.0e+00 ) ( -1.0e+00 , 1.0e+00 ) ( 0.0e+00 , 1.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 1.0e+00 ) ( 2.0e+00 , 1.0e+00 ) ( 3.0e+00 , 1.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 1.0e+00 ) ( 5.0e+00 , 1.0e+00 ) ( 6.0e+00 , 1.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "( -2.0e+00 , 2.0e+00 ) ( -1.0e+00 , 2.0e+00 ) ( 0.0e+00 , 2.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 2.0e+00 ) ( 2.0e+00 , 2.0e+00 ) ( 3.0e+00 , 2.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 2.0e+00 ) ( 5.0e+00 , 2.0e+00 ) ( 6.0e+00 , 2.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "( -2.0e+00 , 3.0e+00 ) ( -1.0e+00 , 3.0e+00 ) ( 0.0e+00 , 3.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 3.0e+00 ) ( 2.0e+00 , 3.0e+00 ) ( 3.0e+00 , 3.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 3.0e+00 ) ( 5.0e+00 , 3.0e+00 ) ( 6.0e+00 , 3.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "( -2.0e+00 , 4.0e+00 ) ( -1.0e+00 , 4.0e+00 ) ( 0.0e+00 , 4.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 4.0e+00 ) ( 2.0e+00 , 4.0e+00 ) ( 3.0e+00 , 4.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 4.0e+00 ) ( 5.0e+00 , 4.0e+00 ) ( 6.0e+00 , 4.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "( -2.0e+00 , 5.0e+00 ) ( -1.0e+00 , 5.0e+00 ) ( 0.0e+00 , 5.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 5.0e+00 ) ( 2.0e+00 , 5.0e+00 ) ( 3.0e+00 , 5.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 5.0e+00 ) ( 5.0e+00 , 5.0e+00 ) ( 6.0e+00 , 5.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 6" << endl;
  of << "( -2.0e+00 , 6.0e+00 ) ( -1.0e+00 , 6.0e+00 ) ( 0.0e+00 , 6.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 6.0e+00 ) ( 2.0e+00 , 6.0e+00 ) ( 3.0e+00 , 6.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 6.0e+00 ) ( 5.0e+00 , 6.0e+00 ) ( 6.0e+00 , 6.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "++++++++BConds object vvbcComponent begin++++++++" << endl;
  of << "BConds:(" << endl;
  of << "ComponentLinearExtrapolateFace, Face=0 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=1 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=0 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=1 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=2 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=3 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=2 , " << endl;
  of << "ComponentLinearExtrapolateFace, Face=3" << endl;
  of << ")" << endl;
  of << "" << endl;
  of << "++++++++BConds object vvbcComponent end++++++++++" << endl;
  of << "" << endl;
  of << "++++++++++vvB+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:6:1, -2:6:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "( -2.0e+00 , -2.0e+00 ) ( -1.0e+00 , -2.0e+00 ) ( 0.0e+00 , -2.0e+00 )"
     << endl;
  of << "( 1.0e+00 , -2.0e+00 ) ( 2.0e+00 , -2.0e+00 ) ( 3.0e+00 , -2.0e+00 )"
     << endl;
  of << "( 4.0e+00 , -2.0e+00 ) ( 5.0e+00 , -2.0e+00 ) ( 6.0e+00 , -2.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "( -2.0e+00 , -1.0e+00 ) ( -1.0e+00 , -1.0e+00 ) ( 0.0e+00 , -1.0e+00 )"
     << endl;
  of << "( 1.0e+00 , -1.0e+00 ) ( 2.0e+00 , -1.0e+00 ) ( 3.0e+00 , -1.0e+00 )"
     << endl;
  of << "( 4.0e+00 , -1.0e+00 ) ( 5.0e+00 , -1.0e+00 ) ( 6.0e+00 , -1.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "( -2.0e+00 , 0.0e+00 ) ( -1.0e+00 , 0.0e+00 ) ( 0.0e+00 , 0.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 0.0e+00 ) ( 2.0e+00 , 0.0e+00 ) ( 3.0e+00 , 0.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 0.0e+00 ) ( 5.0e+00 , 0.0e+00 ) ( 6.0e+00 , 0.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "( -2.0e+00 , 1.0e+00 ) ( -1.0e+00 , 1.0e+00 ) ( 0.0e+00 , 1.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 1.0e+00 ) ( 2.0e+00 , 1.0e+00 ) ( 3.0e+00 , 1.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 1.0e+00 ) ( 5.0e+00 , 1.0e+00 ) ( 6.0e+00 , 1.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "( -2.0e+00 , 2.0e+00 ) ( -1.0e+00 , 2.0e+00 ) ( 0.0e+00 , 2.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 2.0e+00 ) ( 2.0e+00 , 2.0e+00 ) ( 3.0e+00 , 2.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 2.0e+00 ) ( 5.0e+00 , 2.0e+00 ) ( 6.0e+00 , 2.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "( -2.0e+00 , 3.0e+00 ) ( -1.0e+00 , 3.0e+00 ) ( 0.0e+00 , 3.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 3.0e+00 ) ( 2.0e+00 , 3.0e+00 ) ( 3.0e+00 , 3.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 3.0e+00 ) ( 5.0e+00 , 3.0e+00 ) ( 6.0e+00 , 3.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "( -2.0e+00 , 4.0e+00 ) ( -1.0e+00 , 4.0e+00 ) ( 0.0e+00 , 4.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 4.0e+00 ) ( 2.0e+00 , 4.0e+00 ) ( 3.0e+00 , 4.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 4.0e+00 ) ( 5.0e+00 , 4.0e+00 ) ( 6.0e+00 , 4.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "( -2.0e+00 , 5.0e+00 ) ( -1.0e+00 , 5.0e+00 ) ( 0.0e+00 , 5.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 5.0e+00 ) ( 2.0e+00 , 5.0e+00 ) ( 3.0e+00 , 5.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 5.0e+00 ) ( 5.0e+00 , 5.0e+00 ) ( 6.0e+00 , 5.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 6" << endl;
  of << "( -2.0e+00 , 6.0e+00 ) ( -1.0e+00 , 6.0e+00 ) ( 0.0e+00 , 6.0e+00 )"
     << endl;
  of << "( 1.0e+00 , 6.0e+00 ) ( 2.0e+00 , 6.0e+00 ) ( 3.0e+00 , 6.0e+00 )"
     << endl;
  of << "( 4.0e+00 , 6.0e+00 ) ( 5.0e+00 , 6.0e+00 ) ( 6.0e+00 , 6.0e+00 )"
     << endl;
  of << "" << endl;
  of << "" << endl;
  of << "++++++++BConds object cbc2 begin++++++++" << endl;
  of << "BConds:(" << endl;
  of << "LinearExtrapolateFace, Face=0 , " << endl;
  of << "LinearExtrapolateFace, Face=1 , " << endl;
  of << "LinearExtrapolateFace, Face=2 , " << endl;
  of << "LinearExtrapolateFace, Face=3" << endl;
  of << ")" << endl;
  of << "" << endl;
  of << "++++++++BConds object cbc2 end++++++++++" << endl;
  of << "" << endl;
  of << "++++++++++cB+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:5:1, -2:5:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "++++++++++cBafter+++++++++++" << endl;
  of << "~~~~~~~~ field slice (-2:5:1, -2:5:1) ~~~~~~~~" << endl;
  of << "--------------------------------------------------J = -2" << endl;
  of << "-1.0e+00 0.0e+00 1.0e+00 2.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = -1" << endl;
  of << "0.0e+00 1.0e+00 2.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 0" << endl;
  of << "1.0e+00 2.0e+00 3.0e+00 4.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 1" << endl;
  of << "2.0e+00 3.0e+00 4.0e+00 5.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 2" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 3" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 4" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "--------------------------------------------------J = 5" << endl;
  of << "3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 3.0e+00 "
     << endl;
  of << "" << endl;
  of << "" << endl;
  
  of.close();
  return;
}

/***************************************************************************
 * $RCSfile: TestBC2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:04:50 $
 * IPPL_VERSION_ID: $Id: TestBC2.cpp,v 1.1.1.1 2000/11/30 21:04:50 adelmann Exp $ 
 ***************************************************************************/
