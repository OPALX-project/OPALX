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

// TestCartesianCentering.cpp
// Various tests of CartesianCentering classes.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/CartesianCentering.h"
#include "AppTypes/Vektor.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <fstream>
using namespace std;
#else
#include <fstream.h>
#endif

// forward declarations
void hardCodedOutput(char* filename); // Prototype of function defined below.
bool thediff(char* filename1, char* filename2);
extern const CenteringEnum zz[2] = {CELL, VERTEX};


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  // For writing file output to compare against hardcoded correct file output:
  Inform fdi(NULL,"text.test.TestCartesianCentering",Inform::OVERWRITE,0);

  const unsigned nx=4, ny=4, nz=4;
  Index I(nx), J(ny), K(nz);

  const unsigned ND3 = 3;
  typedef UniformCartesian<ND3> M3;
  M3 m3(I,J,K);

  const unsigned ND2 = 2;
  typedef UniformCartesian<ND2> M2;
  M2 m2(I,J);

#ifdef __MWERKS__
  // Work around CW4 default-template-parameter bug
  typedef CommonCartesianCenterings<ND2,1U,0U>::allCell CA;
#else
  typedef CommonCartesianCenterings<ND2,1U>::allCell CA;
#endif // __MWERKS__
  CenteredFieldLayout<ND2,M2,CA> clA(m2);
  Field<double, ND2, M2, CA > A(clA);
  A.print_Centerings(fdi.getStream());

  typedef CartesianCentering<zz,ND2,1U> CB;
  CenteredFieldLayout<ND2,M2,CB> clB(m2);
  Field<double, ND2, M2, CB > B(clB);
  B.print_Centerings(fdi.getStream());

  typedef CommonCartesianCenterings<ND3,1U,1U>::allFace CC;
  CenteredFieldLayout<ND3,M3,CC> clC(m3);
  Field<double, ND3, M3, CC> C(clC);
  C.print_Centerings(fdi.getStream());

#ifdef __MWERKS__
  // Work around CW4 default-template-parameter bug
  typedef CommonCartesianCenterings<ND3,3U,0U>::allVertex CD;
#else
  typedef CommonCartesianCenterings<ND3,3U>::allVertex CD;
#endif // __MWERKS__
  CenteredFieldLayout<ND3,M3,CD> clD(m3);
  Field<Vektor<double,ND3>, ND3, M3, CD> D(clD);
  D.print_Centerings(fdi.getStream());

#ifdef __MWERKS__
  // Work around CW4 default-template-parameter bug
  typedef CommonCartesianCenterings<ND3,3U,0U>::vectorFace CE;
#else
  typedef CommonCartesianCenterings<ND3,3U>::vectorFace CE;
#endif // __MWERKS__
  CenteredFieldLayout<ND3,M3,CE> clE(m3);
  Field<Vektor<double,ND3>, ND3, M3, CE> E(clE);
  E.print_Centerings(fdi.getStream());

  fdi << endl ; // Needed to flush output to file

  // Write out "by hand" into another file what the previous field-printing
  // functions should have produced; this will be compared with what they
  // actually did produce:
  hardCodedOutput("text.correct.TestCartesianCentering");

  // Compare the two files by mocking up the Unix "diff" command:
  passed = thediff("text.test.TestCartesianCentering",
		   "text.correct.TestCartesianCentering");

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
#ifdef __MWERKS__
  of 
#else
  of << "CartesianCentering: no specialized name (yet) for this case" << endl
#endif
     << "Dim = 2 ; NComponents = 1" << endl
     << "centering[dim=0][component=0] = CELL  " << endl
     << "centering[dim=1][component=0] = CELL  " << endl
#ifndef __MWERKS__
     << "CartesianCentering: no specialized name (yet) for this case" << endl
#endif // __MWERKS__
     << "Dim = 2 ; NComponents = 1" << endl
     << "centering[dim=0][component=0] = CELL  " << endl
     << "centering[dim=1][component=0] = VERTEX" << endl
#ifndef __MWERKS__
     << "CartesianCentering: no specialized name (yet) for this case" << endl
#endif // __MWERKS__
     << "Dim = 3 ; NComponents = 1" << endl
     << "centering[dim=0][component=0] = CELL  " << endl
     << "centering[dim=1][component=0] = VERTEX" << endl
     << "centering[dim=2][component=0] = CELL  " << endl
#ifndef __MWERKS__
     << "CartesianCentering: no specialized name (yet) for this case" << endl
#endif // __MWERKS__
     << "Dim = 3 ; NComponents = 3" << endl
     << "centering[dim=0][component=0] = VERTEX" << endl
     << "centering[dim=0][component=1] = VERTEX" << endl
     << "centering[dim=0][component=2] = VERTEX" << endl
     << "centering[dim=1][component=0] = VERTEX" << endl
     << "centering[dim=1][component=1] = VERTEX" << endl
     << "centering[dim=1][component=2] = VERTEX" << endl
     << "centering[dim=2][component=0] = VERTEX" << endl
     << "centering[dim=2][component=1] = VERTEX" << endl
     << "centering[dim=2][component=2] = VERTEX" << endl
#ifndef __MWERKS__
     << "CartesianCentering: no specialized name (yet) for this case" << endl
#endif // __MWERKS__
     << "Dim = 3 ; NComponents = 3" << endl
     << "centering[dim=0][component=0] = VERTEX" << endl
     << "centering[dim=0][component=1] = CELL  " << endl
     << "centering[dim=0][component=2] = CELL  " << endl
     << "centering[dim=1][component=0] = CELL  " << endl
     << "centering[dim=1][component=1] = VERTEX" << endl
     << "centering[dim=1][component=2] = CELL  " << endl
     << "centering[dim=2][component=0] = CELL  " << endl
     << "centering[dim=2][component=1] = CELL  " << endl
     << "centering[dim=2][component=2] = VERTEX" << endl
     << endl;
  of.close();
  return;
}

/***************************************************************************
 * $RCSfile: TestCartesianCentering.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:04:59 $
 * IPPL_VERSION_ID: $Id: TestCartesianCentering.cpp,v 1.1.1.1 2000/11/30 21:04:59 adelmann Exp $ 
 ***************************************************************************/
