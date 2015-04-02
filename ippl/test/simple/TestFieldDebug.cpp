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

// TestFieldDebug.cpp , Tim Williams 10/21/1996
// This tests the functions [e,s]fp[1,2,3]() and the function setFormat() from
// Utility/FieldDebug.[h,cpp] . These are meant to be called from the debugger,
// but this function tests whether they work (for a couple of possible calls)
// within a program. It also includes specialized function definitions like the
// user of FieldDebug must have in his own source code in order to be able to
// access callable functions from the debugger, as an example for users.  This
// function also tests the setInform() function, to specify the Inform object
// used internally by FieldDebug functions. Constructing an Inform object that
// writes into a file makes it easy to do the comparson with correct output.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Field/BCond.h"
#include "Field/GuardCellSizes.h"
#include "AppTypes/Vektor.h"
#include "Utility/FieldDebug.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <fstream>
using namespace std;
#else
#include <fstream>
#endif

//-----------------------------------------------------------------------------
// User-inserted prototypes to get debugger access (examples for user ref.):
// Scalar (double):
void dfp3(BareField<double,3U> f) {fp3(f);}
void defp3(BareField<double,3U> f, int i, int j, int k) {efp3(f,i,j,k);}
void dsfp3(BareField<double,3U>& f,
	   int base1, int bound1, int stride1,
	   int base2, int bound2, int stride2,
	   int base3, int bound3, int stride3) {
  sfp3(f,base1,bound1,stride1,base2,bound2,stride2,base3,bound3,stride3);}
// Vektor (double):
void vdfp3(BareField<Vektor<double,3U>,3U> f) {fp3(f);}
void vdefp3(BareField<Vektor<double,3U>,3U> f, int i, int j, int k) {
  efp3(f,i,j,k);}
void vdsfp3(BareField<Vektor<double,3U>,3U>& f,
	    int base1, int bound1, int stride1,
	    int base2, int bound2, int stride2,
	    int base3, int bound3, int stride3) {
  sfp3(f,base1,bound1,stride1,base2,bound2,stride2,base3,bound3,stride3);}
//-----------------------------------------------------------------------------

// forward declarations
void hardCodedOutput(char* filename); // Prototype of function defined below.
bool thediff(char* filename1, char* filename2);


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  bool docomm = true; // Should try the test both ways, really....

  const unsigned Dim3=3;

  int nCells[Dim3];       // Number of grid cells in each direction
  unsigned nVNodes[Dim3]; // Number of vnodes (subdomains) in each direction.

  // Hardwired values for automated test, as in regression testing:
  for (int d=0; d<Dim3; d++) nCells[d] = 4;
  for (int d=0; d<Dim3; d++) nVNodes[d] = 2;

#ifdef TESTFIELDDEBUG_MANUAL_TEST // Set this flag for non-automated test

  // Read in the command-line arguments, appropriate number for Dim:
  int argcForDim = 1 + 2*Dim3 + 1;
  if (argc != argcForDim) {
    testmsg << "Usage: " << argv[0] 
	    << " <nx> [<ny> [<nz>]]"
	    << " <vnodesX> [<vnodesY> [<vnodesZ>]]"
	    << " <docomm (1=true,0=false)>" << endl;
    exit(1);
  }
  int av = 0;
  for (int d=0; d<Dim3; d++) nCells[d] = atoi(argv[++av]);
  for (int d=0; d<Dim3; d++) nVNodes[d] = atoi(argv[++av]);
  int docommInt = atoi(argv[++av]);
  if (docommInt == 1) {
    docomm = true;
  } else {
    docomm = false;
  }

#endif  // TESTFIELDDEBUG_MANUAL_TEST

  int nx, ny, nz;
  nx = nCells[0]; ny = nCells[1]; nz = nCells[2];
  
  Index I(nx); Index J(ny); Index K(nz);
  // Specify multipple vnodes (8) to make sure this works right:
  //  FieldLayout<Dim3> layout3(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);
  NDIndex<Dim3> ndi; ndi[0] = I; ndi[1] = J; ndi[2] = K;
  e_dim_tag serialParallelSpec[Dim3];
  for (int d=0; d<Dim3; d++) serialParallelSpec[d] = PARALLEL;
  FieldLayout<Dim3> layout3(ndi, serialParallelSpec, nVNodes);

  // New Inform-based version (tjw):
  Inform* fdip = 
    new Inform(NULL,"text.test.TestFieldDebug",Inform::OVERWRITE,0);
  Inform& fdi = *fdip;
  setInform(fdi);

  // Put guard cells and red-flag (value = -999) boundary conditions on
  // Fields, to make sure nothing funny is happening:
  GuardCellSizes<Dim3> gc(2);
  BConds<double,Dim3> sbc;
  for (int face=0; face < 2*Dim3; face++) {
    sbc[face] = new ConstantFace<double,Dim3>(face,-999.0);
  }
  BConds<Vektor<double,Dim3>,Dim3> vbc;
  for (int face=0; face < 2*Dim3; face++) {
    vbc[face] = new ConstantFace<Vektor<double,Dim3>,Dim3>(face,-999.0);
  }

  // Scalar Field -------------------------------------------------------------
  Field<double,Dim3> A3(layout3,sbc,gc);
  assign(A3[I][J][K], I + J + K);

  fdi << endl << "--------fp3(A3)-------" << endl;
  fp3(A3,docomm);

  fdi << endl << "--------sfp3(A3,nx-1,1,0,ny-1,1,0,nz-1,1)-------" << endl;
  sfp3(A3,0,nx-1,1,0,ny-1,1,0,nz-1,1,docomm);
  
  fdi << endl << "--------sfp3(A3,nx-1,1,0,ny-1,2,0,nz-1,2)-------" << endl;
  sfp3(A3,0,nx-1,1,0,ny-1,1,0,nz-1,1,docomm);


  // Vector Field--------------------------------------------------------------
  Field<Vektor<double,Dim3>,Dim3> B3(layout3,vbc,gc);
  Vektor<double, Dim3 > Vinit3(1.0,2.0,3.0);
  assign(B3,Vinit3);

  fdi << endl << "--------setFormat(1,8)-------" << endl;
  setFormat(1,8);

  fdi << endl << "--------fp3(B3)-------" << endl;
  fp3(B3,docomm);

  // Write out "by hand" into another file what the previous field-printing
  // functions should have produced; this will be compared with what they
  // actually did produce:
  hardCodedOutput("text.correct.TestFieldDebug");

  // Compare the two files by mocking up the Unix "diff" command:
  delete fdip;
  passed = thediff("text.test.TestFieldDebug","text.correct.TestFieldDebug");

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
    } else {
      same = false;                 // If file 2 ends before file 1, different
    }
  }
  return(same);
}

//-----------------------------------------------------------------------------
void hardCodedOutput(char* filename)
{
  ofstream of(filename);
  of << endl 
     << "--------fp3(A3)-------" << endl
     << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << endl
     << "==================================================K = 0" << endl
     << "--------------------------------------------------J = 0" << endl
     << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "==================================================K = 1" << endl
     << "--------------------------------------------------J = 0" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "==================================================K = 2" << endl
     << "--------------------------------------------------J = 0" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "==================================================K = 3" << endl
     << "--------------------------------------------------J = 0" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << endl
     << "" << endl
     << "" << endl
     << "--------sfp3(A3,nx-1,1,0,ny-1,1,0,nz-1,1)-------" << endl
     << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << endl
     << "==================================================K = 0" << endl
     << "--------------------------------------------------J = 0" << endl
     << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "==================================================K = 1" << endl
     << "--------------------------------------------------J = 0" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "==================================================K = 2" << endl
     << "--------------------------------------------------J = 0" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "==================================================K = 3" << endl
     << "--------------------------------------------------J = 0" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << endl
     << "" << endl
     << "" << endl
     << "--------sfp3(A3,nx-1,1,0,ny-1,2,0,nz-1,2)-------" << endl
     << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << endl
     << "==================================================K = 0" << endl
     << "--------------------------------------------------J = 0" << endl
     << "0.000e+00 1.000e+00 2.000e+00 3.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "==================================================K = 1" << endl
     << "--------------------------------------------------J = 0" << endl
     << "1.000e+00 2.000e+00 3.000e+00 4.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "==================================================K = 2" << endl
     << "--------------------------------------------------J = 0" << endl
     << "2.000e+00 3.000e+00 4.000e+00 5.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "==================================================K = 3" << endl
     << "--------------------------------------------------J = 0" << endl
     << "3.000e+00 4.000e+00 5.000e+00 6.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "4.000e+00 5.000e+00 6.000e+00 7.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "5.000e+00 6.000e+00 7.000e+00 8.000e+00 " << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "6.000e+00 7.000e+00 8.000e+00 9.000e+00 " << endl
     << "" << endl
     << "" << endl
     << "--------setFormat(1,8)-------" << endl
     << "" << endl
     << "--------fp3(B3)-------" << endl
     << "~~~~~~~~ field slice (0:3:1, 0:3:1, 0:3:1) ~~~~~~~~" << endl
     << "==================================================K = 0" << endl
     << "--------------------------------------------------J = 0" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "==================================================K = 1" << endl
     << "--------------------------------------------------J = 0" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "==================================================K = 2" << endl
     << "--------------------------------------------------J = 0" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "==================================================K = 3" << endl
     << "--------------------------------------------------J = 0" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 1" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 2" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl
     << "--------------------------------------------------J = 3" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "( 1.00000000e+00 , 2.00000000e+00 , 3.00000000e+00 )" << endl
     << "" << endl
     << "" << endl;
  of.close();
  return;
}

/***************************************************************************
 * $RCSfile: TestFieldDebug.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:38 $
 * IPPL_VERSION_ID: $Id: TestFieldDebug.cpp,v 1.1.1.1 2003/01/23 07:40:38 adelmann Exp $ 
 ***************************************************************************/
