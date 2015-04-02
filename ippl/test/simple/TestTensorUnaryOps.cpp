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

// TestTensorUnaryOps.cpp , Tim Williams 10/19/1999
// Tests trace(), det(), transpose(), cofactors() with Tenzor, SymTenzor,
// AntiSymTenzor operands.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "AppTypes/AntiSymTenzor.h"

// define helper class for performing tests in 1D, 2D, and 3D
template <int Dim>
class testTensorUnaryOps
{
public:
  static bool apply(Inform& testmsg)
  {
    bool passedAllTests = true;

    NDIndex<Dim> domain;
    int nCellsTot = 1;
    for (int d = 0; d < Dim; d++) { 
      domain[d] = Index(5);
      nCellsTot *= 5;
    }
    FieldLayout<Dim> layout(domain);
    typedef UniformCartesian<Dim> Mesh;
    typedef Tenzor<double,Dim> FT_t;
    typedef AntiSymTenzor<double,Dim> AT_t;
    typedef SymTenzor<double,Dim> ST_t;
    Field<FT_t,Dim,Mesh> tff(layout);
    Field<AT_t,Dim,Mesh> tfa(layout);
    Field<ST_t,Dim,Mesh> tfs(layout);

    // Assign values:
    Tenzor<double,Dim> tf, tfTranspose;
    AntiSymTenzor<double,Dim> ta, taTranspose;
    SymTenzor<double,Dim> ts, tsTranspose;
    double fullSymTrace = 0.0;
    for (int i = 0; i < Dim; i++) {
      for (int j = 0; j < Dim; j++) {
	tf(i,j) = (i+1)*(i+1) + (j+1)*(j+1) + (i+4)*(j+4) + i;
	if (i == j) fullSymTrace += tf(i,j);
	tfTranspose(j,i) = tf(i,j);
      }
    }
    ta = tf;
    ts = tf;
    tff = tf;
    tfa = ta;
    tfs = ts;
    for (int i = 0; i < Dim; i++) {
      for (int j = 0; j < Dim; j++) {
	taTranspose(j,i) = ta(i,j);
	tsTranspose(j,i) = ts(i,j);
      }
    }

    // Test determinant of Tenzor:
    PInsist(Dim<4, "[Sym]Tenzor det() function not implemented for Dim>3!");
    double detValue = sum(det(tff));
    //  testmsg << "detValue = " << detValue << endl;
    switch (Dim) {
    case 1:
      if (detValue != 18*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tff)) = " << detValue
		<< " != 18*nCellsTot = " << 18*nCellsTot << endl;
      }
      break;
    case 2:
      if (detValue != -38*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tff)) = " << detValue
		<< " != -38*nCellsTot = " << -38*nCellsTot << endl;
      }
      break;
    case 3:
      if (detValue != -4*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tff)) = " << detValue
		<< " != -4*nCellsTot = " << -4*nCellsTot << endl;
      }
      break;
    default:
      ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
      break;
    }

    // Test determinant of AntiSymTenzor
    double detValueA = sum(det(tfa));
    //  testmsg << "detValueA = " << detValueA << endl;
    switch (Dim) {
    case 1:
      if (detValueA != 0) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfa)) = " << detValueA << " != 0."
		<< endl;
      }
      break;
    case 2:
      if (detValueA != -ta(1,0)*ta(0,1)*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfa)) = " << detValueA
		<< " != -ta(1,0)*ta(0,1)*nCellsTot = " 
		<< -ta(1,0)*ta(0,1)*nCellsTot << endl;
      }
      break;
    case 3:
      if (detValueA != 0) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfa)) = " << detValueA << " != 0."
		<< endl;
      }
      break;
    default:
      ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
      break;
    }

    // Test determinant of SymTenzor
    double detValueS = sum(det(tfs));
    switch (Dim) {
    case 1:
      if (detValueS != 18*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfs)) = " << detValueS
		<< " != 18*nCellsTot = " << 18*nCellsTot << endl;
      }
      break;
    case 2:
      if (detValueS != -38.25*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfs)) = " << detValueS
		<< " != -38.25*nCellsTot = " << -38.25*nCellsTot << endl;
      }
      break;
    case 3:
      if (detValueS != -4*nCellsTot) {
	passedAllTests = false;
	testmsg << Dim << "D, sum(det(tfs)) = " << detValueS
		<< " != -4*nCellsTot = " << -4*nCellsTot << endl;
      }
      break;
    default:
      ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
      break;
    }

    // Test trace of Tenzor:
    double traceValue;
    traceValue = sum(trace(tff));
    if (traceValue != fullSymTrace*nCellsTot) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(trace(tff) = " << traceValue
	      << " != fullSymTrace*nCellsTot = " << fullSymTrace*nCellsTot
	      << endl;
    }

    // Test trace of AntiSymTenzor:
    traceValue = sum(trace(tfa));
    if (traceValue != 0) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(trace(tfa) = " << traceValue
	      << " != 0." << endl;
    }

    // Test trace of SymTenzor:
    traceValue = sum(trace(tfs));
    if (traceValue != fullSymTrace*nCellsTot) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(trace(tfs) = " << traceValue
	      << " != fullSymTrace*nCellsTot = " << fullSymTrace*nCellsTot
	      << endl;
    }

    // Test transpose of Tenzor:
    Tenzor<double,Dim> transposeValue;
    transposeValue = sum(transpose(tff));
    if (transposeValue != tfTranspose*nCellsTot) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(transpose(tff) = " << transposeValue
	      << " != tfTranspose*nCellsTot = " << tfTranspose*nCellsTot
	      << endl;
    }

    // Test transpose of AntiSymTenzor:
    AntiSymTenzor<double,Dim> transposeValueA;
    transposeValueA = sum(transpose(tfa));
    if (transposeValueA != taTranspose*nCellsTot) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(transpose(tfa) = " << transposeValueA
	      << " != taTranspose*nCellsTot = " << taTranspose*nCellsTot
	      << endl;
    }

    // Test transpose of SymTenzor:
    SymTenzor<double,Dim> transposeValueS;
    transposeValueS = sum(transpose(tfs));
    if (transposeValueS != tsTranspose*nCellsTot) {
      passedAllTests = false;
      testmsg << Dim << "D, sum(transpose(tfs) = " << transposeValueS
	      << " != tsTranspose*nCellsTot = " << tsTranspose*nCellsTot
	      << endl;
    }

    // Test cofactors of Tenzor:
    Tenzor<double,Dim> cofactorsValue;
    cofactorsValue = sum(cofactors(tff))/nCellsTot;
    // Check results by computing det using all possible Laplace expansions,
    // and comparing to directly-computed det value:
    Tenzor<double,Dim> sumValue;
    sumValue = sum(tff);
    double altDetValue;
    // Laplace expansions using rows:
    for (int i = 0; i < Dim; i++) {
      altDetValue = 0.0;
      for (int j = 0; j < Dim; j++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValue) {
	passedAllTests = false;
	testmsg << Dim << "D, i=" << i << ", sum(cofactors(tff) = " 
		<< cofactorsValue << " ; Laplace expansion det = " 
		<< altDetValue << " != det = " << detValue << endl;
      }
    }
    // Laplace expansions using columns:
    for (int j = 0; j < Dim; j++) {
      altDetValue = 0.0;
      for (int i = 0; i < Dim; i++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValue) {
	passedAllTests = false;
	testmsg << Dim << "D, j=" << j << ", sum(cofactors(tff) = " 
		<< cofactorsValue << " ; Laplace expansion det = "
		<< altDetValue << " != det = " << detValue << endl;
      }
    }

    // Test cofactors of AntiSymTenzor:
    cofactorsValue = sum(cofactors(tfa))/nCellsTot;
    // Check results by computing det using all possible Laplace expansions,
    // and comparing to directly-computed det value:
    sumValue = sum(tfa);
    // Laplace expansions using rows:
    for (int i = 0; i < Dim; i++) {
      altDetValue = 0.0;
      for (int j = 0; j < Dim; j++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValueA) {
	passedAllTests = false;
	testmsg << Dim << "D, i=" << i << ", sum(cofactors(tfa) = " 
		<< cofactorsValue << " ; Laplace expansion det = " 
		<< altDetValue << " != det = " << detValueA << endl;
      }
    }
    // Laplace expansions using columns:
    for (int j = 0; j < Dim; j++) {
      altDetValue = 0.0;
      for (int i = 0; i < Dim; i++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValueA) {
	passedAllTests = false;
	testmsg << Dim << "D, j=" << j << ", sum(cofactors(tfa) = " 
		<< cofactorsValue << " ; Laplace expansion det = "
		<< altDetValue << " != det = " << detValueA << endl;
      }
    }

    // Test cofactors of SymTenzor:
    cofactorsValue = sum(cofactors(tfs))/nCellsTot;
    // Check results by computing det using all possible Laplace expansions,
    // and comparing to directly-computed det value:
    sumValue = sum(tfs);
    // Laplace expansions using rows:
    for (int i = 0; i < Dim; i++) {
      altDetValue = 0.0;
      for (int j = 0; j < Dim; j++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValueS) {
	passedAllTests = false;
	testmsg << Dim << "D, i=" << i << ", sum(cofactors(tfs) = " 
		<< cofactorsValue << " ; Laplace expansion det = " 
		<< altDetValue << " != det = " << detValueS << endl;
      }
    }
    // Laplace expansions using columns:
    for (int j = 0; j < Dim; j++) {
      altDetValue = 0.0;
      for (int i = 0; i < Dim; i++) {
	altDetValue += sumValue(i,j)*cofactorsValue(i,j);
      }
      if (altDetValue != detValueS) {
	passedAllTests = false;
	testmsg << Dim << "D, j=" << j << ", sum(cofactors(tfs) = " 
		<< cofactorsValue << " ; Laplace expansion det = "
		<< altDetValue << " != det = " << detValueS << endl;
      }
    }

    return passedAllTests;
  }
};


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  bool passed = true;
  if (!(testTensorUnaryOps<1>::apply(testmsg))) {
    passed = false;
    testmsg << "problem in 1D run" << endl;
  }
  if (!(testTensorUnaryOps<2>::apply(testmsg))) {
    passed = false;
    testmsg << "problem in 2D run" << endl;
  }
  if (!(testTensorUnaryOps<3>::apply(testmsg))) {
    passed = false;
    testmsg << "problem in 3D run" << endl;
  }

  testmsg << ( passed ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: TestTensorUnaryOps.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: TestTensorUnaryOps.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $ 
 ***************************************************************************/
