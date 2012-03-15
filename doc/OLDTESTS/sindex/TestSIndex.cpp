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
  A benchmark program which creates a Field of given size, assigns it
  random values, and compares the time to do an operation on those Field
  elements with values above some threshold value using the standard IPPL
  'where' statement and using the SIndex mechanism.
 ***************************************************************************/

int main(int argc, char *argv[]) {
  TAU_PROFILE("main()", "int (int, char **)", TAU_DEFAULT);

  Ippl ippl(argc,argv);
  Inform testmsg(argv[0], INFORM_ALL_NODES);

  Timer t; t.start(); t.stop(); t.clear();
  float tcpu, tclk;
  float time1, time2;

  // Various counters, constants, etc:
  const unsigned D=3U;
  double hightail;
  int d;
  int tag = Ippl::Comm->next_tag(tag);
  unsigned ngrid[D];   // grid sizes sizes
  int iters = 100;

  // Layout information:
  unsigned vnodes;             // number of vnodes; input at cmd line
  e_dim_tag serialParallel[D]; // Specifies SERIAL, PARALLEL dims
  for (d=0; d<D; d++) serialParallel[d] = PARALLEL;

  // read in vnodes etc. from command line
  if (argc != (D+4)) {
    testmsg << "Usage: " << argv[0] << " <nx> <ny> <nz> <vnodes> <percent> <iterations>";
    testmsg << "\n, where <percent> = percentage of field not used in where";
    testmsg << endl;
  } else {
    for (d=0; d < D; ++d)
      ngrid[d] = atoi(argv[d+1]);
    vnodes = atoi(argv[d+1]);
    hightail = 0.01 * atof(argv[d+2]);
    iters = atoi(argv[d+3]);
    double diters = 1.0 / (double)iters;
    Index I(ngrid[0]), J(ngrid[1]), K(ngrid[2]);
    Index I2(ngrid[0]/2), J2(ngrid[1]/2), K2(ngrid[2]/2);
    NDIndex<3U> ndi(I,J,K);
    NDIndex<3U> ndi2(I2,J2,K2);
    FieldLayout<D> layout(ndi, serialParallel, vnodes);
    Field<double,D> A(layout), B(layout);

    SIndex<D> SI1(layout), SI2(layout), SI3(layout);
    assign(A, IpplRandom);
    assign(B, 0.0);

    cout << "sum of A at start = " << sum(A) << endl;
    cout << "sum of B at start = " << sum(B) << endl;

    cout << "------------------------------------------------" << endl;

    // time how long the SIndex assignment takes
    t.start();
    // SI2.reserve(1.0 - hightail);
    for (d=0; d < iters; ++d) {
      SI2 = gt(A, hightail);
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    
    cout << "time/iter for SIndex setup (SI = gt(A,percentage)) = " << tclk;
    cout << " clock, " << tcpu << " cpu " << endl;

    // time how long a regular where-statement assignment takes
    t.start();
    for (d=0; d < iters; ++d) {
      B[I][J][K] = where((gt(A[I][J][K],hightail)), A[I][J][K], B[I][J][K]);
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time1 = tclk;

    cout << "time for gt-where assign = " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after where assign = " << sum(B) << endl;

    // time how long a sparse index field assignment takes
    assign(B, 0.0);
    t.start();
    for (d=0; d < iters; ++d) {
      B[SI2] = A[SI2];
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time2 = tclk;

    cout << "SIndex equivalent " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after SIndex assign = " << sum(B) << endl;

    cout << "Ratio: " << (time2 > 0.0 ? (time1/time2) : 0.0) << endl;

    cout << "------------------------------------------------" << endl;

    // time how long the subsetted SIndex assignment takes
    t.start();
    // SI2.reserve(1.0 - hightail);
    for (d=0; d < iters; ++d) {
      SI2[ndi2] = gt(A[ndi2], hightail);
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();

    cout << "time for subsetted SIndex setup (SI[ndi2] = gt(A[ndi2],%)) = " << tclk;
    cout << " clock, " << tcpu << " cpu " << endl;

    // time how long a regular subsetted where-statement assignment takes
    assign(B, 0.0);
    t.start();
    for (d=0; d < iters; ++d) {
      B[ndi2] = where((gt(A[ndi2],hightail)), A[ndi2], B[ndi2]);
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time1 = tclk;

    cout << "time for subsetted gt-where assign = " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after where assign = " << sum(B) << endl;

    // time how long a subsetted sparse index field assignment takes
    assign(B, 0.0);
    t.start();
    for (d=0; d < iters; ++d) {
      B[SI2] = A[SI2];
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time2 = tclk;

    cout << "Subsetted SIndex equivalent " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after subsetted SIndex assign = " << sum(B) << endl;

    cout << "Ratio: " << (time2 > 0.0 ? (time1/time2) : 0.0) << endl;

    cout << "------------------------------------------------" << endl;

    // time how long it takes to do an expression assigment for the
    // different schemes
    SI2 = gt(A, hightail);

    // first for regular where
    B = 1.0;
    B[SI2] = A[SI2];
    double val1 = 1.01;
    double val2 = 1.0;
    cout << "sum of B before where expression = " << sum(B) << endl;
    t.start();
    for (d=0; d < iters; ++d) {
      B[I][J][K] *= where((gt(A[I][J][K],hightail)), val1, val2);
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time1 = tclk;

    cout << "time for gt-where expression = " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after where expression = " << sum(B) << endl;

    // time how long a sparse index field assignment takes
    B = 1.0;
    B[SI2] = A[SI2];
    cout << "sum of B before SIndex expression = " << sum(B) << endl;
    t.start();
    for (d=0; d < iters; ++d) {
      B[SI2] *= 1.01;
    }
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time2 = tclk;

    cout << "SIndex equivalent " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after SIndex expression = " << sum(B) << endl;

    cout << "Ratio: " << (time2 > 0.0 ? (time1/time2) : 0.0) << endl;

    cout << "------------------------------------------------" << endl;

    // time how long for a sindex -> particleattrib
    B = 1.0;
    B[SI2] = A[SI2];
    ParticleAttrib<double> PA;
    cout << "sum of B before SIndex->PA expression = " << sum(B) << endl;
    t.start();
    PA[SI2] = 1.01 * B[SI2];
    for (d=1; d < iters; ++d) {
      PA *= 1.01;
    }
    B[SI2] = PA[SI2];
    t.stop(); tcpu=diters*t.cpu_time(); tclk=diters*t.clock_time(); t.clear();
    time2 = tclk;
    
    cout << "SIndex->ParticleAttrib equivalent = " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    cout << "sum of B after SIndex->PA expression = " << sum(B) << endl;

    cout << "Ratio: " << (time2 > 0.0 ? (time1/time2) : 0.0) << endl;

    cout << "------------------------------------------------" << endl;

    // make sure we can do a simple PA[S] = PA[S] assignment
    cout << "Testing PAtemp[S] = PA[S] ..." << endl;
    ParticleAttrib<double> PAtemp;
    PAtemp[SI2] = PA[SI2];
    cout << "Sum(PA) = " << sum(PA) << endl;
    cout << "Sum(PAtemp) = " << sum(PAtemp) << endl;
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: TestSIndex.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:48 $
 * IPPL_VERSION_ID: $Id: TestSIndex.cpp,v 1.1.1.1 2000/11/30 21:06:48 adelmann Exp $ 
 ***************************************************************************/
