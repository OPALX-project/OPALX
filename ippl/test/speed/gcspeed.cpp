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
  a new constant value, and then calculates the time needed to fill
  guard cells.
 ***************************************************************************/

int main(int argc, char *argv[]) {
  // profiling macros
  

  unsigned int i, d;

  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Timer t; t.start(); t.stop(); t.clear();
  float time1, time2;

  // Various counters, constants, etc:
  const unsigned D=3U;
  unsigned ngrid[D];   // grid sizes sizes
  unsigned gcsize=0;

  // Layout information:
  unsigned vnodes;             // number of vnodes; input at cmd line
  e_dim_tag serialParallel[D]; // Specifies SERIAL, PARALLEL dims

  // read in vnodes etc. from command line
  if (argc < (D+3)) {
    testmsg << "Usage: " << argv[0] << " <nx> <ny> <nz> <vnodes> <gcsize>";
    testmsg << endl;
  } else {
    NDIndex<D> ndi;
    for (d=0; d < D; ++d) {
      ngrid[d] = atoi(argv[d+1]);
      ndi[d] = Index(ngrid[d]);
      serialParallel[d] = PARALLEL;
    }
    vnodes = atoi(argv[d+1]);
    gcsize = atoi(argv[d+2]);
    FieldLayout<D> layout(ndi, serialParallel, vnodes);
    Field<double,D> A(layout, GuardCellSizes<D>(gcsize));
    // Field<double,D> B(layout, GuardCellSizes<D>(gcsize));
    Field<double,D> B(layout);

    testmsg << "Assigning A constant values 3 ... 102 ..." << endl;
    t.start();
    for (i=3; i < 103; ++i) {
     assign(A, IpplRandom);
    }
    t.stop(); time1 = t.clock_time(); t.clear();

    testmsg << "Assigning B constant values 3 ... 102 ..." << endl;
    t.start();
    for (i=3; i < 103; ++i) {
     assign(B, IpplRandom);
    }
    t.stop(); time2 = t.clock_time(); t.clear();

    testmsg << "Time for A loop (w/ gcsize=" << gcsize << " : ";
    testmsg << time1 << endl;
    testmsg << "Time for B loop (w/ gcsize=" << 0 << " : ";
    testmsg << time2 << endl;

    float avgratio = time1/time2;
    reduce(avgratio, avgratio, OpAddAssign());
    testmsg << "Average Ratio: ";
    testmsg << avgratio / (float)(Ippl::getNodes()) << endl;
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

