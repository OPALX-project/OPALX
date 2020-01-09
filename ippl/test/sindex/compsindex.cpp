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

#include "Ippl.h"

// Benchmarks where() against SIndex for a particular example, where sparsness
// and compression should be high. For testing Bill's new code making SIndex
// compression-smarter.

//-----------------------------------------------------------------------------
// User-inserted prototypes to get debugger access (examples for user ref.):
// Scalar (double):
void dfp3(BareField<double,3U> &f) {fp3(f);}
void defp3(BareField<double,3U> &f, int i, int j, int k) {efp3(f,i,j,k);}
void dsfp3(BareField<double,3U> &f,
	   int base1, int bound1, int stride1,
	   int base2, int bound2, int stride2,
	   int base3, int bound3, int stride3) {
  sfp3(f,base1,bound1,stride1,base2,bound2,stride2,base3,bound3,stride3);}
//-----------------------------------------------------------------------------


int main(int argc, char *argv[]) {

  Ippl ippl(argc,argv);
  Inform pout(NULL,0);

  Timer t; t.start(); t.stop(); t.clear();
  float tcpu, tclk;
  float tcreatecpu, tcreateclk;
  float tsetupcpu, tsetupclk;
  float tassigncpu, tassignclk;
  float time1, time2;

  // Various counters, constants, etc:
  const unsigned D=3U;
  double hightail;
  int d;
  int tag = Ippl::Comm->next_tag(tag);
  unsigned ngrid[D];   // grid sizes sizes

  // Layout information:
  unsigned vnodes;             // number of vnodes; input at cmd line
  e_dim_tag serialParallel[D]; // Specifies SERIAL, PARALLEL dims
  for (d=0; d<D; d++) serialParallel[d] = PARALLEL;

  // read in vnodes etc. from command line
  if (argc != (D+3)) {
    pout << "Usage: " << argv[0] << " <nx> <ny> <nz> <vnodes> <percent>";
    pout << "\n, where <percent> = percentage of field not used in where";
    pout << endl;
  } else {
    for (d=0; d < D; ++d)
      ngrid[d] = atoi(argv[d+1]);
    vnodes = atoi(argv[d+1]);
    hightail = 0.01 * atof(argv[d+2]);
    Index I(ngrid[0]), J(ngrid[1]), K(ngrid[2]);
    NDIndex<3U> ndi(I,J,K);
    FieldLayout<D> layout(ndi, serialParallel, vnodes);
    Field<double,D> A(layout), B(layout);

    // time how long the SIndex creation takes
    t.start();
    SIndex<D> SI(layout);
    t.stop(); tcreatecpu=t.cpu_time(); tcreateclk=t.clock_time(); t.clear();

    // Assign A to be random in the lower eighth of the 3D box, zero elsewhere,
    // which gives good compression with enough vnodes:
    A = 0.0;
    NDIndex<D> lowerEighth;
    for (d=0; d<D; d++) lowerEighth[d] = Index(ngrid[d]/8);
    A[lowerEighth] = IpplRandom;
    pout << endl << "A.CompressedFraction() = " 
	 << A.CompressedFraction() << endl;

    // Assign B to zero, so it starts fully compressed; should next uncompress
    // only lower eighth of it:
    B = 0.0;

    //    pout << "sum of A at start = " << sum(A) << endl;
    //    pout << "sum of B at start = " << sum(B) << endl;

    // time how long the SIndex assignment takes
    t.start();
    SI = gt(A, hightail);
    t.stop(); tsetupcpu=t.cpu_time(); tsetupclk=t.clock_time(); t.clear();

    // time how long a regular where-statement assignment takes
    t.start();
    B[I][J][K] = where((gt(A[I][J][K],hightail)), A[I][J][K], B[I][J][K]);
    t.stop(); tcpu= t.cpu_time(); tclk=t.clock_time(); t.clear();
    time1 = tclk;

    pout << endl << "time for where(gt(A[I][J][K]," << hightail 
	 << ")) assign = " << tclk << " clock, "
         << tcpu << " cpu " << endl;
    pout << "sum of B after where assign = " << sum(B) << endl;

    // Reassign B to zero, so it starts fully compressed; should next
    // uncompress only lower eighth of it:
    B = 0.0;

    // time how long a sparse index field assignment takes
    t.start();
    B[SI] = A[SI];
    t.stop(); tassigncpu=t.cpu_time(); tassignclk=t.clock_time(); t.clear();
    time2 = tcreateclk + tsetupclk + tassignclk;

    pout << endl << "time for SIndex creation = " << tcreateclk;
    pout << " clock, " << tcreatecpu << " cpu " << endl;
    pout << "time for SIndex setup (SI = gt(A," << hightail << ")) = " 
	 << tsetupclk;
    pout << " clock, " << tsetupcpu << " cpu " << endl;
    pout << "time for [SIndex] assign " << tassignclk << " clock, "
         << tassigncpu << " cpu " << endl;
    pout << "sum of B after SIndex assign = " << sum(B) << endl;

    pout << endl << "========================================" << endl << endl;

    pout << "where() timing:                     " <<time1 << " clock" << endl;

    pout << "SIndex timing: create+setup+assign: " <<time2 << " clock" << endl;

    pout << "where()/SIndex Ratio:               " 
	 << (time2 > 0.0 ? (time1/time2) : 0.0) << endl << endl;
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

