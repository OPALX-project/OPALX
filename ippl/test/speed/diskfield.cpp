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

int main(int argc, char *argv[]) {
  Ippl ippl(argc,argv);
  Inform msg(argv[0]);

  int vnodes = -1;
  int size = 100;
  int i, iterations = 10;

  if (argc > 4) {
    msg << "Usage: " << argv[0] << " [<size>] [<vnodes>] [<iters>]" << endl;
    exit(1);
  }

  if (argc > 1)
    size = atoi(argv[1]);
  if (argc > 2)
    vnodes = atoi(argv[2]);
  if (argc > 3)
    iterations = atoi(argv[3]);

  const unsigned Dim = 3;

  double ratebase = 2 * iterations * size * size * size * sizeof(int);

  Index I(size);
  Index J(size);
  Index K(size);
  FieldLayout<Dim> layout(I,J,K, PARALLEL, PARALLEL, PARALLEL, vnodes);
  Field<int,Dim>  A(layout);
//  Field<int,Dim>  B(layout, GuardCellSizes<Dim>(1));

  msg << ">>> DiscField benchmark: 3D fields, " << size << "^3 elements, on ";
  msg << vnodes << " vnodes <<<" << endl;

  // initialize A and B
  A[I][J][K] = I * J - K;
//  B = 1 - A;

  // set up a DiscField
  DataConnect *dc =new FileDataConnect("DiscField.config",1,"Field<double,3>");
  dc->connect("A", A, DataSource::OUTPUT);
//  dc->connect("B", B, DataSource::OUTPUT);

  // write A & B out to disk 100 times, and compute bytes/usec rate
  msg << "Writing " << iterations << " iterations of A and B:" << endl;
  Timer mytimer;
  mytimer.clear();
  mytimer.start();
  for (i=0; i < iterations; ++i) {
    A.updateConnection();
    A.updateConnection();
//    B.updateConnection();
  }
  mytimer.stop();
  delete dc;
  double rate = (ratebase/1000000.0) / mytimer.clock_time();
  double prate = (DiscBuffer::writebytes/1000000.0) / DiscBuffer::writetime;
  msg << "diskfield write rate = " << rate << " MB/s" << endl;
  msg << "    IPPL write rate = " << prate << " MB/s" << endl;

  // set up a DiscField for reading
  dc = new FileDataConnect("DiscField.config");
  dc->connect("A", A, DataSource::INPUT);
//  dc->connect("B", B, DataSource::INPUT);

  // read A & B from disk iteration times, and compute bytes/usec rate
  msg << "Reading " << iterations << " iterations of A and B:" << endl;
  mytimer.clear();
  mytimer.start();
  for (i=0; i < iterations; ++i) {
    A.updateConnection();
    A.updateConnection();
//    B.updateConnection();
  }
  mytimer.stop();
  delete dc;
  rate = (ratebase/1000000.0) / mytimer.clock_time();
  prate = (DiscBuffer::readbytes/1000000.0) / DiscBuffer::readtime;
  msg << "DiscField read rate = " << rate << " MB/s" << endl;
  msg << "    IPPL read rate = " << prate << " MB/s" << endl;


  // sanity check on A and B
  /*
  msg << "Checking accuracy of A ..." << endl;
  B[I][J][K] = I * J - K;
  int MSD = sum((A-B)*(A-B));
  msg << "A ok? " << (MSD == 0) << " (MSD = " << MSD << ")" << endl;
  */

  return 0;
}

/***************************************************************************
 * $RCSfile: diskfield.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:40 $
 * IPPL_VERSION_ID: $Id: diskfield.cpp,v 1.1.1.1 2003/01/23 07:40:40 adelmann Exp $ 
 ***************************************************************************/
