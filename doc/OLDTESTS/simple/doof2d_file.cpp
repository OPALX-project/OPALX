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

  int sizeX      = 100, sizeY   = 100;
  int centerX    =  5, centerY = 5;
  int iterations =  10;

  const unsigned Dim=2;
  Index I(sizeX);
  Index J(sizeY);
  Index PrintI(4);
  Index PrintJ(4);
  FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL,2*Ippl::getNodes());
  Field<double,Dim>  A(layout,GuardCellSizes<Dim>(1));
  Field<double,Dim> A2(layout,GuardCellSizes<Dim>(1));
  Field<double,Dim>  B(layout);

  // set up a DiscField for A
  DiscField<Dim> df("Ack", "Ack.out.config", 2, "Field<double,2>");

  // make other connections, using default connect method
  DataConnect *Bcon = DataConnectCreator::create("B Connection");
  B.connect("B", Bcon);

  A = 0.0;
  A[centerX][centerY] = 1.0*iterations;
  double fact = 1.0/9.0;

  for(int iter = 0 ; iter < iterations ; iter++ ) {
    msg << "Iteration " << iter << ": Computing new A ..." << endl;
    assign(B[I][J],
	   fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] + 
		 A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] + 
		 A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1])
	   );
    assign(A,B);
    msg << "  iter = " << iter << ", sum = " << sum(A) << endl;
    msg << "  new A = " << A[PrintI][PrintJ] << endl;
    msg << "  Writing A to file ..." << endl;
    df.write(A,0);
    B = -2.0;
    df.write(B,1);
    msg << "  Updating connection ..." << endl;
    B.updateConnection();
    if ((iter+1) % 10 == 0)
      A.interact();
  }

  // now read the data back into a Field
  msg << "....................................................." << endl;
  DiscField<Dim> readf("Ack", "Ack.in.config");
  unsigned int nrecords = readf.get_NumRecords();
  unsigned int nfields  = readf.get_NumFields();
  NDIndex<Dim> nsize    = readf.get_Domain();
  msg << "Reading back data, records = " << nrecords << ", fields = ";
  msg << nfields << ", total domain = " << nsize << endl;

  if (nfields > 0 && nrecords > 0 && nsize.size() > 0) {
    // read all records and print them out
    FieldLayout<Dim> rdlayout(nsize);
    Field<double,Dim> rdfield(rdlayout);
    for (unsigned int r=0; r < nrecords; ++r) {
      readf.read(rdfield, 0, r);
      msg << "Read record " << r << ": " << rdfield[PrintI][PrintJ] << endl;
    }
  } else {
    // there was an error reading the file
    msg << "An error was encountered while trying to read the file." <<endl;
  }

  delete Bcon;
  return 0;
}

/***************************************************************************
 * $RCSfile: doof2d_file.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:05:49 $
 * IPPL_VERSION_ID: $Id: doof2d_file.cpp,v 1.1.1.1 2000/11/30 21:05:49 adelmann Exp $ 
 ***************************************************************************/
