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
// Test code to send some fields via PAWS

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform tester(argv[0]);

  tester << "Paws field send/receive test A" << endl;
  tester << "-------------------------------" << endl;

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;
  int iters = 10;

  // Some Field's to send and receive

  Index I(2), J(4), K(8);
  Index I2(1), J2(2), K2(2);
  e_dim_tag serpar[3];
  serpar[0] = SERIAL;
  serpar[1] = PARALLEL;
  serpar[2] = SERIAL;
  NDIndex<3> ND1(I,J,K), ND2(I2, J2, K2);
  FieldLayout<3> layout1(ND1, serpar);
  FieldLayout<3> layout2(ND2, serpar);
  Field<float, 3> a1(layout1);
  Field<int, 3> a2(layout1);
  Field<float, 3> a3(layout2);

  // Initialize the Fields

  a1[I][J][K] = 100 * (K + 1) + 10 * (J + 1) + I + 1;
  a2 = a1 + 1000;
  a3[I2][J2][K2] = a1[I2][J2][K2];

  // Create a Paws connection

  tester << "Creating Paws connection object ..." << endl;
  DataConnect *paws = new PawsDataConnect("field-send");
  tester << "Finished creating Paws connection object." << endl;

  // Connect up the scalars

  tester << "Connecting s1 = " << s1 << " for output ..." << endl;
  ScalarDataSource<int> ps1(s1);
  paws->connect("s1", ps1, DataSource::OUTPUT);
  tester << "Connecting s2 = " << s2 << " for input ..." << endl;
  ScalarDataSource<double> ps2(s2);
  paws->connect("s2", ps2, DataSource::INPUT);
  tester << "Connecting iters = " << iters << " for output ..." << endl;
  ScalarDataSource<int> piters(iters);
  paws->connect("iters", piters, DataSource::OUTPUT);

  // Connect up the Fields

  tester << "Connecting a1 for output ..." << endl;
  paws->connect("a1", a1, DataSource::OUTPUT);
  tester << "Connecting a2 for output ..." << endl;
  paws->connect("a2", a2, DataSource::OUTPUT);
  tester << "Connecting a1view for output ..." << endl;
  paws->connect("a1view", a3, DataSource::OUTPUT);

  // Wait for everything to be ready to proceed

  tester << "Waiting for ready signal ..." << endl;
  paws->ready();
  tester << "Ready complete, moving on." << endl;

  // Modify s2, and update

  s2 *= 2;
  tester << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester << ", plus arrays ..." << endl;
  paws->updateConnections();

  // Report the results

  tester << "Received update.  New values:" << endl;
  tester << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester << "  iters = " << iters << "\n";
  tester << endl;

  // Do, in a loop, updates of the receiver. Add one to the arrays each time.

  int myiters = iters;
  while (myiters-- > 0)
    {
      a1 = a1 + 1.0;
      a2 = a2 + 1;
      a3 = a3 + 1.0;
      tester << "Sending for iters = " << myiters << endl;
      paws->updateConnections();
    }

  // Delete PAWS connection, disconnecting us from the other code.

  tester << "Closing paws connection ..." << endl;
  delete paws;

  // Finish up

  tester << "-------------------------------------------" << endl;
  return 0;
}


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

