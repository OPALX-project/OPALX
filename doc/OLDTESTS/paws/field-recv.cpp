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
// Test code to receive some fields via PAWS

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform tester(argv[0]);

  tester << "Paws field send/receive test B" << endl;
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
  Field<float, 3> ca1(layout1);
  Field<int, 3> ca2(layout1);
  Field<float, 3> ca3(layout2);

  // Initialize the Fields

  a1 = 0;
  a2 = 0;
  a3 = 0;
  ca1[I][J][K] = 100 * (K + 1) + 10 * (J + 1) + I + 1;
  ca2 = ca1 + 1000;
  ca3[I2][J2][K2] = ca1[I2][J2][K2];

  // Create a Paws connection

  tester << "Creating Paws connection object ..." << endl;
  DataConnect *paws = new PawsDataConnect("field-send");
  tester << "Finished creating Paws connection object." << endl;

  // Connect up the scalars

  tester << "Connecting s1 = " << s1 << " for input ..." << endl;
  ScalarDataSource<int> ps1(s1);
  paws->connect("s1", ps1, DataSource::INPUT);
  tester << "Connecting s2 = " << s2 << " for output ..." << endl;
  ScalarDataSource<double> ps2(s2);
  paws->connect("s2", ps2, DataSource::OUTPUT);
  tester << "Connecting iters = " << iters << " for input ..." << endl;
  ScalarDataSource<int> piters(iters);
  paws->connect("iters", piters, DataSource::INPUT);

  // Connect up the Fields

  tester << "Connecting a1 for input ..." << endl;
  paws->connect("a1", a1, DataSource::INPUT);
  tester << "Connecting a2 for input ..." << endl;
  paws->connect("a2", a2, DataSource::INPUT);
  tester << "Connecting a1view for input ..." << endl;
  paws->connect("a1view", a3, DataSource::INPUT);

  // Wait for everything to be ready to proceed

  tester << "Waiting for ready signal ..." << endl;
  paws->ready();
  tester << "Ready complete, moving on." << endl;

  // Modify s1, and update

  s1 *= 2;
  tester << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester << ", plus arrays ..." << endl;
  paws->updateConnections();

  // Report the results

  tester << "Received update.  New values:" << endl;
  tester << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester << "  iters = " << iters << "\n";
  tester << endl;

  // Do, in a loop, updates of the receiver. Add one to the arrays each time
  // and compare to the results.

  int myiters = iters;
  while (myiters-- > 0)
    {
      ca1 = ca1 + 1;
      ca2 = ca2 + 1;
      ca3 = ca3 + 1;
      tester << "Receiving for iters = " << myiters << endl;
      paws->updateConnections();
      float sum1 = sum(a1);
      int sum2 = sum(a2);
      float sum3 = sum(a3);
      float csum1 = sum(ca1);
      int csum2 = sum(ca2);
      float csum3 = sum(ca3);
      tester << "a1 OK: " << sum1 << " == " << csum1 << endl;
      tester << "a2 OK: " << sum2 << " == " << csum2 << endl;
      tester << "a3 OK: " << sum3 << " == " << csum3 << endl;
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

