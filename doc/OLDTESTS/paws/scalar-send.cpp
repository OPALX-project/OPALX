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
// Test code to send (and receive back) some scalars via PAWS

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform tester(argv[0]);

  tester << "Paws scalar send/receive test A" << endl;
  tester << "-------------------------------" << endl;

  // Some scalars to send and receive

  int s1 = 1, origs1 = 1;
  double s2 = 2.5, origs2 = 2.5;

  // Create a Paws connection

  tester << "Creating Paws connection object ..." << endl;
  DataConnect *paws = new PawsDataConnect("scalar-send");
  tester << "Finished creating Paws connection object." << endl;

  // Connect up the scalars

  tester << "Connecting s1 = " << s1 << " for output ..." << endl;
  ScalarDataSource<int> ps1(s1);
  paws->connect("s1", ps1, DataSource::OUTPUT);
  tester << "Connecting s2 = " << s2 << " for input ..." << endl;
  ScalarDataSource<double> ps2(s2);
  paws->connect("s2", ps2, DataSource::INPUT);

  // Wait for everything to be ready to proceed

  tester << "Waiting for ready signal ..." << endl;
  paws->ready();
  tester << "Ready complete, moving on." << endl;

  // Modify s2, and update

  s2 *= 2;
  tester << "Updating current s1 = " << s1 << " and s2 = " << s2;
  tester << " ..." << endl;
  paws->updateConnections();

  // Report the results

  tester << "Received update.  New values:" << endl;
  tester << "  s1 = " << s1 << " (should be " << origs1 << ")\n";
  tester << "  s2 = " << s2 << " (should be " << origs2 << ")\n";
  tester << endl;

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

