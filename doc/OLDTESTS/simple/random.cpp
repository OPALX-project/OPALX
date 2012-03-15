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
#include <stdlib.h>


int main(int argc, char *argv[])
{
  int i;
  Ippl ippl(argc,argv);

  Inform testmsg(argv[0], INFORM_ALL_NODES);
  testmsg << "Testing output of RNGRand ..." << endl;

  RNGRand rng;
  for (i = 0; i < 10; i++)
    testmsg << "  " << rng() << endl;

  testmsg << "Compare with output of rand() ..." << endl;
  srand(1);
  for (i = 0; i < 10; i++) {
    double num = double(rand()) / double(RAND_MAX+1);
    testmsg << "  " << num << endl;
  }

  testmsg << "Compare with output of IpplRandom ..." << endl;
  for (i = 0; i < 10; i++)
    testmsg << "  " << IpplRandom() << endl;


  testmsg << "Sum of IpplRandom values ..." << endl;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (i = 0; i < 100000; i++){ 
    sum1 += 1.0 - (2.0*IpplRandom());
    sum2 += 1.0 - (2.0*rng());
  }
  testmsg << "sum1 = " << sum1/100000 << " sum2 = " << sum2/100000 << endl;
  

  return 0;
}

/***************************************************************************
 * $RCSfile: random.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:19 $
 * IPPL_VERSION_ID: $Id: random.cpp,v 1.1.1.1 2000/11/30 21:06:19 adelmann Exp $ 
 ***************************************************************************/
