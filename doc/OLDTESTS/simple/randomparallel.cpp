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

/***************************************************************************
 *
 * The IPPL Framework
 * 
 * Test parallel random generation in the easiest way
 * 
 * Usage: randomparallel np
 ***************************************************************************/

#include "Ippl.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);

  Inform msg2all(argv[0], INFORM_ALL_NODES);
  msg2all << "Parallel random generation: Nodes= " << Ippl::getNodes() << " myNode= " << Ippl::myNode() << endl;
  Ippl::Comm->barrier();

  unsigned int p = 0;
  for (int n=0; n<atoi(argv[1]); n++) {
    if (p == Ippl::myNode()) {
      // create and init particle
      // msg2all << "processor " << Ippl::myNode() << " on particle " << n << " p=" << p << endl;
    }
    p++;
    if (p == Ippl::getNodes())
      p=0;
  }
  return 0;
}

/***************************************************************************
 * $RCSfile: random.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2000/11/30 21:06:19 $
 * IPPL_VERSION_ID: $Id: random.cpp,v 1.1.1.1 2000/11/30 21:06:19 adelmann Exp $ 
 ***************************************************************************/
