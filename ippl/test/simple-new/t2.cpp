/***************************************************************************
 *
 * 
 * 
 * 
 ***************************************************************************/

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc, argv);
  Inform msg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);

  int localV = Ippl::myNode();
  int globalV = 0;
  
  reduce(localV, globalV, OpAddAssign());
  
  msg2all << "nyNode= " << localV << endl;	
  msg << "Sum of all nodes = " << globalV << endl;
  return 0;

}
/***************************************************************************
 * $RCSfile: t2.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t2.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $
 ***************************************************************************/
