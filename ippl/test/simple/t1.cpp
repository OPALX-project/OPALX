/***************************************************************************
 *
 * The EXPDE Framework
 * 
 *
 ***************************************************************************/

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  Inform msg2all(argv[0],INFORM_ALL_NODES);
  
  testmsg << "Test" << endl;
  msg2all << "From all nodes " << endl;
  return 0;
}
/***************************************************************************
 * $RCSfile: t1.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $ 
 * IPPL_VERSION_ID: $Id: t1.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $
 ***************************************************************************/
