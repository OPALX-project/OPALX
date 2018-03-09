/***************************************************************************
 *
 * 
 * 
 *
 ***************************************************************************/

#include "Ippl.h"

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  
  // random number check
  
  IpplRandom.SetSeed(static_cast<unsigned long>(23131719));
  double dummy = IpplRandom();
  unsigned int N=Ippl::getNodes();
  unsigned int K = 1000;
  
  for (unsigned i=0; i<K; i++) {
    double localData = IpplRandom();
    double globalRes = 0.0;
    
    reduce(localData, globalRes, OpAddAssign());
    
    if (abs(globalRes-(N*localData)) > 1e-9) {
      ERRORMSG("Parallel random check failed: abs(globalRes-(N*localData)) > 1e-9 " << endl);
    }   
  }
  return 0;
}
/***************************************************************************
 * $RCSfile: t4.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t4.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $ 
 ***************************************************************************/
