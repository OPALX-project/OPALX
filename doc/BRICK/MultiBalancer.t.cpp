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
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "MultiBalancer.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void multiBalancerTester()
{
  const int procs = 8;
  const int vnodes = 32;
  const int materials = 10;
  int m,v,p;
  double weights[materials][vnodes];
  
  MultiBalancer balancer(procs,vnodes);

  for (m=0; m<materials; ++m)
    {
      for (v=0; v<vnodes; ++v)
        weights[m][v] = ( rand()<(RAND_MAX/4) ? 1.0 : 0.0 );
      balancer.newMaterial();
      balancer.appendWeights(weights[m],weights[m]+vnodes);
    }
  balancer.distribute();
  
  printf("\nInput weights:\n");
  for (v=0; v<vnodes; ++v)
    printf("%7d",v);
  printf("\n");
  for (m=0; m<materials; ++m)
    {
      for (v=0; v<vnodes; ++v)
        printf("%7.3f",balancer.m_inputWeights[m][0][v]);
      printf("\n");
    }
  
  printf("\nProcessor destinations:\n");
  for (v=0; v<vnodes; ++v)
    printf("%7d",balancer.m_vnodeProcs[v]);
  printf("\n");

  printf("\nProcessor weights:\n");
  for (m=0; m<materials; ++m)
    {
      for (p=0; p<procs; ++p)
        printf("%7.3f",balancer.m_materialProcWeights[m][0][p]);
      printf("\n");
    }
  printf("Totals:\n");
  for (p=0; p<procs; ++p)
    printf("%7.3f",balancer.m_procWeights[p]);
  printf("\n");

  printf("\nMaterial Max weights:\n");
  for (m=0; m<materials; ++m)
    printf("%7.3f",balancer.m_materialMaxWeights[m]);
  printf("\n");
}

int main()
{
  multiBalancerTester();
  return 0;
}


/***************************************************************************
 * $RCSfile: MultiBalancer.t.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:27 $
 * IPPL_VERSION_ID: $Id: MultiBalancer.t.cpp,v 1.1.1.1 2003/01/23 07:40:27 adelmann Exp $ 
 ***************************************************************************/
