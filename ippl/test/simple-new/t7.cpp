/***************************************************************************
 *
 * 
 * 
 *
 ***************************************************************************/

// include files
#include "Ippl.h"

// set dimensionality and problem size
const unsigned Dim = 3;

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Vektor<double,Dim> boxMin(-1.0,-1.0,-1.0);
  Vektor<double,Dim> boxMax( 1.0, 1.0, 1.0);
  
  Vektor<double,Dim> p1( .5, .5, .5);
  Vektor<double,Dim> p2( 1.5, 1.5, 1.5);

  if ( dot(p1,p1) <=  dot(boxMax,boxMax))
    testmsg << p1 << " is in the box" << endl;
  else
    testmsg << p1 << " is out of the box" << endl;
  
  if ( dot(p2,p2) <=  dot(boxMax,boxMax))
    testmsg << p2 << " is in the box" << endl;
  else
    testmsg << p2 << " is out of the box" << endl;
  
  return 0;
}

/***************************************************************************
 * $RCSfile: t7.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t7.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $ 
 ***************************************************************************/
