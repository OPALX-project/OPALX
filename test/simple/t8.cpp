/***************************************************************************
 *
 * The EXPDE Framework
 * 
 *
 ***************************************************************************/
#include "Ippl.h"

// set dimensionality and problem size
const unsigned Dim = 3;

int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);

  Vektor<double,Dim> boxMin(atof(argv[1]),atof(argv[2]),atof(argv[3]));
  Vektor<double,Dim> boxMax(atof(argv[4]),atof(argv[5]),atof(argv[6]));
  Vektor<double,Dim> h(1.0e-2, 1.0e-2, 0.5e-2);

  Vektor<double,Dim> p1(atof(argv[7]),atof(argv[8]),atof(argv[9]));
  
  Vektor<unsigned,Dim> N(static_cast<unsigned> (ceil( (abs(boxMin[0])+boxMax[0])/h[0])),
		       static_cast<unsigned> (ceil( (abs(boxMin[1])+boxMax[1])/h[1])),
		       static_cast<unsigned> (ceil( (abs(boxMin[2])+boxMax[2])/h[2])));

  testmsg << "orig= " << boxMin << " maxext= " << boxMax << endl;
  testmsg << "h=  " << h << " N= " << N << endl;
  testmsg << "p1= " << p1 << endl;
  
  /*
  p1 -> n
  */

  testmsg << "n(p1)= " << -(boxMin-p1)/h << endl;
  return 0;
}

/***************************************************************************
 * $RCSfile: t8.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/03/28 12:47:26 $
 * IPPL_VERSION_ID: $Id: t8.cpp,v 1.1.1.1 2003/03/28 12:47:26 adelmann Exp $ 
 ***************************************************************************/
