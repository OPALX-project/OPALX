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
 *
 ***************************************************************************/

#include "Ippl.h"
#include "Particle/ParticleAttrib.h"
#include "AppTypes/Vektor.h"
#include <vector>
using namespace std;


typedef ParticleAttrib<Vektor<double,3> > ParticlePos_t;
typedef ParticleAttrib<unsigned>       ParticleIndex_t;

typedef vector<ParticleAttribBase *>      attrib_container_t;


int main()
{

ParticlePos_t R;

attrib_container_t aContainer;


aContainer.push_back(&R);



	return 0;
}

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

