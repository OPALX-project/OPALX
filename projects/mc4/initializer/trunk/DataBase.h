/*
   Initializer:
   DataBase.h
 
      Defines GlobalStuff class. (main makes a globaly visible 
      instance of that class -- DataBase.) 
      This class will hold simulation 
      and cosmology parameters needed by most of routines. 
 
                     Zarija Lukic, February 2009
                          zarija@lanl.gov
 */

#ifndef DataBase_Header_Included
#define DataBase_Header_Included

#include "TypesAndDefs.h"

#ifdef INITIALIZER_NAMESPACE
namespace initializer {
#endif

class ParallelTools;
class Basedata;

class GlobalStuff {
public :
	void GetParameters(Basedata& bdata, ParallelTools& Parallel);

	// Simulation parameters:
	integer ngrid;       // Number of mesh points in any direction
	real box_size;       // Lenght of the box size in Mpc/h
	integer dim;         // dimensionality of domain decomposition
	unsigned long seed;  // Master seed for random numbers
	real z_in;           // Starting redshift
	
	// Cosmological parameters:
	real Omega_m;    // Total matter fraction (today)
	real Omega_bar;  // Baryon fraction (today)
	real Hubble;     // Hubble constant
	real Sigma_8;    // Mass RMS in 8Mpc/h spheres
	real n_s;        // Spectral index
	real w_de;       // Dark energy EOS parameter
	integer TFFlag;  /* Transfer function used: 0=CMBFast, 1=KH, 
						2=HS, 3=PD, 4=BBKS */
	
	// Code stuff:
	integer PrintFlag;  /* 0 - no output, 1 - serial ASCII, 
						   2 - serial binary, 3 - parallel binary */
};

#ifdef INITIALIZER_NAMESPACE
}
#endif

#endif
