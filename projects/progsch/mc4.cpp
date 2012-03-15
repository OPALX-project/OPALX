/***************************************************************************

Units:

R, position in units of mega paresec
V, velocity in units of redshift
F, force fiels in units of ???

Field layout (3D domain decomposition)

parallel  I,J,K

In consequence to look over K use the local domain, loops over
I and J can go over the "global" domain

gDomain_m = layout_m->getDomain();        // global domain
lDomain_m = layout_m->getLocalNDIndex();  // local domain

Usage: mc4 <inputfile> [<distrInputFile>]

mpirun -np 2 mc4 <inputfile> [<distrInputFile>] --commlib mpi

And in real:

mpirun -np 1 mc4 indat fort.66 --commlib mpi

And for lesser output try:

mpirun -np 1 mc4 indat fort.66 --info 0 --commlib mpi

 ***************************************************************************/

#include "ChargedParticles/ChargedParticles.hh"

#include <ctime>

#include "FieldSolvers/FFTPoissonSolverPeriodic.hh"
#include <iostream>
#include <sstream>

#include "Parser/FileStream.h"
#include "SimData.hh"

#include "config.h"

/**
   Initializer stuff
*/

#include "IPPLInitializer.h"
#include "Cosmology.h"
#include "InputParser.h"

#ifdef MC4HALOFINDER
/**
   Halo Finder stuff"
*/
#include "CosmoHaloFinderP.h"
#include "HaloProperties.h"
#endif

using namespace std;

#ifdef IPPL_USE_SINGLE_PRECISION
#define T float
#else
#define T double
#endif

/**
   Defines NMAXPERCORE how many particles 
   will be breated at on time i.e. until
   a update() takes place. With this I 
   hope to reduce the load on the MPI
   buffers. 
*/

#define NMAXPERCORE 500000

/************************************/
const string version = "1.1";
/************************************/

const unsigned int DimA = 3;


typedef ChargedParticles<T,DimA>::Mesh_t Mesh_t;
typedef ChargedParticles<T,DimA>::Layout_t Layout_t;
typedef ChargedParticles<T,DimA>::Vector_t Vector_t;  
typedef ChargedParticles<T,DimA>::FieldLayout_t FieldLayout_t;

#ifdef MC4HALOFINDER
void findHalo(ChargedParticles<TT,3> *univ, SimData<T,DimA> simData) 
{
    POSVEL_T boxSize = (POSVEL_T) simData.rL;
    POSVEL_T hubbleConstant = (POSVEL_T) simData.hubble;
    POSVEL_T rL = boxSize / hubbleConstant;

    // Physical coordinate dead zone area   
    POSVEL_T deadSize = (POSVEL_T) (5.);
    
    // Superimposed grid on physical box used to determine wraparound
    int np = simData.ng_comp;
    
    // BB parameter for distance between particles comprising a halo
    POSVEL_T bb = (POSVEL_T) (0.3);
    
    // Minimum number of particles to make a halo
    int pmin = 10;

    string hOutFile("myHalo.out");

    CosmoHaloFinderP haloFinder(univ);
    haloFinder.setParameters(string(hOutFile), rL, deadSize, np, pmin, bb);
    
    /** 
	haloFinder.setParticles(univ);
	 
	At this point the particles which are
	alive status==-1 are in the particle container i.e.
	in univ. This should works with np=1
	 
	ToDO: What is left is a gather of the particles from
	neighbor pe's, and add them to the particle container 
	WITHOUT doing and UPDATE 
	and do a propper haloFinder.setParticles.
	 
    */     
    // set interact radius to 5Mpc
    double radius = 5.0/simData.rL*simData.ng_comp;
    univ->getLayout().setInteractionRadius(radius);
    // ensure we have ghost particles ready
    univ->getLayout().getCashedParticles(*univ);
    haloFinder.setParticles(univ);
     
    haloFinder.executeHaloFinder();
    // Collect the serial halo finder results
    haloFinder.collectHalos();
     
    /** Merge the halos so that only one copy of each is written 
	Parallel halo finder must consult with each of the 26 possible neighbor 
	halo finders to see who will report a particular halo
    */                                                                                               
    haloFinder.mergeHalos();

    // Write halo results on each processor
    haloFinder.writeHalos();
}
#endif


int main(int argc, char *argv[]) {

    Ippl ippl(argc, argv);
    Inform msg ("mc4 ");

    return 0;
	
}

