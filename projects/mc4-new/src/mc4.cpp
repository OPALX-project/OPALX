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
#include <ctime>
#include <iostream>
#include <sstream>

#include "ChargedParticles/ChargedParticles.hh"
#include "FieldSolvers/FFTPoissonSolverPeriodic.hh"
#include "Parser/FileStream.h"
#include "SimData.hh"
#include "config.h"

/**
   Initializer stuff
*/

#include "Initializer/IPPLInitializer.h"
#include "Initializer/Cosmology.h"
#include "Initializer/InputParser.h"

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

const string version = PACKAGE_VERSION;

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

    Inform hmsg("");
    string mySpace("            ");

    /*    hmsg << mySpace << "  .___  ___.  ______  _  _    "<< endl;
    hmsg << mySpace << "  |   \/   |  /      || || |   "<< endl;
    hmsg << mySpace << "  |  \  /  | |  ,----'| || |_  "<< endl;
    hmsg << mySpace << "  |  |\/|  | |  |     |__   _| "<< endl;
    hmsg << mySpace << "  |  |  |  | |  `----.   | |   "<< endl;
    hmsg << mySpace << "  |__|  |__|  \______|   |_|   "<< endl;
    */
    hmsg << "This is MC 4 Version " << PACKAGE_VERSION << endl;
    hmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) to " << PACKAGE_BUGREPORT << endl;

#ifdef IPPL_USE_SINGLE_PRECISION
    msg << "Running in SINGLE precision mode" << endl;
#else
    msg << "Running in DOUBLE precision mode" << endl;
#endif

    const T Pi = 4.0*atan(1.0);
    IpplTimings::TimerRef TWall = IpplTimings::getTimer("Wall");
    IpplTimings::TimerRef TCDist = IpplTimings::getTimer("Create Distribution");
    IpplTimings::TimerRef TStep = IpplTimings::getTimer("Step Total");
    IpplTimings::TimerRef TPwrSpec = IpplTimings::getTimer("PowerSpectra");

    const int memoryPerNode = 2; // GB
    static IpplMemoryUsage::MemRef mainMemWatch = IpplMemoryUsage::getMemObserver("main", memoryPerNode);
    static IpplMemoryUsage::MemRef initiMemWatch = IpplMemoryUsage::getMemObserver("beforePartCrea", memoryPerNode);
    static IpplMemoryUsage::MemRef initfMemWatch = IpplMemoryUsage::getMemObserver("afterPartCrea ", memoryPerNode);

    static IpplMemoryUsage::MemRef afterSolverInitWatch = IpplMemoryUsage::getMemObserver("afterSolverInit", memoryPerNode);     
    static IpplMemoryUsage::MemRef afterPwrSpecWatch = IpplMemoryUsage::getMemObserver("afterPwrSpec", memoryPerNode);

 
    IpplMemoryUsage::sample(mainMemWatch,"");
    IpplTimings::startTimer(TWall);
    size_t chunksize = 0;

    //handle chunking
    if(argc > 4) {
        chunksize = atoi(argv[4]);
        msg << "Chunking enabled with size = " << chunksize << endl;
    }

    /**
       univ is used in the calculation
    */
    ChargedParticles<T,DimA>* univ;

    /// here we read in all data from the input
    /// file and compute simple derifed quantities

    SimData<T,DimA> simData(argc, argv);
    msg << simData << endl;
    const int nginit = simData.ng_dist;

    // create layout objects
    // Changed from ng to np by LD - right for initial conditions, wrong for force calc.

    Index I(nginit), J(nginit), K(nginit);
    NDIndex<DimA> domainI(I,J,K); //Initial domain

    const int ng = simData.ng_comp;
    Index A(ng), B(ng), C(ng);
    NDIndex<DimA> domainF(A,B,C); //Force solver domain

    e_dim_tag decomp[DimA];
    for (int d=0; d < DimA; ++d)
       	decomp[d] = PARALLEL;

    if(decomp[0] == PARALLEL && decomp[1] == PARALLEL && decomp[2] == PARALLEL)
	msg << "Using 3D domain decomposition" << endl;
    else
	msg << "Using 1D or 2D  domain decomposition" << endl;
    
    // create mesh and layout objects for this problem domain
    Mesh_t *meshI, *meshF;
    FieldLayout_t *FLI, *FLF;
    meshI = new Mesh_t(domainI);
    meshF = new Mesh_t(domainF);
    FLI = new FieldLayout_t(*meshI, decomp);
    FLF = new FieldLayout_t(*meshF, decomp);
    Layout_t * PLI = new Layout_t(*FLI, *meshI);  
    Layout_t * PLF = new Layout_t(*FLF, *meshF);

    /**
       simData.jinit == 0 (generate)
       simData.jinit == 1 ASCII (fort.66)
       simData.jinit == 3 H5Part
       simData.jinit == 4 Uniform 
    */

    IpplTimings::startTimer(TCDist);
    IpplMemoryUsage::sample(initiMemWatch,"");

    if ((simData.jinit == 1) || (simData.jinit == 3)) {
        univ = new ChargedParticles<T,DimA>(PLF,simData,Vector_t(1.0*simData.ng_comp));
        univ->readInputDistribution(string(argv[2]));   
    }
    else if (simData.jinit == 0) {
        //use initialize wrapper to create universe
        string indatName = argv[1];
        string tfName = argv[2];

        initializer::InputParser bdata("init.in");

        size_t Npart;
        int NumProcs, rank, i;

        NumProcs = Ippl::Comm->getNodes();
        rank = Ippl::Comm->myNode();

        //DETERMINE HOWMANY PROCS TO USE
        int np = 0;
        bdata.getByName("np", np);
        
        ParticleAttrib<Vector_t> pos;
        ParticleAttrib<Vector_t> vel;
        ParticleAttrib<T> mass;

        Npart = (np*np*np)/NumProcs;

	pos.create(Npart);
	vel.create(Npart);
	mass.create(Npart);
	initializer::IPPLInitializer init;
	initializer::InputParser par("init.in");
	init.init_particles(pos, vel, FLI, meshI, par, tfName.c_str(), MPI_COMM_WORLD);
	mass = 1.0;
        	
	size_t rest = 0;
	int iterations = 0;
	int tag = IPPL_APP_TAG1;
	int parent = 0;

	if( rank == parent && Npart != 0 ) {
	    rest = Npart % chunksize;
	    iterations = (Npart - rest) / chunksize;

	    Message *mess1 = new Message();
	    putMessage( *mess1, rest );
	    putMessage( *mess1, iterations );
	    Ippl::Comm->broadcast_all(mess1, tag);
	}

	Message *mess = Ippl::Comm->receive_block(parent, tag);
	PAssert(mess);
	getMessage( *mess, rest );
	getMessage( *mess, iterations );
	if (mess)
	  delete mess;

	msg << "about to copy ics over to univ junksize= " << chunksize << " iterations " << iterations << " rest " << rest << endl;
        //pos, vel and mass particle attributes are destroyed in copy constructor
        univ = new ChargedParticles<T,DimA>(PLF, simData, pos, vel, mass, Npart, Vector_t(simData.ng_comp), chunksize, iterations, rest);
        msg << "done coping ics to univ" << endl;	 	 
        
        //XXX; make sure R is scaled from [0, ng_comp]
        univ->R *= Vector_t((double)simData.ng_comp/simData.np);

        //if (simData.doNeutrinos) 
        //injectNeutrinos(bdata,simData,usedNumProcs,tfName,gridgencomm,univ);

    }
    else if (simData.jinit == 4) {
        /**
          univ0 is used to generate uniform initial conditions
          for testing only
          */
        ChargedParticles<T,DimA>* univ0 = new ChargedParticles<T,DimA>(PLI,simData);
        univ0->set_FieldSolver(new FFTPoissonSolverPeriodic<T,DimA>(univ0,simData));
        univ0->setBoxSize(Vector_t(1.0,1.0,1.0));
        //T gpscal = (1.0*simData.np);
        T gpscal = (1.0*simData.ng_comp);
        univ0->genLocalUniformICiter(gpscal);
        
        // create an empty ChargedParticles object and set the BC's
        univ = univ0;
            //new ChargedParticles<T,DimA>(PLF,simData,
                //univ0->R,
                //univ0->V,
                //univ0->M,
                //univ0->getLocalNum(),
                //Vector_t(simData.ng_comp),
                //chunksize);
    }

    IpplTimings::stopTimer(TCDist);
    IpplMemoryUsage::sample(initfMemWatch,"");

    /*
      FIXME: At the moment deleting univ0 causes a core dump, no idea why 
      delete univ0;
      delete FLI;
      delete meshI;
      delete PLI;
    */

    T ain  = 1.0 /(1+simData.zin);
    T ainv = 1.0 /(simData.alpha);
    T pp   = pow(ain,simData.alpha);
     
    //     DataSink *itsDataSink = NULL;
    DataSink *itsDataSink = new DataSink("universe.h5", univ);
    if (simData.iprint < 99)
        itsDataSink->writePhaseSpace(pow(pp, ainv), 1.0/pow(pp,ainv)-1.0, 0);

    univ->setTime(0.0);  
    univ->set_FieldSolver(new FFTPoissonSolverPeriodic<T,DimA>(univ,simData));
    msg << "field solver ready " << endl;
    IpplMemoryUsage::sample(afterSolverInitWatch,"");

    T omegab   = simData.deut / simData.hubble / simData.hubble;
    T omegatot = simData.omegadm + omegab;

    /// Convert to `symplectic velocities' at a=a_in 
    univ->V *= Vector_t(ain*ain,ain*ain,ain*ain);
    
    size_t localNp = univ->getLocalNum();
    reduce(localNp, localNp, OpAddAssign());
    msg << "Sanity check: total Np= " << localNp << endl; 

/** 
	Initial dump
    */
    if ((simData.np < 128) && (Ippl::getNodes()==1))
	univ->dumpParticles(string("mc4-out/step-univ"));      
     
    /** 
	Initial P(k)
    */
    IpplTimings::startTimer(TPwrSpec);
    msg << "Calculate power spectra of initial step "  << endl;
    univ->fs_m->calcPwrSpecAndSave(univ,string("mc4-1D-initial.spec"));
    IpplTimings::stopTimer(TPwrSpec);

    IpplMemoryUsage::sample(afterPwrSpecWatch,"");

    // timesteps loop
    IpplTimings::startTimer(TStep);
    univ->tStep(itsDataSink);
    IpplTimings::stopTimer(TStep);
     
    //if (Ippl::getNodes()==1)
	//findHalo(univ,simData);

    IpplTimings::startTimer(TPwrSpec);
    msg << "Calculate power spectra of last step "  << endl;
    univ->fs_m->calcPwrSpecAndSave(univ,string("mc4-1D.spec"));
    //if (Ippl::getNodes()==1)
    //univ->calcAndDumpCorrelation(63);
    IpplTimings::stopTimer(TPwrSpec);

    T pmass = 2.77536627e11*std::pow(simData.rL,(T)3.0)*(omegatot)/simData.hubble/univ->getTotalNum();

    msg << "Total Mass = " << pmass << endl;

    msg << "All nodes done dude ..." << endl;

    IpplTimings::stopTimer(TWall);
    IpplTimings::print();
    IpplTimings::print(string("timing.dat"));
    IpplMemoryUsage::print();
    IpplMemoryUsage::print(string("memory.dat"));

    //if(Ippl::getNodes() == 1)
    //    univ->dumpCosmo(10);

    if(univ)
	delete univ;
     
    if (itsDataSink)
	delete itsDataSink;

    Ippl::Comm->barrier();
    return 0;
	
}

