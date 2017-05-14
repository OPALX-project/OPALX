/*!
 * @file testTagging.cpp
 * @author Matthias Frey
 * @date February 2017
 * 
 * Domain:  [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: 1e-10\n
 * Gaussian particle distribution N(0.0, 0.07)
 * 
 * @details Test the different tagging schemes.
 * 
 * Call:\n
 *  mpirun -np [#cores] testTagging [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                     [#particles] [#levels] [max. box size] [tagging]
 * 
 * @brief Perturbes particles randomly in space for several time steps.
 */

#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>


#include <ParmParse.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;

void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, AmrOpal::TaggingCriteria criteria, double factor, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("regrid");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, BL_SPACEDIM> lower = {{-0.5, -0.5, -0.5}}; // m
    std::array<double, BL_SPACEDIM> upper = {{ 0.5,  0.5,  0.5}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    Array<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    unsigned long int nloc = nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.gaussian(0.0, 0.07, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam( *(bunch.get()) );
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->qm[i] = 1.0e-10;  // in [C]
    }
    
    // redistribute on single-level
    bunch->update();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl;
    
    
    myAmrOpal.setBunch(bunch.get());
    
    const Array<Geometry>& geoms = myAmrOpal.Geom();
    
    // tagging using potential strength
    myAmrOpal.setTagging(criteria);
    
    if ( criteria == AmrOpal::kChargeDensity )
        myAmrOpal.setCharge(factor);
    else
        myAmrOpal.setScalingFactor(factor);
    
    IpplTimings::startTimer(regridTimer);
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    IpplTimings::stopTimer(regridTimer);
    
//     // single core
//     unsigned int lev = 0, nLocParticles = 0;
//     for (unsigned int ip = 0; ip < bunch->getLocalNum(); ++ip) {
//         while ( lev != bunch->m_lev[ip] ) {
//             std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
//             nLocParticles = 0;
//             lev++;
//         }
//         nLocParticles++;
//     }
//     
//     std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
    
    bunch->python_format(0);
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] "
         << "[#levels] [max. box size] [tagging (i.e. charge, efield, potential)] [factor (charge or 0..1]";
    
    if ( argc < 9 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    size_t nParticles = std::atoi(argv[4]);
    std::string tagging = argv[7];
    double factor = std::atof(argv[8]);
    
    AmrOpal::TaggingCriteria criteria = AmrOpal::kChargeDensity;
    if ( !tagging.compare("efield") )
        criteria = AmrOpal::kEfieldGradient;
    else if ( !tagging.compare("potential") )
        criteria = AmrOpal::kPotentialStrength;
    else
        tagging = "charge"; // take default method: kChargeDensity
    
    msg << "Particle test running with" << endl
        << "- #particles     = " << nParticles << endl
        << "- grid           = " << nr << endl
        << "- tagging        = " << tagging << endl
        << "- charge/scaling = " << factor << endl;
    
    BoxLib::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    doBoxLib(nr, nParticles, nLevels, maxBoxSize, criteria, factor, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
