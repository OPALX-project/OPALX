/*!
 * @file testTagging.cpp
 * @author Matthias Frey
 * @date February 2017
 * 
 * Domain:  [-1, 1] x [-1, 1] x [-1, 1]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: 1e-10\n
 * Gaussian particle distribution N(0.0, 0.25)
 * 
 * @details Test the different tagging schemes.
 * 
 * Call:\n
 *  ./testTagging --help
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


#include <AMReX_ParmParse.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

#include <getopt.h>

using namespace amrex;

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    size_t nParticles;
    AmrOpal::TaggingCriteria criteria;
    double tagfactor;
    bool isHelp;
    double sigma;
    double length;
};

bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isHelp = false;
    params.criteria = AmrOpal::kChargeDensity;
    params.tagfactor = 1.0e-14;
    params.sigma = 0.25;
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 8;
    
    while ( true ) {
        static struct option long_options[] = {
            { "gridx",          required_argument, 0, 'x' },
            { "gridy",          required_argument, 0, 'y' },
            { "gridz",          required_argument, 0, 'z' },
            { "level",          required_argument, 0, 'l' },
            { "maxgrid",        required_argument, 0, 'm' },
            { "boxlength",      required_argument, 0, 'b' },
            { "nparticles",     required_argument, 0, 'n' },
            { "help",           no_argument,       0, 'h' },
            { "tagging",        required_argument, 0, 't' },
            { "tagging-factor", required_argument, 0, 'f' },
            { "sigma",          required_argument, 0, 's' },
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "x:y:z:l:m:r:b:n:ht:f:s:", long_options, &option_index);
        
        if ( c == -1 )
            break;
        
        switch ( c ) {
            case 'x':
                params.nr[0] = std::atoi(optarg); ++cnt; break;
            case 'y':
                params.nr[1] = std::atoi(optarg); ++cnt; break;
            case 'z':
                params.nr[2] = std::atoi(optarg); ++cnt; break;
            case 'l':
                params.nLevels = std::atoi(optarg) + 1; ++cnt; break;
            case 'm':
                params.maxBoxSize = std::atoi(optarg); ++cnt; break;
            case 'b':
                params.length = std::atof(optarg); ++cnt; break;
            case 'n':
                params.nParticles = std::atoi(optarg); ++cnt; break;
            case 't':
                if ( std::strcmp("efield", optarg) == 0 )
                    params.criteria = AmrOpal::kEfieldStrength;
                else if ( std::strcmp("potential", optarg) == 0 )
                    params.criteria = AmrOpal::kPotentialStrength;
                else if ( std::strcmp("min-npart", optarg) == 0 )
                    params.criteria = AmrOpal::kMinNumParticles;
                else if ( std::strcmp("max-npart", optarg) == 0 )
                    params.criteria = AmrOpal::kMaxNumParticles;
                else if ( std::strcmp("momentum", optarg) == 0 )
                    params.criteria = AmrOpal::kMomentum;
                break;
            case 'f':
                params.tagfactor = std::atof(optarg);
                break;
            case 's':
                params.sigma = std::atof(optarg); ++cnt; break;
            case 'h':
                msg << "Usage: " << argv[0]
                    << endl
                    << "--gridx [#gridpoints in x]" << endl
                    << "--gridy [#gridpoints in y]" << endl
                    << "--gridz [#gridpoints in z]" << endl
                    << "--level [#levels]" << endl
                    << "--maxgrid [max. grid]" << endl
                    << "--nparticles [#particles]" << endl
                    << "--boxlength [cube side length]" << endl
                    << "--sigma [Gaussian sigma]" << endl
                    << "--tagging charge (default) / efield / potential / min-npart (optional)" << endl
                    << "--tagging-factor [charge value / 0 ... 1] (optional)" << endl;
                params.isHelp = true;
                break;
            case '?':
                break;
            
            default:
                break;
            
        }
    }

    return ( cnt == required );
}

void doAMReX(const param_t& params, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("regrid");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double halflength = 0.5 * params.length;
    
    std::array<double, AMREX_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
    std::array<double, AMREX_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, params.nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    Vector<int> error_buf(params.nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    Vector<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    AmrOpal myAmrOpal(&domain, params.nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Vector<BoxArray>& ba = myAmrOpal.boxArray();
    const Vector<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Vector<Geometry>& geom = myAmrOpal.Geom();
    
    Vector<int> rr(params.nLevels);
    for (uint i = 0; i < params.nLevels; ++i)
        rr[i] = 2;
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    unsigned long int nloc = params.nParticles / ParallelDescriptor::NProcs();
    Distribution dist;
    dist.gaussian(0.0, params.sigma, nloc, ParallelDescriptor::MyProc());
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam( *(bunch.get()) );
    
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->qm[i] = 1.0e-10;  // in [C]
        
        // increase momentum with distance from center
        Vector_t r2 = bunch->R[i] * bunch->R[i];
        bunch->P[i] = Vector_t(std::sqrt(r2(0)),
                               std::sqrt(r2(1)),
                               std::sqrt(r2(2)));
    }
    
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    // redistribute on single-level
    bunch->update();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "#Particles: " << params.nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << params.nParticles * bunch->qm[0] << " C" << endl;
    
    
    
    myAmrOpal.setBunch(bunch.get());
    
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    
    // tagging using potential strength
    myAmrOpal.setTagging(params.criteria);
    
    switch ( params.criteria )
    {
        case AmrOpal::kChargeDensity:
            myAmrOpal.setCharge(params.tagfactor); break;
        case AmrOpal::kMinNumParticles:
            myAmrOpal.setMinNumParticles(params.tagfactor); break;
        case AmrOpal::kMaxNumParticles:
            myAmrOpal.setMaxNumParticles(params.tagfactor); break;
        case AmrOpal::kEfieldStrength:
            ;
        case AmrOpal::kPotentialStrength:
            ;
        case AmrOpal::kMomentum:
            myAmrOpal.setScalingFactor(params.tagfactor); break;
        default:
            break;
    }
    
    IpplTimings::startTimer(regridTimer);
    
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, scale /* 0.0 time*/);
    
    IpplTimings::stopTimer(regridTimer);
    
    bunch->python_format(0);
    
    
    domainMapping(*bunch, scale, true);
    
    bunch->update();
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    param_t params;
    
    amrex::Initialize(argc,argv, false);
    
    try {
        if ( !parseProgOptions(argc, argv, params, msg) && !params.isHelp )
            throw std::runtime_error("\033[1;31mError: Check the program options.\033[0m");
        else if ( params.isHelp )
            return 0;
        
        
        std::string tagging = "charge";
        if ( params.criteria == AmrOpal::kEfieldStrength )
            tagging = "efield";
        else if ( params.criteria == AmrOpal::kPotentialStrength )
            tagging = "potential";
        else if ( params.criteria == AmrOpal::kMinNumParticles )
            tagging = "min-npart";
        else if ( params.criteria == AmrOpal::kMaxNumParticles )
            tagging = "max-npart";
        else if ( params.criteria == AmrOpal::kMomentum )
            tagging = "momentum";
    
        msg << "Particle test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- cube side length [m]  = " << params.length << endl
            << "- #particles            = " << params.nParticles << endl
            << "- sigma  [Gaussian]     = " << params.sigma << endl
            << "- tagging               = " << tagging << endl
            << "- tagging factor        = " << params.tagfactor << endl;

        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    return 0;
}
