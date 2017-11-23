/*!
 * @file testLayout.cpp
 * @author Matthias Frey
 * @date 2. July 2017
 * @brief Test to check the particle level counter
 * @details In this example we initialize particles
 *          randomly and destroy again.
 */
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>

#include <AMReX_ParmParse.H>

#include "../Solver.h"
#include "../MGTSolver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include "../writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"
#include <random>

#include <getopt.h>

using namespace amrex;

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

struct param_t {
    Vektor<size_t, 3> nr;
    size_t nLevels;
    size_t maxBoxSize;
    double length;
    size_t nParticles;
    bool isHelp;
    AmrOpal::TaggingCriteria criteria;
    double tagfactor;
};

void printHelp(Inform& msg, char* argv[]) {
    msg << "Usage: " << argv[0]
        << endl
        << "--gridx [#gridpoints in x]" << endl
        << "--gridy [#gridpoints in y]" << endl
        << "--gridz [#gridpoints in z]" << endl
        << "--level [#levels]" << endl
        << "--maxgrid [max. grid]" << endl
        << "--boxlength [cube side length]" << endl
        << "--nparticles [#particles >= #procs * 4]" << endl
        << "--tagging charge (default) / efield / potential (optional)" << endl
        << "--tagfactor [charge value / 0 ... 1] (optiona)" << endl;
}


bool parseProgOptions(int argc, char* argv[], param_t& params, Inform& msg) {
    /* Parsing Command Line Arguments
     * 
     * 26. June 2017
     * https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
     */
    
    params.isHelp = false;
    params.criteria = AmrOpal::kChargeDensity;
    params.tagfactor = 1.0e-14; 
    
    int c = 0;
    
    int cnt = 0;
    
    int required = 7;
    
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
            { 0,                0,                 0,  0  }
        };
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "x:y:z:l:m:b:n:ht:f:", long_options, &option_index);
        
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
                params.nParticles = std::atoi(optarg);
                
                if ( (int)params.nParticles < Ippl::getNodes() * 4 ) {
                    printHelp(msg, argv);
                    params.isHelp = true;
                    break;
                }
                
                ++cnt;
                break;
            case 't':
                if ( std::strcmp("efield", optarg) == 0 )
                    params.criteria = AmrOpal::kEfieldStrength;
                else if ( std::strcmp("potential", optarg) == 0 )
                    params.criteria = AmrOpal::kPotentialStrength;
                break;
            case 'f':
                params.tagfactor = std::atof(optarg);
                break;
            case 'h':
                printHelp(msg, argv);
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
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    double halflength = 0.5 * params.length;
    
    std::array<double, AMREX_SPACEDIM> lower = {{-halflength, -halflength, -halflength}}; // m
    std::array<double, AMREX_SPACEDIM> upper = {{ halflength,  halflength,  halflength}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, params.nr, lower, upper);
    
    msg << "Domain: " << domain << endl;
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(params.maxBoxSize));
    
    Array<int> error_buf(params.nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    ParmParse pgeom("geometry");
    Array<int> is_per = { 0, 0, 0};
    pgeom.addarr("is_periodic", is_per);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = params.nr[i];
    
    
    std::vector<int> rr(params.nLevels);
    Array<int> rrr(params.nLevels);
    for (unsigned int i = 0; i < params.nLevels; ++i) {
        rr[i] = 2;
        rrr[i] = 2;
    }
    
    AmrOpal myAmrOpal(&domain, params.nLevels - 1, nCells, 0 /* cartesian */, rr);
    
    myAmrOpal.setTagging(params.criteria);
    
    if ( params.criteria == AmrOpal::kChargeDensity )
        myAmrOpal.setCharge(params.tagfactor);
    else
        myAmrOpal.setScalingFactor(params.tagfactor);
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rrr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    bunch->setAllowParticlesNearBoundary(true);
    
    // initialize a particle distribution
    int nLocParticles = params.nParticles / Ippl::getNodes();
    bunch->create(nLocParticles);
    std::mt19937_64 eng[3];
    for (int i = 0; i < 3; ++i) {
        eng[i].seed(42 + 3 * i);
        eng[i].discard( nLocParticles * Ippl::myNode());
    }
    
    std::uniform_real_distribution<> dist(-halflength, halflength);
    for (uint i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->R[i] = Vector_t( dist(eng[0]), dist(eng[1]), dist(eng[2]) );    // m
        bunch->qm[i] = 1.0e-9;   // C
    }
    
    msg << "==============================================" << endl
        << "                BEFORE UPDATE" << endl;
    auto LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    

    // map particles
    double scale = 1.0;
    
    scale = domainMapping(*bunch, scale);
    
    bunch->update();
    
    msg << "                AFTER UPDATE" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    
    bunch->update();
    
    msg << "             Remove a particle per core " << endl;
    bunch->destroy(1, 0, true);
    bunch->update();
    
    msg << "                AFTER REMOVAL" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    msg << "==============================================" << endl;
    
    
    msg << "             Create a particle per core " << endl;
    bunch->create(1);
    
    for (uint i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->R[i] = Vector_t( dist(eng[0]), dist(eng[1]), dist(eng[2]) );    // m
        bunch->qm[i] = 1.0e-9;   // C
    }
    
    bunch->update();
    
    msg << "                AFTER CREATION" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    msg << "==============================================" << endl;
    
    msg << "#Particles: " << bunch->getTotalNum() << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << bunch->getTotalNum() * bunch->qm[0] << " C" << endl;
    
    myAmrOpal.setBunch(bunch.get());
    
    msg << "//\n//  Process Statistic BEFORE regrid\n//" << endl;
    
    bunch->gatherStatistics();
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, scale/*0.0*/ /*time*/);
    
    bunch->gatherLevelStatistics();
    
    msg << "//\n//  Process Statistic AFTER regrid\n//" << endl;
    
    msg << "==============================================" << endl
        << "                LocalNumPerLevel" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    
    
    msg << "             Remove two particles per core " << endl;
    bunch->destroy(2, 0, true);
    bunch->update();
    
    msg << "                AFTER REMOVAL" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    msg << "==============================================" << endl;
    
    
    msg << "             Create a particle per core " << endl;
    bunch->create(1);
    
    for (uint i = 0; i < bunch->getLocalNum(); ++i) {
        bunch->R[i] = Vector_t( dist(eng[0]), dist(eng[1]), dist(eng[2]) );    // m
        bunch->qm[i] = 1.0e-9;   // C
    }
    
    bunch->update();
    
    msg << "                AFTER CREATION" << endl;
    LocalNumPerLevel = bunch->getLocalNumPerLevel();
    for (int rank = 0; rank < Ippl::getNodes(); ++rank) {
        if ( rank == Ippl::myNode() ) {
            for (uint i = 0; i < LocalNumPerLevel.size(); ++i)
                std::cout << "Rank " << rank
                          << " Level " << i << ": "
                          << LocalNumPerLevel[i] << std::endl;
        }
        Ippl::Comm->barrier();
    }
    msg << "==============================================" << endl;
    
    bunch->gatherStatistics();
    
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
    
        msg << "Particle test running with" << endl
            << "- grid                  = " << params.nr << endl
            << "- max. grid             = " << params.maxBoxSize << endl
            << "- #level                = " << params.nLevels - 1 << endl
            << "- cube side length [m]  = " << params.length << endl
            << "- #particles            = " << params.nParticles << endl
            << "- tagging               = " << tagging << endl
            << "- tagging factor        = " << params.tagfactor << endl;

        doAMReX(params, msg);
        
    } catch(const std::exception& ex) {
        msg << ex.what() << endl;
    }
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();

    return 0;
}
