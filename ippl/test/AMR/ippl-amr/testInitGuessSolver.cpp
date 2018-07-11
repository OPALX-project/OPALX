/*!
 * @file testInitGuessSolver.cpp
 * @author Matthias Frey
 * @date February 2017
 * 
 * Domain:  [-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5]\n
 * BC:      Dirichlet boundary conditions\n
 * Charge/particle: elementary charge\n
 * Gaussian particle distribution N(0.0, 0.01)
 * 
 * @details Compare the solution of the Poisson solver without and with smart initial guess.
 *          It computes the l2-norm of the difference of both solutions.
 *          In every time step each particle is randomly perturbed by a small \f$\delta x\f$.
 * 
 * Call:\n
 *  mpirun -np [#cores] testInitGuessSolver [#gridpoints x] [#gridpoints y] [#gridpoints z]
 *                                          [#particles] [#levels] [max. box size] [#steps] [reuse]
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

#include <Interpolater.H>


#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

#include "../boxlib-amr/writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, AMREX_SPACEDIM> Vector_t;

void precondition(AmrOpal& myAmrOpal,
                  container_t& rhs,
                  container_t& phi,
                  container_t& grad_phi,
                  int nLevels,
                  bool isRegrid)
{
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, grad_phi, myAmrOpal.boxArray()[lev], lev);
    }
    
    if ( isRegrid ) {
        for (int i = myAmrOpal.finestLevel()-1; i >= 0; --i) {
            MultiFab tmp(phi[i].boxArray(), 1, 0, phi[i].DistributionMap());
            tmp.setVal(0.0);
            BoxLib::average_down(phi[i+1], tmp, 0, 1, myAmrOpal.refRatio(i));
            MultiFab::Add(phi[i], tmp, 0, 0, 1, 0);
        }
        
        for (int lev = 1; lev < nLevels; ++lev) {
            phi.clear(lev);
            phi.set(lev, new MultiFab(myAmrOpal.boxArray()[lev],1          ,1));
            phi[lev].setVal(0.0);
        }
        
        PCInterp mapper;
        Vector<BCRec> bc(1);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            bc[0].setLo(i, EXT_DIR);
            bc[0].setHi(i, EXT_DIR);
        }
        
        
        IntVect fine_ratio(2, 2, 2);
        const InterpolaterBoxCoarsener& coarsener = mapper.BoxCoarsener(fine_ratio);
        
        for (int lev = 0; lev < nLevels - 1; ++lev) {
            const BoxArray& ba = phi[lev + 1].boxArray();
            const DistributionMapping& dm = phi[lev + 1].DistributionMap();
            int ngrow = phi[lev + 1].nGrow();
    
            const IndexType& typ = ba.ixType();
    
            BL_ASSERT(typ == phi[lev].boxArray().ixType());
            
            Geometry& fgeom = myAmrOpal.Geom(lev + 1);
            Geometry& cgeom = myAmrOpal.Geom(lev);
    
            Box fdomain = fgeom.Domain();
            fdomain.convert(typ);
    
            Box fdomain_g(fdomain);
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                if (fgeom.isPeriodic(i)) {
                    fdomain_g.grow(i,ngrow);
                }
            }
    
            BoxArray ba_crse_patch(ba.size());
            {  // TODO: later we might want to cache this
                for (int i = 0, N = ba.size(); i < N; ++i)
                {
                    Box bx = BoxLib::convert(BoxLib::grow(ba[i],ngrow), typ);
                    bx &= fdomain_g;
                    ba_crse_patch.set(i, coarsener.doit(bx));
                }
            }
            
            int ncomp = 1;
            int scomp = 0;
            int dcomp = 0;
            MultiFab mf_crse_patch(ba_crse_patch, ncomp, 0, dm);
    
            mf_crse_patch.copy(phi[lev], scomp, 0, ncomp, cgeom.periodicity());
    
//             cbc.FillBoundary(mf_crse_patch, 0, ncomp, time);
    
            int idummy1=0, idummy2=0;
    
#ifdef _OPENMP
            #pragma omp parallel
#endif
            for (MFIter mfi(mf_crse_patch); mfi.isValid(); ++mfi)
            {
                FArrayBox& dfab = phi[lev+1][mfi];
                const Box& dbx = dfab.box() & fdomain_g;

                Vector<BCRec> bcr(ncomp);
                BoxLib::setBC(dbx,fdomain,0,0,1,bc,bcr);
//                 BoxLib::setBC(dbx,fdomain,scomp,0,ncomp,bcs,bcr);
                
                mapper.interp(mf_crse_patch[mfi],
                            0,
                            phi[lev+1][mfi],
                            dcomp,
                            ncomp,
                            dbx,
                            fine_ratio,
                            cgeom,
                            fgeom,
                            bcr,
                            idummy1, idummy2);	    
            }
        }
            
//             for (MFIter fmfi(phi[lev + 1], false); fmfi.isValid(); ++fmfi) {
//                 
//                 const Box& bx = fmfi.fabbox(); //validbox();
//                 FArrayBox& fab = phi[lev + 1][fmfi];
//                 
//                 std::cout << bx << std::endl;
//                 
// //                 FArrayBox finefab(bx, 1);
// //                 FArrayBox crsefab(mapper.CoarseBox(bx, fine_ratio), 1);
// //                 
// //                 for (MFIter cmfi(phi[lev], false); cmfi.isValid(); ++cmfi) {
// //                     crsefab.copy(phi[lev][cmfi], 0, 0);
// //                 }
// //                 
// //                 mapper.interp(crsefab,
// //                               0, // comp
// //                               finefab,
// //                               0, // comp
// //                               1, // ncomp
// //                               bx,
// //                               fine_ratio,
// //                               gcrse,
// //                               gfine,
// //                               bc,
// //                               0,
// //                               0);
// //                 
// // //                 std::cout << "Before: " << fab << std::endl; std::cin.get();
// //                 fab.copy(finefab, 0, 0);
// // //                 break;
// // //                 std::cout << "After: " << fab << std::endl; std::cin.get();
//             }
//         }
    }
    
    Vector<int> rr(nLevels);
    for (int i = 0; i < nLevels; ++i)
        rr[i] = 2;
    std::string plotnormal = BoxLib::Concatenate("init_p", 10, 4);
    writePlotFile(plotnormal, rhs, phi, grad_phi, rr, myAmrOpal.Geom(), 10);
}
    
void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Vector<Geometry>& geom,
             const Vector<int>& rr,
             int nLevels,
             bool reuse,
             bool isRegrid,
             Inform& msg,
             IpplTimings::TimerRef& assignTimer,
             IpplTimings::TimerRef& actualSolveTimer,
             IpplTimings::TimerRef& initGuessTimer)
{
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    if ( reuse ) {
        IpplTimings::startTimer(initGuessTimer);
        precondition(myAmrOpal, rhs, phi, grad_phi, nLevels, isRegrid);
        IpplTimings::stopTimer(initGuessTimer);
    } else {
        IpplTimings::startTimer(initGuessTimer);
        for (int lev = 0; lev < nLevels; ++lev) {
            initGridData(rhs, phi, grad_phi, myAmrOpal.boxArray()[lev], lev);
        }
        IpplTimings::stopTimer(initGuessTimer);
    }

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

    IpplTimings::startTimer(assignTimer);
    bunch->AssignDensity(bunch->qm, false, rhs, base_level, finest_level);
    IpplTimings::stopTimer(assignTimer);
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(constant, 0, 1);
#endif
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;

    // solve                                                                                                                                                                                                     
    Solver sol;
    IpplTimings::startTimer(actualSolveTimer);
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset,
                        false);
    IpplTimings::stopTimer(actualSolveTimer);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
}

void doBoxLib(const Vektor<size_t, 3>& nr, size_t nParticles,
              int nLevels, size_t maxBoxSize, int nSteps, Inform& msg)
{
    static IpplTimings::TimerRef regridTimer = IpplTimings::getTimer("tracking-regrid");
    static IpplTimings::TimerRef solveWithTimer = IpplTimings::getTimer("solve-smart");
    static IpplTimings::TimerRef solveWithoutTimer = IpplTimings::getTimer("solve-normal");
    static IpplTimings::TimerRef assignTimer = IpplTimings::getTimer("assign-charge");
    static IpplTimings::TimerRef initGuessWithTimer = IpplTimings::getTimer("init-solver-guess-smart");
    static IpplTimings::TimerRef initGuessWithoutTimer = IpplTimings::getTimer("init-solver-guess-normal");
    static IpplTimings::TimerRef actualSolveWithTimer = IpplTimings::getTimer("actual-solve-smart");
    static IpplTimings::TimerRef actualSolveWithoutTimer = IpplTimings::getTimer("actual-solve-normal");
    // ========================================================================
    // 1. initialize physical domain (just single-level)
    // ========================================================================
    
    std::array<double, AMREX_SPACEDIM> lower = {{-0.5, -0.5, -0.5}}; // m
    std::array<double, AMREX_SPACEDIM> upper = {{ 0.5,  0.5,  0.5}}; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Vector<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
//     ParmParse pgeom("geometry");
//     Vector<int> is_per = { 1, 1, 1};
//     pgeom.addarr("is_periodic", is_per);
    
    Vector<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Vector<BoxArray>& ba = myAmrOpal.boxArray();
    const Vector<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Vector<Geometry>& geom = myAmrOpal.Geom();
    
    Vector<int> rr(nLevels);
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
        bunch->qm[i] = 1.0e-19;  // in [C]
    }
    
    // redistribute on single-level
    bunch->update();
    
    msg << "Single-level statistics" << endl;
    bunch->gatherStatistics();
    
    msg << "#Particles: " << nParticles << endl
        << "Charge per particle: " << bunch->qm[0] << " C" << endl
        << "Total charge: " << nParticles * bunch->qm[0] << " C" << endl;
    
    
    myAmrOpal.setBunch(bunch.get());
    
    myAmrOpal.setTagging(AmrOpal::kChargeDensity);
    
    const Vector<Geometry>& geoms = myAmrOpal.Geom();
    
    // ========================================================================
    // 3. multi-level redistribute
    // ========================================================================
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    // container for "without smart guess"
    container_t rhs_wo(PArrayManage);
    container_t phi_wo(PArrayManage);
    container_t grad_phi_wo(PArrayManage);
    rhs_wo.resize(nLevels);
    phi_wo.resize(nLevels);
    grad_phi_wo.resize(nLevels);
    
    // container for "with smart guess"
    container_t rhs_w(PArrayManage);
    container_t phi_w(PArrayManage);
    container_t grad_phi_w(PArrayManage);
    rhs_w.resize(nLevels);
    phi_w.resize(nLevels);
    grad_phi_w.resize(nLevels);
    
    std::mt19937_64 eng(0);
    std::uniform_real_distribution<> perturbation(-0.01, 0.01);
    
    bool isRegrid = false;
    for (int t = 0; t < nSteps; ++t) {
        
        IpplTimings::startTimer(solveWithTimer);
        doSolve(myAmrOpal, bunch.get(), rhs_w, phi_w, grad_phi_w, geoms,
                rr, nLevels, bool(t), isRegrid, msg, assignTimer, actualSolveWithTimer, initGuessWithTimer);
        IpplTimings::stopTimer(solveWithTimer);
        
//         std::string plotsmart = BoxLib::Concatenate("smart", t, 4);
//         writePlotFile(plotsmart, rhs_w, phi_w, grad_phi_w, rr, geoms, t);
        
        IpplTimings::startTimer(solveWithoutTimer);
        doSolve(myAmrOpal, bunch.get(), rhs_wo, phi_wo, grad_phi_wo, geoms,
                rr, nLevels, false, isRegrid, msg, assignTimer, actualSolveWithoutTimer, initGuessWithoutTimer);
        IpplTimings::stopTimer(solveWithoutTimer);
        
//         std::string plotnormal = BoxLib::Concatenate("normal", t, 4);
//         writePlotFile(plotnormal, rhs_wo, phi_wo, grad_phi_wo, rr, geoms, t);
        
        /*
         * compare results of potential
         */
        msg << "Timestep " << t << endl;
        for (int lev = 0; lev < nLevels; ++lev) {
            /* we can overwrite the solution "without smart guess"
            * since it will be reset to zero anyway
            */
            MultiFab::Subtract(phi_wo[lev], phi_w[lev], 0, 0, 1, 0);
//             phi_wo[lev].minus(phi_w[lev], 0, 1, 0 /* nghost */);
            double l0error = phi_wo[lev].norm0();
            double l1error = phi_wo[lev].norm1();
            double l2error = phi_wo[lev].norm2();
            msg << "L0-error level " << lev << ": " << l0error << endl;
            msg << "L1-error level " << lev << ": " << l1error << endl;
            msg << "L2-error level " << lev << ": " << l2error << endl;
        }
        
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
            bunch->R[i] += Vector_t(perturbation(eng),
                                    perturbation(eng),
                                    perturbation(eng)
                                   );
        }
        
        isRegrid = false;
        
        IpplTimings::startTimer(regridTimer);
        if ( myAmrOpal.maxLevel() > 0 ) {
            isRegrid = true;
            for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
                myAmrOpal.regrid(i /*lbase*/, t /*time*/);
        } else
            bunch->update();
        IpplTimings::stopTimer(regridTimer);
        
    }
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    
    Inform msg("Solver");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#particles] "
         << "[#levels] [max. box size] [#steps]";
    
    if ( argc < 8 ) {
        msg << call.str() << endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<size_t, 3> nr(std::atoi(argv[1]),
                         std::atoi(argv[2]),
                         std::atoi(argv[3]));
    
    
    size_t nParticles = std::atoi(argv[4]);
    
    
        
    amrex::Initialize(argc,argv, false);
    size_t nLevels = std::atoi(argv[5]) + 1; // i.e. nLevels = 0 --> only single level
    size_t maxBoxSize = std::atoi(argv[6]);
    int nSteps = std::atoi(argv[7]);
    
    
    msg << "Particle test running with" << endl
        << "- #particles      = " << nParticles << endl
        << "- grid            = " << nr << endl
        << "- max grid        = " << maxBoxSize << endl
        << "- #steps          = " << nSteps << endl
        << "- #level          = " << nLevels << endl;
    
    
    doBoxLib(nr, nParticles, nLevels, maxBoxSize, nSteps, msg);
    
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    IpplTimings::print(timefile.str());
    
    amrex::Finalize(true);
    
    return 0;
}
