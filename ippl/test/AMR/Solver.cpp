#include "Solver.h"

#include "Ippl.h"

void 
Solver::solve_for_accel(container_t& rhs,
                        container_t& phi,
                        container_t& grad_phi, 
                        const Array<Geometry>& geom,
                        int base_level,
                        int finest_level,
                        Real offset,
                        bool timing)
{
    static IpplTimings::TimerRef edge2centerTimer;
    
    if ( timing )
        edge2centerTimer = IpplTimings::getTimer("grad-edge2center");
    
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-14;

#ifdef UNIQUE_PTR
    Array< Array<std::unique_ptr<MultiFab> > > grad_phi_edge(rhs.size());
#else
    Array< PArray<MultiFab> > grad_phi_edge(rhs.size());
#endif
    for (int lev = base_level; lev <= finest_level ; lev++)
    {
        grad_phi_edge[lev].resize(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n) {
#ifdef UNIQUE_PTR
            BoxArray ba = rhs[lev]->boxArray();
            grad_phi_edge[lev][n].reset(new MultiFab(ba.surroundingNodes(n), 1, 1));
#else
            BoxArray ba = rhs[lev].boxArray();
            grad_phi_edge[lev].set(n, new MultiFab(ba.surroundingNodes(n), 1, 1));
#endif
        }
    }

    Real     strt    = ParallelDescriptor::second();

    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************
    
    solve_with_f90  (rhs,
                     phi,
                     grad_phi_edge,
                     geom,
                     base_level,
                     finest_level,
                     tol,
                     abs_tol,
                     timing);

    // Average edge-centered gradients to cell centers and fill the values in ghost cells.
    if ( timing )
        IpplTimings::startTimer(edge2centerTimer);
    
    for (int lev = base_level; lev <= finest_level; lev++)
    {
#ifdef UNIQUE_PTR
        BoxLib::average_face_to_cellcenter(*grad_phi[lev],
                                           BoxLib::GetArrOfConstPtrs(grad_phi_edge[lev]),
                                           geom[lev]);
        
        grad_phi[lev]->FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
#else
        BoxLib::average_face_to_cellcenter(grad_phi[lev],
                                           grad_phi_edge[lev],
                                           geom[lev]);
        
        grad_phi[lev].FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
#endif
    }

    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;
    }
    
#ifdef UNIQUE_PTR
    for (int lev = base_level; lev <= finest_level; ++lev) {
        grad_phi[lev]->mult(-1.0, 0, 3);
    }
#else
    for (int lev = base_level; lev <= finest_level; ++lev) {
        grad_phi[lev].mult(-1.0, 0, 3);
    }
#endif
    if ( timing )
        IpplTimings::stopTimer(edge2centerTimer);
}


void 
Solver::solve_with_f90(container_t& rhs,
                       container_t& phi,
                       Array< container_t >& grad_phi_edge,
                       const Array<Geometry>& geom,
                       int base_level,
                       int finest_level,
                       Real tol,
                       Real abs_tol,
                       bool timing)
{
    static IpplTimings::TimerRef initSolverTimer;
    static IpplTimings::TimerRef doSolveTimer;
    static IpplTimings::TimerRef gradientTimer;
    
    if ( timing ) {
        initSolverTimer = IpplTimings::getTimer("init-solver");
        doSolveTimer = IpplTimings::getTimer("do-solve");
        gradientTimer = IpplTimings::getTimer("gradient");
    }
    
    if ( timing )
        IpplTimings::startTimer(initSolverTimer);
    
    int nlevs = finest_level - base_level + 1;

    int mg_bc[2*BL_SPACEDIM];

    // This tells the solver that we are using Dirichlet bc's
    if (Geometry::isAllPeriodic()) {
//         if ( ParallelDescriptor::IOProcessor() )
//             std::cerr << "Periodic BC" << std::endl;
        
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            // periodic BC
            mg_bc[2*dir + 0] = MGT_BC_PER;
            mg_bc[2*dir + 1] = MGT_BC_PER;
        }
    } else if ( Geometry::isAnyPeriodic() ) {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            if ( Geometry::isPeriodic(dir) ) {
                mg_bc[2*dir + 0] = MGT_BC_PER;
                mg_bc[2*dir + 1] = MGT_BC_PER;
            } else {
                mg_bc[2*dir + 0] = MGT_BC_DIR;
                mg_bc[2*dir + 1] = MGT_BC_DIR;
            }
        }
    } else {
//         if ( ParallelDescriptor::IOProcessor() )
//             std::cerr << "Dirichlet BC" << std::endl;
        
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            // Dirichlet BC
            mg_bc[2*dir + 0] = MGT_BC_DIR;
            mg_bc[2*dir + 1] = MGT_BC_DIR;
        }
    }

    // Have to do some packing because these arrays does not always start with base_level
#ifdef UNIQUE_PTR
    Array<Geometry> geom_p(nlevs);
    container_pt rhs_p(nlevs);
    container_pt phi_p(nlevs);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p[ilev] = geom[ilev+base_level];
        rhs_p[ilev] = rhs[ilev+base_level].get();
        phi_p[ilev] = phi[ilev+base_level].get();
    }
    
    // Refinement ratio is hardwired to 2 here.
    IntVect crse_ratio = (base_level == 0) ? 
	IntVect::TheZeroVector() : IntVect::TheUnitVector() * 2;

    FMultiGrid fmg(geom_p, base_level, crse_ratio);

    if (base_level == 0) {
	fmg.set_bc(mg_bc, *phi[base_level]);
    } else {
	fmg.set_bc(mg_bc, *phi[base_level-1], *phi[base_level]);
    }
#else
    PArray<Geometry> geom_p(nlevs);
    container_t rhs_p(nlevs);
    container_t phi_p(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p.set(ilev, &geom[ilev+base_level]);
        rhs_p.set(ilev, &rhs[ilev+base_level]);
        phi_p.set(ilev, &phi[ilev+base_level]);
    }
    
    // Refinement ratio is hardwired to 2 here.
    IntVect crse_ratio = (base_level == 0) ? 
	IntVect::TheZeroVector() : IntVect::TheUnitVector() * 2;

    FMultiGrid fmg(geom_p, base_level, crse_ratio);

    if (base_level == 0) {
	fmg.set_bc(mg_bc, phi[base_level]);
    } else {
	fmg.set_bc(mg_bc, phi[base_level-1], phi[base_level]);
    }
#endif
    
    /* (alpha * a - beta * (del dot b grad)) phi = rhs
     * (b = 1)
     * 
     * The function call set_const_gravity_coeffs() sets alpha = 0.0
     * and beta = -1 (in MGT_Solver::set_const_gravity_coeffs)
     * 
     * --> (del dot grad) phi = rhs
     */
    fmg.set_const_gravity_coeffs();
    fmg.set_maxorder(3);

    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    fmg.set_verbose(1);
    
    if ( timing )
        IpplTimings::stopTimer(initSolverTimer);
    
    if ( timing )
        IpplTimings::startTimer(doSolveTimer);
    fmg.solve(phi_p, rhs_p, tol, abs_tol, always_use_bnorm, need_grad_phi);
    if ( timing )
        IpplTimings::stopTimer(doSolveTimer);
    
    if ( timing )
        IpplTimings::startTimer(gradientTimer);
#ifdef UNIQUE_PTR
    const Array<Array<MultiFab*> >& g_phi_edge = BoxLib::GetArrOfArrOfPtrs(grad_phi_edge);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + base_level;
        fmg.get_fluxes(g_phi_edge[amr_level], ilev);
    }
#else
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + base_level;
        fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
    }
#endif
    if ( timing )
        IpplTimings::stopTimer(gradientTimer);
}

#ifdef USEHYPRE
// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays
void Solver::solve_with_hypre(MultiFab& soln, MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
    int  verbose       = 2;
    Real tolerance_rel = 1.e-8;
    Real tolerance_abs = 0.0;
    int  maxiter       = 100;
    BL_PROFILE("solve_with_hypre()");
    BndryData bd(bs, 1, geom);
    set_boundary(bd, rhs, 0);
    
    Real a = 0.0;
    Real b = 1.0;
    
    // Set up the Helmholtz operator coefficients.
    MultiFab alpha(bs, 1, 0);
    alpha.setVal(0.0);
    
    PArray<MultiFab> beta(BL_SPACEDIM, PArrayManage);
    for ( int n=0; n<BL_SPACEDIM; ++n ) {
        BoxArray bx(bs);
        beta.set(n, new MultiFab(bx.surroundingNodes(n), 1, 0, Fab_allocate));
        beta[n].setVal(1.0);
    }
    
    HypreABecLap hypreSolver(bs, geom);
    hypreSolver.setScalars(a, b);
    hypreSolver.setACoeffs(alpha);
    hypreSolver.setBCoeffs(beta);
    hypreSolver.setVerbose(verbose);
    hypreSolver.solve(soln, rhs, tolerance_rel, tolerance_abs, maxiter, bd);
}


void Solver::set_boundary(BndryData& bd, const MultiFab& rhs, int comp)
{
  BL_PROFILE("set_boundary()");
  Real bc_value = 0.0;

  for (int n=0; n<BL_SPACEDIM; ++n) {
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
      int i = mfi.index(); 
      
      const Box& bx = mfi.validbox();
      
      // Our default will be that the face of this grid is either touching another grid
      //  across an interior boundary or a periodic boundary.  We will test for the other
      //  cases below.
      {
	// Define the type of boundary conditions to be Dirichlet (even for periodic)
	bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	
	// Set the boundary conditions to the cell centers outside the domain
	bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
	bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
      }

      // Now test to see if we should override the above with Dirichlet or Neumann physical bc's
//       if (bc_type != Periodic) { 
	int ibnd = static_cast<int>(LO_DIRICHLET);
	const Geometry& geom = bd.getGeom();

	// We are on the low side of the domain in coordinate direction n
	if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
	  // Set the boundary conditions to live exactly on the faces of the domain
	  bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  
	  // Set the Dirichlet/Neumann boundary values 
	  bd.setValue(Orientation(n, Orientation::low) ,i, bc_value);
	  
	  // Define the type of boundary conditions 
	  bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,ibnd);
	}
	
	// We are on the high side of the domain in coordinate direction n
	if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
	  // Set the boundary conditions to live exactly on the faces of the domain
	  bd.setBoundLoc(Orientation(n, Orientation::high) ,i,0.0 );
	  
	  // Set the Dirichlet/Neumann boundary values
	  bd.setValue(Orientation(n, Orientation::high) ,i, bc_value);

	  // Define the type of boundary conditions 
	  bd.setBoundCond(Orientation(n, Orientation::high) ,i,comp,ibnd);
	}
//       }
    }
  }
}
#endif