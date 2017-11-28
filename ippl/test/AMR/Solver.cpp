#include "Solver.h"

void 
Solver::solve_for_accel(const container_t& rhs,
                        const container_t& phi,
                        const container_t& grad_phi, 
                        const amrex::Array<amrex::Geometry>& geom,
                        int base_level,
                        int finest_level,
                        amrex::Real offset,
                        bool timing,
                        bool doGradient)
{
    using amrex::Array;

    static IpplTimings::TimerRef edge2centerTimer;
    
    if ( timing )
        edge2centerTimer = IpplTimings::getTimer("grad-edge2center");
    
//     amrex::Real reltol = 1.0e-14;
//     amrex::Real abstol = 1.0e-12;
    
//     amrex::Real reltol = 1.0e-12;
    amrex::Real reltol = 1.0e-14;
    amrex::Real abstol = 1.0e-10;

    amrex::Array<container_t> grad_phi_edge(rhs.size());
    
    if ( doGradient ) {
        for (int lev = base_level; lev <= finest_level ; lev++)
        {
            const amrex::DistributionMapping& dm = rhs[lev]->DistributionMap();
            grad_phi_edge[lev].resize(AMREX_SPACEDIM);
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                amrex::BoxArray ba = rhs[lev]->boxArray();
                grad_phi_edge[lev][n].reset(new amrex::MultiFab(ba.surroundingNodes(n), dm, 1, 1));
            }
        }
    }
    
    // ***************************************************
    // Make sure the RHS sums to 0 if fully periodic
    // ***************************************************
    for (int lev = base_level; lev <= finest_level; lev++)
        rhs[lev]->plus(+offset, 0, 1, 0);
    
    
    
    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************
    
    solve_with_f90  (amrex::GetArrOfPtrs(rhs),
                     amrex::GetArrOfPtrs(phi),
                     amrex::GetArrOfArrOfPtrs(grad_phi_edge),
                     geom,
                     base_level,
                     finest_level,
                     reltol,
                     abstol,
                     timing,
                     doGradient);

    // Average edge-centered gradients to cell centers and fill the values in ghost cells.
    if ( timing )
        IpplTimings::startTimer(edge2centerTimer);
    
    if ( doGradient ) {
        for (int lev = base_level; lev <= finest_level; lev++)
        {
            amrex::average_face_to_cellcenter(*(grad_phi[lev].get()),
                                              amrex::GetArrOfConstPtrs(grad_phi_edge[lev]),
                                              geom[lev]);
        
            grad_phi[lev]->FillBoundary(0,AMREX_SPACEDIM,geom[lev].periodicity());
        }
        
        for (int lev = base_level; lev <= finest_level; ++lev) {
            grad_phi[lev]->mult(-1.0, 0, 3);
        }
    }

    if ( timing )
        IpplTimings::stopTimer(edge2centerTimer);
}


void 
Solver::solve_with_f90(const container_pt& rhs,
                       const container_pt& phi,
                       const amrex::Array<container_pt>& grad_phi_edge,
                       const amrex::Array<amrex::Geometry>& geom,
                       int base_level,
                       int finest_level,
                       amrex::Real reltol,
                       amrex::Real abstol,
                       bool timing,
                       bool doGradient)
{
    static IpplTimings::TimerRef initSolverTimer;
    static IpplTimings::TimerRef doSolveTimer;
    static IpplTimings::TimerRef gradientTimer;
    
    if ( timing ) {
        initSolverTimer = IpplTimings::getTimer("init-solver");
        doSolveTimer = IpplTimings::getTimer("do-solve");
        if ( doGradient )
            gradientTimer = IpplTimings::getTimer("gradient");
    }
    
    if ( timing )
        IpplTimings::startTimer(initSolverTimer);
    
    int nlevs = finest_level - base_level + 1;

    int mg_bc[2*AMREX_SPACEDIM];

    // This tells the solver that we are using Dirichlet bc's
    if (amrex::Geometry::isAllPeriodic()) {
//         if ( ParallelDescriptor::IOProcessor() )
//             std::cerr << "Periodic BC" << std::endl;
        
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // periodic BC
            mg_bc[2*dir + 0] = MGT_BC_PER;
            mg_bc[2*dir + 1] = MGT_BC_PER;
        }
    } else if ( amrex::Geometry::isAnyPeriodic() ) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if ( amrex::Geometry::isPeriodic(dir) ) {
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
        
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Dirichlet BC
            mg_bc[2*dir + 0] = MGT_BC_DIR;
            mg_bc[2*dir + 1] = MGT_BC_DIR;
        }
    }

    // Have to do some packing because these arrays does not always start with base_level
    amrex::Array<amrex::Geometry> geom_p(nlevs);
    container_pt rhs_p(nlevs);
    container_pt phi_p(nlevs);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p[ilev] = geom[ilev+base_level];
        rhs_p[ilev]  = rhs[ilev+base_level];
        phi_p[ilev]  = phi[ilev+base_level];
    }
    
    //FIXME Refinement ratio is hardwired to 2 here.
    amrex::IntVect crse_ratio = (base_level == 0) ? 
	amrex::IntVect::TheZeroVector() : amrex::IntVect::TheUnitVector() * 2;

    amrex::FMultiGrid fmg(geom_p, base_level, crse_ratio);

    if (base_level == 0) {
	fmg.set_bc(mg_bc, *phi_p[base_level]);
    } else {
	fmg.set_bc(mg_bc, *phi_p[base_level-1], *phi_p[base_level]);
    }
    
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
    int need_grad_phi = (doGradient) ? 1 : 0;
    fmg.set_verbose(0);
    
    if ( timing )
        IpplTimings::stopTimer(initSolverTimer);
    
    if ( timing )
        IpplTimings::startTimer(doSolveTimer);
    amrex::Real final_resnorm = fmg.solve(phi_p, rhs_p, reltol, abstol, always_use_bnorm, need_grad_phi);
    
    if ( final_resnorm > abstol ) {
        std::stringstream ss;
        ss << "Residual norm: " << std::setprecision(16) << final_resnorm
           << " > " << abstol << " (absolute tolerance)";
        throw std::runtime_error("\033[1;31mError: The solver did not converge: " +
                                 ss.str() + "\033[0m");
    }
    
    if ( timing )
        IpplTimings::stopTimer(doSolveTimer);
    
    if ( doGradient ) {
        
        if ( timing )
            IpplTimings::startTimer(gradientTimer);
        for (int ilev = 0; ilev < nlevs; ++ilev) {
            int amr_level = ilev + base_level;
            fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
        }
        if ( timing )
            IpplTimings::stopTimer(gradientTimer);
    }
}

#ifdef USEHYPRE
// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays
void Solver::solve_with_hypre(amrex::MultiFab& soln, amrex::MultiFab& rhs, const amrex::BoxArray& bs, const amrex::Geometry& geom)
{
    int  verbose       = 0;
    amrex::Real tolerance_rel = 1.e-8;
    amrex::Real tolerance_abs = 0.0;
    int  maxiter       = 100;
    BL_PROFILE("solve_with_hypre()");
    BndryData bd(bs, 1, geom);
    set_boundary(bd, rhs, 0);
    
    amrex::Real a = 0.0;
    amrex::Real b = 1.0;
    
    // Set up the Helmholtz operator coefficients.
    amrex::MultiFab alpha(bs, 1, 0);
    alpha.setVal(0.0);
    
    PArray<amrex::MultiFab> beta(AMREX_SPACEDIM, PArrayManage);
    for ( int n=0; n<AMREX_SPACEDIM; ++n ) {
        amrex::BoxArray bx(bs);
        beta.set(n, new amrex::MultiFab(bx.surroundingNodes(n), 1, 0, Fab_allocate));
        beta[n].setVal(1.0);
    }
    
    HypreABecLap hypreSolver(bs, geom);
    hypreSolver.setScalars(a, b);
    hypreSolver.setACoeffs(alpha);
    hypreSolver.setBCoeffs(beta);
    hypreSolver.setVerbose(verbose);
    hypreSolver.solve(soln, rhs, tolerance_rel, tolerance_abs, maxiter, bd);
}


void Solver::set_boundary(BndryData& bd, const amrex::MultiFab& rhs, int comp)
{
  BL_PROFILE("set_boundary()");
  amrex::Real bc_value = 0.0;

  for (int n=0; n<AMREX_SPACEDIM; ++n) {
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
	const amrex::Geometry& geom = bd.getGeom();

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
