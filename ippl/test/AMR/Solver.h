#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <VisMF.H>
#include <FMultiGrid.H>


class Solver {

public:
    void solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi,
			 const Array<Geometry>& geom, int base_level, int finest_level, Real offset);

    void solve_with_f90(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, Array< PArray<MultiFab> >& grad_phi_edge, 
			const Array<Geometry>& geom, int base_level, int finest_level, Real tol, Real abs_tol);
};


void 
Solver::solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
		const Array<Geometry>& geom, int base_level, int finest_level, Real offset)
{
 
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-14;

    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(rhs.size());

    for (int lev = base_level; lev <= finest_level ; lev++)
    {
        grad_phi_edge[lev].resize(BL_SPACEDIM, PArrayManage);
        for (int n = 0; n < BL_SPACEDIM; ++n)
            grad_phi_edge[lev].set(n, new MultiFab(BoxArray(rhs[lev].boxArray()).surroundingNodes(n), 1, 1));
    }

    Real     strt    = ParallelDescriptor::second();

    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************
    
    solve_with_f90  (rhs,phi,grad_phi_edge,geom,base_level,finest_level,tol,abs_tol);

    // Average edge-centered gradients to cell centers and fill the values in ghost cells.
    for (int lev = base_level; lev <= finest_level; lev++)
    {
        BoxLib::average_face_to_cellcenter(grad_phi[lev], grad_phi_edge[lev], geom[lev]);
        grad_phi[lev].FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
    }

    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;
    }
}


void 
Solver::solve_with_f90(PArray<MultiFab>& rhs, PArray<MultiFab>& phi,
		       Array< PArray<MultiFab> >& grad_phi_edge, 
		       const Array<Geometry>& geom, int base_level,
		       int finest_level, Real tol, Real abs_tol)
{
    int nlevs = finest_level - base_level + 1;

    int mg_bc[2*BL_SPACEDIM];

    // This tells the solver that we are using Dirichlet bc's
    if (Geometry::isAllPeriodic()) {
        if ( ParallelDescriptor::IOProcessor() )
            std::cerr << "Periodic BC" << std::endl;
        
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            // periodic BC
            mg_bc[2*dir + 0] = MGT_BC_PER;
            mg_bc[2*dir + 1] = MGT_BC_PER;
        }
    } else {
        if ( ParallelDescriptor::IOProcessor() )
            std::cerr << "Dirichlet BC" << std::endl;
        
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            // Dirichlet BC
            mg_bc[2*dir + 0] = MGT_BC_DIR;
            mg_bc[2*dir + 1] = MGT_BC_DIR;
        }
    }

    // Have to do some packing because these arrays does not always start with base_level
    PArray<Geometry> geom_p(nlevs);
    PArray<MultiFab> rhs_p(nlevs);
    PArray<MultiFab> phi_p(nlevs);
    for (int ilev = 0; ilev < nlevs; ++ilev) {
	geom_p.set(ilev, &geom[ilev+base_level]);
	rhs_p.set (ilev,  &rhs[ilev+base_level]);
	phi_p.set (ilev,  &phi[ilev+base_level]);
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

    fmg.set_const_gravity_coeffs();

    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    fmg.set_verbose(0);
    fmg.solve(phi_p, rhs_p, tol, abs_tol, always_use_bnorm, need_grad_phi);
   
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + base_level;
        fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
    }
}

#endif
