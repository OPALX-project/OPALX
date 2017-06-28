#include "MGTSolver.h"

#include <AMReX_MacBndry.H>
#include <AMReX_stencil_types.H>
#include <AMReX_MultiFabUtil.H>

#include <iomanip>

void MGTSolver::solve(const container_t& rho,
                      container_t& phi,
                      container_t& efield,
                      const amrex::Array<amrex::Geometry>& geom)
{
//     Real reltol = 1.0e-14;
//     Real abstol = 1.0e-12;
    double abstol = 0.0;
    double reltol = 1.0e-12;
    
    int baseLevel = 0;
    int finestLevel = rho.size() - 1;
    
    const int num_levels = finestLevel - baseLevel + 1;
    
    
    amrex::Array<amrex::BoxArray> bav(num_levels);
    amrex::Array<amrex::DistributionMapping> dmv(num_levels);
    
    for (int lev = 0; lev < num_levels; lev++)
    {
        bav[lev]          = rho[baseLevel+lev]->boxArray();
        dmv[lev]          = rho[baseLevel+lev]->DistributionMap();
        
        phi[lev]          = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(bav[lev], dmv[lev], 1          , 1));
        efield[lev]       = std::unique_ptr<amrex::MultiFab>(new amrex::MultiFab(bav[lev], dmv[lev], BL_SPACEDIM, 1));
        efield[lev]->setVal(0.0, 1);
    }
    
    amrex::Array< amrex::Array<std::unique_ptr<amrex::MultiFab> > > grad_phi_prev;
    grad_phi_prev.resize(num_levels);
    
    for (int lev = baseLevel; lev <= finestLevel; lev++)
    {
        grad_phi_prev[lev].resize(BL_SPACEDIM);
        BL_ASSERT(grad_phi_prev[lev].size() == BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = amrex::BoxArray(bav[lev]).surroundingNodes(n);
            grad_phi_prev[lev][n].reset(new amrex::MultiFab(eba, dmv[lev], 1, 1));
        }
    }
    
    amrex::Array<amrex::Geometry> fgeom(num_levels);
    for (int i = 0; i < num_levels; i++)
        fgeom[i] = geom[baseLevel+i];
    
    
    Array< Array<Real> > xa(num_levels);
    Array< Array<Real> > xb(num_levels);

    for (int lev = 0; lev < num_levels; lev++)
    {
        xa[lev].resize(BL_SPACEDIM);
        xb[lev].resize(BL_SPACEDIM);
        if (baseLevel + lev == 0)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0;
                xb[lev][i] = 0;
            }
        }
        else
        {
            const Real* dx_crse = geom[baseLevel + lev - 1].CellSize();
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0.5 * dx_crse[i];
                xb[lev][i] = 0.5 * dx_crse[i];
            }
        }
    }
    
    
    amrex::Array<amrex::MultiFab*> phi_p(num_levels);
    amrex::Array<amrex::MultiFab*> rho_p(num_levels);
//     amrex::Array<std::unique_ptr<amrex::MultiFab> > Rhs_p(num_levels);
    
    for (int lev = 0; lev < num_levels; lev++)
    {
        phi_p[lev] = (phi[lev].get());
        
        phi_p[lev]->setVal(0);
        
        rho_p[lev] = rho[lev].get();

//         // Need to set the boundary values before "bndry" is defined so they get copied in
//         if (dirichlet_bcs) set_dirichlet_bcs(baseLevel+lev,phi_p[lev]);

//         Rhs_p[lev].reset(new amrex::MultiFab(bav[baseLevel+lev], dmv[baseLevel+lev], 1, 0));
//         Rhs_p[lev]->setVal(0.0);
    }
    
//     Array<MultiFab*> rho_p = { &rho };
    
    amrex::IntVect crse_ratio = (baseLevel == 0) ? 
        amrex::IntVect::TheZeroVector() : amrex::IntVect::TheUnitVector() * 2;
    
    //
    // Store the Dirichlet boundary condition for phi in bndry.
    //
    MacBndry bndry(bav[baseLevel], dmv[baseLevel], 1, geom[baseLevel]);
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;
    
    amrex::BCRec phys_bc;
    
    // Get boundary conditions
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i, 0);
        phys_bc.setHi(i, 0);
    }
    
    
    bndry.setBndryValues(*phi_p[0], src_comp, dest_comp, num_comp, phys_bc);
    
    int stencil_type = amrex::CC_CROSS_STENCIL;
    int mg_bc[2*BL_SPACEDIM];
    
    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        if (geom[0].isPeriodic(dir))
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
        else
        {
            mg_bc[2*dir + 0] = MGT_BC_DIR;
            mg_bc[2*dir + 1] = MGT_BC_DIR;
        }
    }
    
    MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type, false, 0, 1, 0);
    
    mgt_solver.set_const_gravity_coeffs(xa, xb);
    
    Real final_resnorm = 1;
    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    mgt_solver.set_maxorder(3);
//     amrex::Array<amrex::MultiFab*> phi_p = { &phi };
//     amrex::Array<amrex::MultiFab*> rho_p = { &rho };
    
    //
    // Call the solver
    //
    mgt_solver.solve(phi_p, rho_p,
                     bndry, reltol, abstol, always_use_bnorm, final_resnorm, need_grad_phi);
    
    if ( final_resnorm > reltol ) {
        std::stringstream ss;
        ss << "Residual norm: " << std::setprecision(16) << final_resnorm
           << " > " << reltol << " (relative tolerance)";
        throw std::runtime_error("\033[1;31mError: The solver did not converge: " +
                                 ss.str() + "\033[0m");
    }

    
        for (int lev = 0; lev < num_levels; lev++)
        {
            const Real* dx = geom[baseLevel+lev].CellSize();
            mgt_solver.get_fluxes(lev, amrex::GetArrOfPtrs(grad_phi_prev[baseLevel+lev]), dx);
        }
        
       IntVect fine_ratio = amrex::IntVect::TheUnitVector() * 2; 
        
        // Average phi from fine to coarse level
    for (int lev = finestLevel; lev > baseLevel; lev--)
    {
        amrex::average_down(*phi_p[lev],
                            *phi_p[lev-1],
                            0, 1, fine_ratio);
    }

    // Average grad_phi from fine to coarse level
    for (int lev = finestLevel - 1; lev > baseLevel; lev--) {
        
        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_gphi_fine_BA(bav[lev+1].size());
        
        IntVect fine_ratio = amrex::IntVect::TheUnitVector() * 2;
        
        for (int i = 0; i < crse_gphi_fine_BA.size(); ++i)
            crse_gphi_fine_BA.set(i, amrex::coarsen(bav[lev+1][i],
                                                    fine_ratio));
            
        Array<std::unique_ptr<MultiFab> > crse_gphi_fine(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            const BoxArray eba = BoxArray(crse_gphi_fine_BA).surroundingNodes(n);
            crse_gphi_fine[n].reset(new MultiFab(eba, dmv[lev+1], 1, 0));
        }
        
        amrex::average_down_faces(amrex::GetArrOfConstPtrs(grad_phi_prev[lev+1]),
                                  amrex::GetArrOfPtrs(crse_gphi_fine), fine_ratio);
        
        const Geometry& cgeom = geom[lev];
        
        for (int n = 0; n < BL_SPACEDIM; ++n) {
            grad_phi_prev[lev][n]->copy(*crse_gphi_fine[n], cgeom.periodicity());
        }
        
        
    }
        
        for (int lev = baseLevel; lev <= finestLevel; lev++)
        {
            amrex::average_face_to_cellcenter(*(efield[lev].get()),
                                              amrex::GetArrOfConstPtrs(grad_phi_prev[lev]),
                                              geom[lev]);
        
            efield[lev]->FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
        }
        
        for (int lev = baseLevel; lev <= finestLevel; ++lev) {
            efield[lev]->mult(-1.0, 0, 3);
        }
//         average_fine_ec_onto_crse_ec(lev-1,is_new);
    
}
