#include "FMGPoissonSolver.h"

#include "Utilities/OpalException.h"

FMGPoissonSolver::FMGPoissonSolver(AmrBoxLib* itsAmrObject_p)
    : AmrPoissonSolver<AmrBoxLib>(itsAmrObject_p),
      reltol_m(1.e-14),
      abstol_m(0.0)
{
    // Dirichlet boundary conditions are default
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        bc_m[2 * d]     = MGT_BC_DIR;
        bc_m[2 * d + 1] = MGT_BC_DIR;
    }
}

void FMGPoissonSolver::solve(AmrFieldContainer_t& rho,
                             AmrFieldContainer_t& phi,
                             AmrFieldContainer_t& efield,
                             unsigned short baseLevel,
                             unsigned short finestLevel)
{
    const GeomContainer_t& geom = itsAmrObject_mp->Geom();
    
    
    if (AmrGeometry_t::isAllPeriodic()) {
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            bc_m[2 * d]     = MGT_BC_PER;
            bc_m[2 * d + 1] = MGT_BC_PER;
        }
    } else if ( AmrGeometry_t::isAnyPeriodic() ) {
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            if ( AmrGeometry_t::isPeriodic(d) ) {
                bc_m[2 * d]     = MGT_BC_PER;
                bc_m[2 * d + 1] = MGT_BC_PER;
            }
        }
    }
    
    amrex::Array< AmrFieldContainer_t > grad_phi_edge(rho.size());
    
    for (int lev = baseLevel; lev <= finestLevel ; ++lev) {
        const AmrProcMap_t& dmap = rho[lev]->DistributionMap();
        grad_phi_edge[lev].resize(BL_SPACEDIM);
        AmrGrid_t ba = rho[lev]->boxArray();
        
        for (int n = 0; n < BL_SPACEDIM; ++n) {
            grad_phi_edge[lev][n].reset(new AmrField_t(ba.surroundingNodes(n), dmap, 1, 1));
        }
    }
    
    // initialize potential and electric field grids on each level
    phi.resize( rho.size() );
    efield.resize( rho.size() );
    initGrids_m(rho, phi, efield, baseLevel, finestLevel);
    
    double residNorm = this->solveWithF90_m(amrex::GetArrOfPtrs(rho),
                                            amrex::GetArrOfPtrs(phi),
                                            amrex::GetArrOfArrOfPtrs(grad_phi_edge),
                                            geom,
                                            baseLevel,
                                            finestLevel);
    
    if ( residNorm > reltol_m )
        throw OpalException("FMGPoissonSolver::solve()",
                            "Multigrid solver did not converge. Residual norm: "
                            + std::to_string(residNorm));
    
    
    for (int lev = baseLevel; lev <= finestLevel; ++lev) {
        amrex::average_face_to_cellcenter(*(efield[lev].get()),
                                          amrex::GetArrOfConstPtrs(grad_phi_edge[lev]),
                                          geom[lev]);
        
        efield[lev]->FillBoundary(0, BL_SPACEDIM,geom[lev].periodicity());
        efield[lev]->mult(-1.0, 0, 3);
    }
}


double FMGPoissonSolver::getXRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(0);
}


double FMGPoissonSolver::getXRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(0);
}


double FMGPoissonSolver::getYRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(1);
}


double FMGPoissonSolver::getYRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(1);
}


double FMGPoissonSolver::getZRangeMin(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbLo(2);
}


double FMGPoissonSolver::getZRangeMax(unsigned short level) {
    return itsAmrObject_mp->Geom(level).ProbHi(2);
}


Inform &FMGPoissonSolver::print(Inform &os) const {
    os << "* ************* F M G P o i s s o n S o l v e r ************************************ " << endl
       << "* relative tolerance " << reltol_m << '\n'
       << "* absolute tolerance " << abstol_m << '\n' << endl
       << "* ******************************************************************** " << endl;
    return os;
}


double FMGPoissonSolver::solveWithF90_m(const AmrFieldContainer_pt& rho,
                                        const AmrFieldContainer_pt& phi,
                                        const amrex::Array< AmrFieldContainer_pt > & grad_phi_edge,
                                        const GeomContainer_t& geom,
                                        int baseLevel,
                                        int finestLevel)
{
    int nlevs = finestLevel - baseLevel + 1;
    
    GeomContainer_t geom_p(nlevs);
    AmrFieldContainer_pt rho_p(nlevs);
    AmrFieldContainer_pt phi_p(nlevs);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p[ilev] = geom[ilev + baseLevel];
        rho_p[ilev]  = rho[ilev + baseLevel];
        phi_p[ilev]  = phi[ilev + baseLevel];
    }
    
    //FIXME Refinement ratio
    amrex::IntVect crse_ratio = (baseLevel == 0) ?
        amrex::IntVect::TheZeroVector() : itsAmrObject_mp->refRatio(0);

    amrex::FMultiGrid fmg(geom_p, baseLevel, crse_ratio);

    if (baseLevel == 0)
        fmg.set_bc(bc_m, *phi[baseLevel]);
    else
        fmg.set_bc(bc_m, *phi[baseLevel-1], *phi[baseLevel]);
    
    /* (alpha * a - beta * (del dot b grad)) phi = rho
     * (b = 1)
     * 
     * The function call set_const_gravity_coeffs() sets alpha = 0.0
     * and beta = -1 (in MGT_Solver::set_const_gravity_coeffs)
     * 
     * --> (del dot grad) phi = rho
     */
    fmg.set_const_gravity_coeffs();
    
    // order of approximation at Dirichlet boundaries
    fmg.set_maxorder(3);

    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    
    double residNorm = fmg.solve(phi_p, rho_p, reltol_m, abstol_m, always_use_bnorm, need_grad_phi);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + baseLevel;
        fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
    }
    
    return residNorm;
}

void FMGPoissonSolver::initGrids_m(const AmrFieldContainer_t& rho,
                                   AmrFieldContainer_t& phi,
                                   AmrFieldContainer_t& efield,
                                   int baseLevel,
                                   int finestLevel)
{
    int nlevs = finestLevel - baseLevel + 1;
    for (int lev = 0; lev < nlevs; ++lev) {
        int ilev = lev + baseLevel;
        const AmrGrid_t& ba = rho[ilev]->boxArray();
        const AmrProcMap_t& dmap = rho[ilev]->DistributionMap();
        
        //                                         #components  #ghost cells                                                                                                                                          
        phi[ilev].reset(   new AmrField_t(ba, dmap, 1,           1) );
        efield[ilev].reset(new AmrField_t(ba, dmap, 3,           1) );
        
        // including nghost = 1
        phi[lev]->setVal(0.0, 1);
        efield[lev]->setVal(0.0, 1);
    }
}
