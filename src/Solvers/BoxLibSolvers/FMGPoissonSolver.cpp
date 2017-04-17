#include "FMGPoissonSolver.h"

FMGPoissonSolver::FMGPoissonSolver(AmrBoxLib* itsAmrObject_p)
    : AmrPoissonSolver<AmrBoxLib>(itsAmrObject_p),
      tol_m(1.e-10),
      abstol_m(1.e-14)//,
//       phi_m(PArrayManage)
{
    // Dirichlet boundary conditions are default
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        bc_m[2 * d]     = MGT_BC_DIR;
        bc_m[2 * d + 1] = MGT_BC_DIR;
    }
}

void FMGPoissonSolver::solve(AmrFieldContainer_t &rho,
                             AmrFieldContainer_t& phi,
                             AmrFieldContainer_t &efield,
                             unsigned short baseLevel,
                             unsigned short finestLevel)
{
    const GeomContainer_t& geom = itsAmrObject_mp->Geom();
    
    
    if (Geometry::isAllPeriodic()) {
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            bc_m[2 * d]     = MGT_BC_PER;
            bc_m[2 * d + 1] = MGT_BC_PER;
        }
    } else if ( Geometry::isAnyPeriodic() ) {
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            if ( Geometry::isPeriodic(d) ) {
                bc_m[2 * d]     = MGT_BC_PER;
                bc_m[2 * d + 1] = MGT_BC_PER;
            }
        }
    }
    
    
    Array<AmrFieldContainer_t> grad_phi_edge(rho.size());
    
    
    for (int lev = baseLevel; lev <= finestLevel ; ++lev) {
        grad_phi_edge[lev].resize(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n) {
            BoxArray ba = rho[lev].boxArray();
            grad_phi_edge[lev].set(n, new MultiFab(ba.surroundingNodes(n), 1, 1));
        }
    }
    
    // initialize potential and electric field grids on each level
    initGrids_m(phi, efield);
    
    this->solveWithF90_m(rho, phi, grad_phi_edge, geom,
                         baseLevel, finestLevel);
    
    
    
    for (int lev = baseLevel; lev <= finestLevel; ++lev) {
        BoxLib::average_face_to_cellcenter(efield[lev],
                                           grad_phi_edge[lev],
                                           geom[lev]);
        
        efield[lev].FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
        efield[lev].mult(-1.0, 0, 3);
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
       << "* tolerance " << tol_m << '\n'
       << "* absolute tolerance " << abstol_m << '\n' << endl
       << "* ******************************************************************** " << endl;
    return os;
}


void FMGPoissonSolver::solveWithF90_m(AmrFieldContainer_t& rho,
                                      AmrFieldContainer_t& phi,
                                      Array<AmrFieldContainer_t>& grad_phi_edge,
                                      const GeomContainer_t& geom,
                                      int baseLevel,
                                      int finestLevel)
{
    int nlevs = finestLevel - baseLevel + 1;
    
    PArray<Geometry> geom_p(nlevs, PArrayManage);
    AmrFieldContainer_t rho_p(nlevs, PArrayManage);
    AmrFieldContainer_t phi_p(nlevs, PArrayManage);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        geom_p.set(ilev, &geom[ilev + baseLevel]);
        rho_p.set(ilev, &rho[ilev + baseLevel]);
        phi_p.set(ilev, &phi[ilev + baseLevel]);
    }
    
    // Refinement ratio is hardwired to 2 here.
    IntVect crse_ratio = (baseLevel == 0) ?
        IntVect::TheZeroVector() : IntVect::TheUnitVector() * 2;

    FMultiGrid fmg(geom_p, baseLevel, crse_ratio);

    if (baseLevel == 0)
        fmg.set_bc(bc_m, phi[baseLevel]);
    else
        fmg.set_bc(bc_m, phi[baseLevel-1], phi[baseLevel]);
    
    /* (alpha * a - beta * (del dot b grad)) phi = rho
     * (b = 1)
     * 
     * The function call set_const_gravity_coeffs() sets alpha = 0.0
     * and beta = -1 (in MGT_Solver::set_const_gravity_coeffs)
     * 
     * --> (del dot grad) phi = rho
     */
    fmg.set_const_gravity_coeffs();
    fmg.set_maxorder(3);

    int always_use_bnorm = 0;
    int need_grad_phi = 1;
    
    fmg.solve(phi_p, rho_p, tol_m, abstol_m, always_use_bnorm, need_grad_phi);
    
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        int amr_level = ilev + baseLevel;
        fmg.get_fluxes(grad_phi_edge[amr_level], ilev);
    }
}

void FMGPoissonSolver::initGrids_m(AmrFieldContainer_t& phi,
                                   AmrFieldContainer_t& efield) {
    for (int lev = 0; lev < efield.size(); ++lev) {
        phi.clear(lev);
        efield.clear(lev);
        
        const BoxArray& ba = itsAmrObject_mp->boxArray()[lev];
        
        //                                      #components  #ghost cells                                                                                                                                          
        phi.set(lev,  new AmrField_t(ba, 1,           1) );
        efield.set(lev, new AmrField_t(ba, 3,           1) );
        
        phi[lev].setVal(0.0);
        efield[lev].setVal(0.0);
    }
}
