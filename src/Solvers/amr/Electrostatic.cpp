#include <cmath>

#include <ParmParse.H>
#include <MacBndry.H>

#include "Electrostatic.H"
#include "Accel.H"
#include "Electrostatic_F.H"
#include "Accel_F.H"
#include "TrilinosSolver.h"


extern Inform *gmsg;

// MAX_LEV defines the maximum number of AMR levels allowed by the parent "Amr" object
#define MAX_LEV 15

// Give this a bogus default value to force user to define in inputs file
int  Electrostatic::verbose       = 2;
int  Electrostatic::show_timings  = 0;
int  Electrostatic::no_sync       = 0;
int  Electrostatic::no_composite  = 0;

// ************************************************************************** //

Electrostatic::Electrostatic (Amr*   Parent,
                  int    _finest_level,
                  BCRec* _phys_bc)
  :
    parent(Parent),
    level_data(MAX_LEV),
    grad_phi_curr(MAX_LEV),
    grad_phi_prev(MAX_LEV),
    phi_flux_reg(MAX_LEV,PArrayManage),
    grids(MAX_LEV),
    level_solver_resnorm(MAX_LEV),
    phys_bc(_phys_bc)
{
     read_params();
     finest_level_allocated = -1;
}


Electrostatic::~Electrostatic ()
{
    // nothing to see here.
}

void
Electrostatic::read_params ()
{
    static bool done = false;

    if (!done)
    {
        const Real strt = ParallelDescriptor::second();

        ParmParse pp("gravity");

        pp.query("v", verbose);
        pp.query("show_timings", show_timings);
        pp.query("no_sync", no_sync);
        pp.query("no_composite", no_composite);

        if (show_timings)
        {
            const int IOProc = ParallelDescriptor::IOProcessorNumber();
            Real end = ParallelDescriptor::second() - strt;

            ParallelDescriptor::ReduceRealMax(end, IOProc);

            if (ParallelDescriptor::IOProcessor())
                *gmsg << "Electrostatic::read_params() time = " << end
                          << '\n';
        }

        done = true;
    }
}

void
Electrostatic::install_level (int       level,
                        AmrLevel* level_data_to_install)
{
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        *gmsg << "Installing Electrostatic level " << level << '\n';

    level_data.clear(level);
    level_data.set(level, level_data_to_install);

    grids[level] = level_data[level].boxArray();

    level_solver_resnorm[level] = 0;

    grad_phi_prev[level].clear();
    grad_phi_prev[level].resize(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        grad_phi_prev[level].set(n, new MultiFab(BoxArray(grids[level]).surroundingNodes(n), 1, 1));
        grad_phi_prev[level][n].setVal(1.e50,grad_phi_prev[level][n].nGrow());
    }

    grad_phi_curr[level].clear();
    grad_phi_curr[level].resize(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        grad_phi_curr[level].set(n, new MultiFab(BoxArray(grids[level]).surroundingNodes(n), 1, 1));
        grad_phi_curr[level][n].setVal(1.e50,grad_phi_curr[level][n].nGrow());
    }

    if (level > 0)
    {
        phi_flux_reg.clear(level);
        IntVect crse_ratio = parent->refRatio(level-1);
        phi_flux_reg.set(level, new FluxRegister(grids[level], crse_ratio,
                                                 level, 1));
    }

    finest_level_allocated = level;
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        *gmsg << "Done installing Electrostatic level " << level << '\n';
}

int
Electrostatic::get_no_sync ()
{
    return no_sync;
}

int
Electrostatic::get_no_composite ()
{
    return no_composite;
}

PArray<MultiFab>&
Electrostatic::get_grad_phi_prev (int level)
{
    return grad_phi_prev[level];
}

PArray<MultiFab>&
Electrostatic::get_grad_phi_curr (int level)
{
    return grad_phi_curr[level];
}

void
Electrostatic::plus_grad_phi_curr (int               level,
                             PArray<MultiFab>& addend)
{
    for (int n = 0; n < BL_SPACEDIM; n++)
        grad_phi_curr[level][n].plus(addend[n], 0, 1, 0);
}

void
Electrostatic::swap_time_levels (int level)
{
    for (int n=0; n < BL_SPACEDIM; n++)
    {
        MultiFab* dummy = grad_phi_curr[level].remove(n);
        grad_phi_prev[level].clear(n);
        grad_phi_prev[level].set(n, dummy);
        grad_phi_curr[level].set(n, new MultiFab(BoxArray(grids[level]).surroundingNodes(n), 1, 1));
        grad_phi_curr[level][n].setVal(1.e50,grad_phi_curr[level][n].nGrow());
    }
}

void
Electrostatic::zero_phi_flux_reg (int level)
{
    phi_flux_reg[level].setVal(0);
}

void
Electrostatic::solve_for_old_phi (int               level,
                                  MultiFab&         phi,
                                  PArray<MultiFab>& grad_phi,
                                  int               fill_interior)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        *gmsg << "Electrostatic ... single level solve for old phi at level " 
                  << level << endl;

    MultiFab Rhs(grids[level], 1, 0);
    Rhs.setVal(0.0);

    AddParticlesToRhs(level,Rhs,1);

    // We shouldn't need to use virtual or ghost particles for old phi solves.

    const Real time  = level_data[level].get_state_data(Elec_Field_Type).prevTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Electrostatic::solve_for_new_phi (int               level,
                                  MultiFab&         phi,
                                  PArray<MultiFab>& grad_phi,
                                  int               fill_interior,
                                  int               field_n_grow)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        *gmsg << "Electrostatic ... single level solve for new phi at level " 
                  << level << endl;

    MultiFab Rhs(grids[level], 1, 0);
    Rhs.setVal(0.0);

    AddParticlesToRhs(level,Rhs,field_n_grow);
    AddVirtualParticlesToRhs(level,Rhs,field_n_grow);
    AddGhostParticlesToRhs(level,Rhs,0);
    
    const Real time = level_data[level].get_state_data(Elec_Field_Type).curTime();
    solve_for_phi(level, Rhs, phi, grad_phi, time, fill_interior);
}

void
Electrostatic::solve_for_phi (int               level,
                              MultiFab&         Rhs,
                              MultiFab&         phi,
                              PArray<MultiFab>& grad_phi,
                              Real              time,
                              int               fill_interior)

{
    if (verbose && ParallelDescriptor::IOProcessor())
        *gmsg << " ... solve for phi at level " << level << '\n';

    const Real* dx = parent->Geom(level).CellSize();
#if 1
    Rhs.mult(-8.9875517879979115e+09);
#endif

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings) 
        ParallelDescriptor::Barrier();
    const Real strt = ParallelDescriptor::second();

//  *****************************************************************************
//  Sanity check
//  *****************************************************************************
#ifndef NDEBUG
    if (Rhs.contains_nan(0,1,0))
    {
        *gmsg << "Rhs in solve_for_phi at level " << level << " has NaNs" << endl;
        BoxLib::Abort("");
    }
#endif

//  *****************************************************************************
//  Set up coarse-fine stencil information
//  *****************************************************************************
    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0;
            xb[0][i] = 0;
        }
    }
    else
    {
        const Real* dx_crse = parent->Geom(level-1).CellSize();
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0.5 * dx_crse[i];
            xb[0][i] = 0.5 * dx_crse[i];
        }
    }

//  *****************************************************************************
//  Set up boundary conditions for phi for solve at level > 0
//  *****************************************************************************
    Accel* cs = dynamic_cast<Accel*>(&parent->getLevel(level));

    BL_ASSERT(cs != 0);

    const Geometry& geom = parent->Geom(level);
    MacBndry bndry(grids[level], 1, geom);

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                     : IntVect::TheZeroVector();
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    if (level == 0)
    {
        bndry.setBndryValues(phi, src_comp, dest_comp, num_comp, *phys_bc);
        phi.setVal(0.,phi.nGrow());
    }
    else
    {
        MultiFab c_phi;
        get_crse_phi(level, c_phi, time);
        BoxArray crse_boxes = BoxArray(grids[level]).coarsen(crse_ratio);
        const int in_rad     = 0;
        const int out_rad    = 1;
        const int extent_rad = 2;
        BndryRegister crse_br(crse_boxes, in_rad, out_rad, extent_rad,
                              num_comp);
        crse_br.copyFrom(c_phi, c_phi.nGrow(), src_comp, dest_comp, num_comp);
        bndry.setBndryValues(crse_br, src_comp, phi, src_comp, dest_comp,
                             num_comp, crse_ratio, *phys_bc);
    }

//  *****************************************************************************
//  Trilinos Solver only
//  *****************************************************************************

    // Need to pass multi-level rather than single-level MultiFabs to init_trilinos
    PArray<MultiFab> phi_p(1);
    PArray<MultiFab> rhs_p(1);

    phi_p.set(0,&phi);
    rhs_p.set(0,&Rhs);

    Solver* trilinos_solver = init_trilinos(rhs_p,phi_p,dx);

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
        ParallelDescriptor::Barrier();
    const Real strt_solve = ParallelDescriptor::second();

//  *****************************************************************************
//  Actual solve
//  *****************************************************************************
    trilinos_solver->Compute();
    trilinos_solver->CopySolution(phi_p);

    phi_p[0].FillBoundary();

    delete trilinos_solver;

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    Real end_solve = ParallelDescriptor::second() - strt_solve;
    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        ParallelDescriptor::ReduceRealMax(end_solve,IOProc);
        if (ParallelDescriptor::IOProcessor())
            *gmsg << "Electrostatic:: time in solve          = " << end_solve << '\n';
    }

//  *****************************************************************************
//  Get the fluxes on cell faces (*not* centers)
//  ***************************************************************************** 
    grad_phi[0].setVal(0.,grad_phi[0].nGrow());
    grad_phi[1].setVal(0.,grad_phi[1].nGrow());
    grad_phi[2].setVal(0.,grad_phi[2].nGrow());
    compute_face_fluxes(grids[level], grad_phi, phi, dx);

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

        ParallelDescriptor::ReduceRealMax(end, IOProc);

        if (ParallelDescriptor::IOProcessor())
            *gmsg << "Electrostatic::solve_for_phi() time = " << end << '\n';
    }
//  *****************************************************************************
//  END TIMINGS
//  *****************************************************************************
}

void
Electrostatic::solve_for_delta_phi (int                        crse_level,
                              int                        fine_level,
                              MultiFab&                  crse_rhs,
                              PArray<MultiFab>&          delta_phi,
                              PArray<PArray<MultiFab> >& grad_delta_phi)
{
  /*
#if 0
    const int num_levels = fine_level - crse_level + 1;

    BL_ASSERT(grad_delta_phi.size() == num_levels);
    BL_ASSERT(delta_phi.size() == num_levels);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "... solving for delta_phi at crse_level = " << crse_level
                  << '\n';
        std::cout << "...                    up to fine_level = " << fine_level
                  << '\n';
    }

    const Geometry& geom = parent->Geom(crse_level);
    MacBndry bndry(grids[crse_level], 1, geom);

    IntVect crse_ratio = crse_level > 0 ? parent->refRatio(crse_level-1)
                                          : IntVect::TheZeroVector();

    // Set homogeneous Dirichlet values for the solve.
    bndry.setHomogValues(*phys_bc, crse_ratio);

    std::vector<BoxArray>            bav(num_levels);
    std::vector<DistributionMapping> dmv(num_levels);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        bav[lev-crse_level] = grids[lev];
        MultiFab& G_new = level_data[lev].get_new_data(Elec_Field_Type);
        dmv[lev-crse_level] = G_new.DistributionMap();
    }
    std::vector<Geometry> fgeom(num_levels);
    for (int lev = crse_level; lev <= fine_level; lev++)
        fgeom[lev-crse_level] = parent->Geom(lev);

    MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type);

    Array< Array<Real> > xa(num_levels);
    Array< Array<Real> > xb(num_levels);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        xa[lev-crse_level].resize(BL_SPACEDIM);
        xb[lev-crse_level].resize(BL_SPACEDIM);
        if (lev == 0)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev-crse_level][i] = 0;
                xb[lev-crse_level][i] = 0;
            }
        }
        else
        {
            const Real* dx_crse = parent->Geom(lev-1).CellSize();
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev-crse_level][i] = 0.5 * dx_crse[i];
                xb[lev-crse_level][i] = 0.5 * dx_crse[i];
            }
        }
    }

    MultiFab** phi_p = new MultiFab*[num_levels];
    MultiFab** Rhs_p = new MultiFab*[num_levels];

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        // Turn over control of pointer for a moment
        phi_p[lev-crse_level] = delta_phi.remove(lev-crse_level);
        phi_p[lev-crse_level]->setVal(0);

        if (lev == crse_level)
        {
            Rhs_p[0] = &crse_rhs;
        }
        else
        {
            Rhs_p[lev-crse_level] = new MultiFab(grids[lev], 1, 0);
            Rhs_p[lev-crse_level]->setVal(0);
        }

    }

    mgt_solver.set_const_gravity_coeffs(xa, xb);

    const Real tol     = delta_tol;
    Real       abs_tol = level_solver_resnorm[crse_level];
    for (int lev = crse_level + 1; lev < fine_level; lev++)
        abs_tol = std::max(abs_tol,level_solver_resnorm[lev]);

    Real final_resnorm = 0;
    mgt_solver.solve(phi_p, Rhs_p, tol, abs_tol, bndry, 1, final_resnorm);

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        PArray<MultiFab>& gdphi = grad_delta_phi[lev-crse_level];
        const Real* dx = parent->Geom(lev).CellSize();
        mgt_solver.get_fluxes(lev-crse_level, gdphi, dx);
    }

    for (int lev = 0; lev < num_levels; lev++)
    {
        delta_phi.set(lev, phi_p[lev]);  // return control of pointer
        if (lev != 0)
            delete Rhs_p[lev]; // Do not delete the [0] Rhs, it is passed in
    }
    delete [] phi_p;
    delete [] Rhs_p;
#endif
  */
}

void
Electrostatic::e_field_sync (int               crse_level,
                       int               fine_level,
                       const MultiFab&   drho_and_drhoU,
                       const MultiFab&   dphi,
                       PArray<MultiFab>& grad_delta_phi_cc)
{
    BL_ASSERT(parent->finestLevel()>crse_level);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        *gmsg << " ... gravity_sync at crse_level " << crse_level << '\n'
              << " ...     up to finest_level     " << fine_level << '\n';
    }
    /*#if 0

    // Build Rhs for solve for delta_phi
    MultiFab crse_rhs(grids[crse_level], 1, 0);
    MultiFab::Copy(crse_rhs, drho_and_drhoU, 0, 0, 1, 0);
    crse_rhs.plus(dphi, 0, 1, 0);

    const Geometry& crse_geom   = parent->Geom(crse_level);
    const Box&      crse_domain = crse_geom.Domain();
    if (crse_geom.isAllPeriodic()
        && (grids[crse_level].numPts() == crse_domain.numPts()))
    {
        Real local_correction = 0;
        for (MFIter mfi(crse_rhs); mfi.isValid(); ++mfi)
            local_correction += crse_rhs[mfi].sum(mfi.validbox(), 0, 1);
        ParallelDescriptor::ReduceRealSum(local_correction);

        local_correction /= grids[crse_level].numPts();

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "WARNING: Adjusting RHS in gravity_sync solve by "
                      << local_correction << '\n';
        for (MFIter mfi(crse_rhs); mfi.isValid(); ++mfi)
            crse_rhs[mfi].plus(-local_correction);
    }

    // delta_phi needs a ghost cell for the solve
    PArray<MultiFab>  delta_phi(fine_level - crse_level + 1, PArrayManage);
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        delta_phi.set(lev - crse_level, new MultiFab(grids[lev], 1, 1));
        delta_phi[lev-crse_level].setVal(0);
    }

    PArray<PArray<MultiFab> > ec_gdPhi(fine_level - crse_level + 1,
                                       PArrayManage);
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        ec_gdPhi.set(lev-crse_level, new PArray<MultiFab>(BL_SPACEDIM,
                                                          PArrayManage));
        for (int n = 0; n < BL_SPACEDIM; ++n)
            ec_gdPhi[lev-crse_level].set(n, new MultiFab(BoxArray(grids[lev]).surroundingNodes(n), 1, 0));
    }

    // Do multi-level solve for delta_phi
    solve_for_delta_phi(crse_level, fine_level, crse_rhs, delta_phi, ec_gdPhi);

    crse_rhs.clear();

    // In the all-periodic case we enforce that delta_phi averages to zero.
    if (crse_geom.isAllPeriodic()
        && (grids[crse_level].numPts() == crse_domain.numPts()))
    {
        Real local_correction = 0;
        for (MFIter mfi(delta_phi[0]); mfi.isValid(); ++mfi)
            local_correction += delta_phi[0][mfi].sum(mfi.validbox(), 0, 1);
        ParallelDescriptor::ReduceRealSum(local_correction);

        local_correction = local_correction / grids[crse_level].numPts();

        for (int lev = crse_level; lev <= fine_level; lev++)
            for (MFIter mfi(delta_phi[lev-crse_level]); mfi.isValid(); ++mfi)
                delta_phi[lev-crse_level][mfi].plus(-local_correction);
    }

    // Add delta_phi to phi_curr, and grad(delta_phi) to grad(delta_phi_curr) on
    // each level
    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        level_data[lev].get_new_data(PhiGrav_Type).plus(delta_phi[lev-crse_level], 0, 1, 0);
        for (int n = 0; n < BL_SPACEDIM; n++)
            grad_phi_curr[lev][n].plus(ec_gdPhi[lev-crse_level][n], 0, 1, 0);
    }

    int is_new = 1;

    // Average phi_curr from fine to coarse level
    for (int lev = fine_level - 1; lev >= crse_level; lev--)
    {
        const IntVect ratio = parent->refRatio(lev);
        average_down(level_data[lev  ].get_new_data(PhiGrav_Type),
                     level_data[lev+1].get_new_data(PhiGrav_Type),ratio);
    }

    // Average the edge-based grad_phi from finer to coarser level
    for (int lev = fine_level-1; lev >= crse_level; lev--)
        average_fine_ec_onto_crse_ec(lev, is_new);

    // Add the contribution of grad(delta_phi) to the flux register below if
    // necessary.
    if (crse_level > 0)
    {
        for (MFIter mfi(delta_phi[0]); mfi.isValid(); ++mfi)
            for (int n = 0; n < BL_SPACEDIM; ++n)
                phi_flux_reg[crse_level].FineAdd(ec_gdPhi[0][n][mfi], n,
                                                 mfi.index(), 0, 0, 1, 1);
    }

    int lo_bc[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        lo_bc[dir] = (phys_bc->lo(dir) == Symmetry);
    }
    int symmetry_type = Symmetry;
    int coord_type    = Geometry::Coord();

    for (int lev = crse_level; lev <= fine_level; lev++)
    {
        const Real* problo = parent->Geom(lev).ProbLo();
        const Real* dx     = parent->Geom(lev).CellSize();
        grad_delta_phi_cc[lev-crse_level].setVal(0);

        for (MFIter mfi(grad_delta_phi_cc[lev-crse_level]); mfi.isValid(); ++mfi)
        {
            const int  index = mfi.index();
            const Box& bx    = grids[lev][index];

            // Average edge-centered gradients of crse dPhi to cell centers
            BL_FORT_PROC_CALL(FORT_AVG_EC_TO_CC, fort_avg_ec_to_cc)
                (bx.loVect(), bx.hiVect(), lo_bc, &symmetry_type,
                 BL_TO_FORTRAN(grad_delta_phi_cc[lev-crse_level][mfi]),
                 D_DECL(BL_TO_FORTRAN(ec_gdPhi[lev-crse_level][0][mfi]),
                        BL_TO_FORTRAN(ec_gdPhi[lev-crse_level][1][mfi]),
                        BL_TO_FORTRAN(ec_gdPhi[lev-crse_level][2][mfi])),
                 problo);
        }
    }
#endif
    */
}

void
Electrostatic::get_crse_phi (int       level,
                       MultiFab& phi_crse,
                       Real      time)
{
    BL_ASSERT(level != 0);
    /*
#if 0

    const Real t_old = level_data[level-1].get_state_data(Elec_Field_Type).prevTime();
    const Real t_new = level_data[level-1].get_state_data(Elec_Field_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);

    phi_crse.clear();
    phi_crse.define(grids[level-1], 1, 1, Fab_allocate);

    // BUT NOTE we don't trust phi's ghost cells.
    FArrayBox phi_crse_temp;

    // Note that we must do these cases separately because it's possible to do a 
    //   new solve after a regrid when the old data on the coarse grid may not yet
    //   be defined.
    for (MFIter mfi(phi_crse); mfi.isValid(); ++mfi)
    {
        phi_crse_temp.resize(phi_crse[mfi].box(),1);

        if (fabs(alpha-1.0) < 1.e-15)
        {
            phi_crse[mfi].copy(level_data[level-1].get_new_data(PhiGrav_Type)[mfi]);
        } 
        else if (fabs(alpha) < 1.e-15)
        {
            phi_crse[mfi].copy(level_data[level-1].get_old_data(PhiGrav_Type)[mfi]);
        } 
        else 
        {
            phi_crse_temp.copy(level_data[level-1].get_old_data(PhiGrav_Type)[mfi]);
            Real omalpha = 1.0 - alpha;
            phi_crse_temp.mult(omalpha);
            phi_crse[mfi].copy(level_data[level-1].get_new_data(PhiGrav_Type)[mfi]);
            phi_crse[mfi].mult(alpha);
            phi_crse[mfi].plus(phi_crse_temp);
        } 
    }
    phi_crse_temp.clear();

    phi_crse.FillBoundary();

    const Geometry& geom = parent->Geom(level-1);
    geom.FillPeriodicBoundary(phi_crse, true);
#endif
    */
}

void
Electrostatic::get_crse_grad_phi (int               level,
                                  PArray<MultiFab>& grad_phi_crse,
                                  Real              time)
{
    BL_ASSERT(level!=0);

    const Real t_old = level_data[level-1].get_state_data(Elec_Field_Type).prevTime();
    const Real t_new = level_data[level-1].get_state_data(Elec_Field_Type).curTime();
    const Real alpha = (time - t_old) / (t_new - t_old);

    FArrayBox grad_phi_crse_temp;

    BL_ASSERT(grad_phi_crse.size() == BL_SPACEDIM);

    for (int i = 0; i < BL_SPACEDIM; ++i)
    {
        BL_ASSERT(!grad_phi_crse.defined(i));
        const BoxArray eba = BoxArray(grids[level-1]).surroundingNodes(i);
        grad_phi_crse.set(i, new MultiFab(eba, 1, 0));

        for (MFIter mfi(grad_phi_crse[i]); mfi.isValid(); ++mfi)
        {
            grad_phi_crse_temp.resize(mfi.validbox(),1);

            grad_phi_crse_temp.copy(grad_phi_prev[level-1][i][mfi]);
            const Real omalpha = 1.0 - alpha;
            grad_phi_crse_temp.mult(omalpha);

            grad_phi_crse[i][mfi].copy(grad_phi_curr[level-1][i][mfi]);
            grad_phi_crse[i][mfi].mult(alpha);
            grad_phi_crse[i][mfi].plus(grad_phi_crse_temp);
        }
        grad_phi_crse[i].FillBoundary();

        const Geometry& geom = parent->Geom(level-1);
//         geom.FillPeriodicBoundary(grad_phi_crse[i], false);
    }
}

void
Electrostatic::multilevel_solve_for_new_phi (int level,
                                             int finest_level,
                                             int use_previous_phi_as_guess)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        *gmsg << "Electrostatic ... multilevel solve for new phi at base level " << level
              << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_curr[lev].size()==BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            grad_phi_curr[lev].clear(n);
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_curr[lev].set(n, new MultiFab(eba, 1, 1));
        }
    }

    int is_new = 1;
    actual_multilevel_solve(level, finest_level, grad_phi_curr,
                            is_new, use_previous_phi_as_guess);
}

void
Electrostatic::multilevel_solve_for_old_phi (int level,
                                             int finest_level,
                                             int use_previous_phi_as_guess)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        *gmsg << "Electrostatic ... multilevel solve for old phi at base level " << level
              << " to finest level " << finest_level << '\n';

    for (int lev = level; lev <= finest_level; lev++)
    {
        BL_ASSERT(grad_phi_prev[lev].size() == BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            grad_phi_prev[lev].clear(n);
            const BoxArray eba = BoxArray(grids[lev]).surroundingNodes(n);
            grad_phi_prev[lev].set(n, new MultiFab(eba, 1, 1));
        }
    }

    int is_new = 0;
    actual_multilevel_solve(level, finest_level, grad_phi_prev,
                            is_new, use_previous_phi_as_guess);
}

void
Electrostatic::multilevel_solve_for_phi(int level, int finest_level, int use_previous_phi_as_guess)
{
    multilevel_solve_for_new_phi(level, finest_level, use_previous_phi_as_guess);
}

void
Electrostatic::actual_multilevel_solve (int                 level,
                                  int                       finest_level,
                                  Array<PArray<MultiFab> >& grad_phi,
                                  int                       is_new,
                                  int                       use_previous_phi_as_guess)
{
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const Real      strt = ParallelDescriptor::second();

    const int num_levels = finest_level - level + 1;

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
        ParallelDescriptor::Barrier();
    const Real strt_rhs = ParallelDescriptor::second();

//  *****************************************************************************
//   Create RHS based on particle locations and charge
//  *****************************************************************************
    PArray<MultiFab> Rhs_particles(num_levels,PArrayManage);
    for (int lev = 0; lev < num_levels; lev++)
    { 
       Rhs_particles.set(lev, new MultiFab(grids[level+lev], 1, 0));
       Rhs_particles[lev].setVal(0.);
    } 

    *gmsg << "Total: " << Accel::thePAPC()->TotalNumberOfParticles() << endl
          << "Level: " << Accel::thePAPC()->NumberOfParticlesAtLevel(0) << endl;

    AddParticlesToRhs(level,finest_level,Rhs_particles);

    AddGhostParticlesToRhs(level,Rhs_particles);

    AddVirtualParticlesToRhs(finest_level,Rhs_particles);

#ifdef AMR_DUMMY_SOLVE
    for (int lev = 0; lev < num_levels; ++lev)
        Rhs_particles[lev].setVal(-1.0);
#endif
    
#ifdef DBG_SCALARFIELD
    for (int l = 0; l < num_levels; ++l) {
        for (MFIter mfi(Rhs_particles[l]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            FArrayBox& fab = Rhs_particles[l][mfi];

            for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
                if ( proc == ParallelDescriptor::MyProc() ) {
                    std::string outfile = "data/amr-rho_scalar-level-" + std::to_string(l);
                    
                    std::ofstream out;
                    
                    if ( proc == 0 )
                        out.open(outfile);
                    else
                        out.open(outfile, std::ios_base::app);
                    
                    if ( !out.is_open() )
                        throw OpalException("Error in Electrostatic::actual_multilevel_solve",
                                            "Couldn't open the file: " + outfile);
                    
                    for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                        for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                            for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                                IntVect ivec(i, j, k);
                                out << i + 1 << " " << j + 1 << " " << k + 1 << " " << fab(ivec, 0) << " " << proc
                                    << std::endl;
                            }
                        }
                    }
                    out.close();
                }
                ParallelDescriptor::Barrier();
            }
        }
    }
#endif

#ifndef AMR_DUMMY_SOLVE
    // rhs: /*-*/ \rho / (4 * pi * epsilon_0 )
    for (int lev = 0; lev < num_levels; lev++)
        Rhs_particles[lev].mult(-8.9875517879979115e+09, 0, 1);  //opal_coupling ( 4 * pi * epsilon_0 )
#endif
    
//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
    {
        ParallelDescriptor::Barrier();
        Real      end    = ParallelDescriptor::second() - strt_rhs;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            *gmsg << "Electrostatic:: time in making rhs        = " << end << '\n';
    }

//  *****************************************************************************
//  Set up coarse-fine stencil information
//  *****************************************************************************
    Array< Array<Real> > xa(num_levels);
    Array< Array<Real> > xb(num_levels);

    for (int lev = 0; lev < num_levels; lev++)
    {
        xa[lev].resize(BL_SPACEDIM);
        xb[lev].resize(BL_SPACEDIM);
        if (level + lev == 0)
        {
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0;
                xb[lev][i] = 0;
            }
        }
        else
        {
            const Real* dx_crse = parent->Geom(level + lev - 1).CellSize();
            for (int i = 0; i < BL_SPACEDIM; ++i)
            {
                xa[lev][i] = 0.5 * dx_crse[i];
                xb[lev][i] = 0.5 * dx_crse[i];
            }
        }
    }

//  *****************************************************************************
//  We will work in result data structure directly
//  *****************************************************************************
    PArray<MultiFab> phi_p(num_levels);

    for (int lev = 0; lev < num_levels; lev++)
    {
        if (is_new == 0)
        {
           phi_p.set(lev, &level_data[level+lev].get_old_data(Elec_Potential_Type));
        }
        else
        {
           // Working in result data structure directly
           phi_p.set(lev, &level_data[level+lev].get_new_data(Elec_Potential_Type));
        }

        if (!use_previous_phi_as_guess)
           phi_p[lev].setVal(0.,phi_p[lev].nGrow());

    }
     
    // Average phi from fine to coarse level before the solve.
    for (int lev = num_levels-1; lev > 0; lev--)
    {
        const IntVect ratio = parent->refRatio(level+lev-1);
        average_down(phi_p[lev-1], phi_p[lev], ratio);
    }

//  *****************************************************************************
//  Trilinos Solver only
//  *****************************************************************************

    if (level != 0 && finest_level !=0)    
       BoxLib::Error("This only works for single level right now");

    const Real* dx = parent->Geom(level).CellSize();
    
    *gmsg << "dx = " << dx[0] << ", dy = " << dx[1] << ", dz = " << dx[2] << endl;
    
    Solver* trilinos_solver = init_trilinos(Rhs_particles,phi_p,dx);
    
//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
        ParallelDescriptor::Barrier();
    const Real strt_solve = ParallelDescriptor::second();
    
//  *****************************************************************************
//  Actual solve
//  *****************************************************************************
    trilinos_solver->Compute();
    
    trilinos_solver->CopySolution(phi_p);
    
    
    // Average phi from fine to coarse level after the solve since the coarse grid
    //   values under the fine grid have not been defined.
    for (int lev = num_levels-1; lev > 0; lev--)
    {
        const IntVect ratio = parent->refRatio(level+lev-1);
        average_down(phi_p[lev-1], phi_p[lev], ratio);
    }
    for (int i = 0; i < num_levels; i++)
        phi_p[i].FillBoundary();

    if (phi_p[0].norm0() > 1.e100)
       for (MFIter mfi(phi_p[0]); mfi.isValid(); ++mfi)
           *gmsg << "BAD PHI " << phi_p[0][mfi] << endl;
       
    delete trilinos_solver;

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    Real end_solve = ParallelDescriptor::second() - strt_solve;
    if (show_timings)
    {
        ParallelDescriptor::ReduceRealMax(end_solve,IOProc);
        if (ParallelDescriptor::IOProcessor())
            *gmsg << "Electrostatic:: time in solve          = " << end_solve << '\n';
    }

//  *****************************************************************************
//  Get the fluxes on cell faces (*not* centers)
//  *****************************************************************************
    for (int lev = 0; lev < num_levels; lev++)
    {
        const Real* dx = parent->Geom(level+lev).CellSize();
        compute_face_fluxes(grids[level+lev], grad_phi[level+lev], phi_p[lev], dx);
    }

    // Average phi from fine to coarse level
    for (int lev = finest_level; lev > level; lev--)
    {
        const IntVect ratio = parent->refRatio(lev-1);
        if (is_new == 1)
        {
            average_down(level_data[lev-1].get_new_data(Elec_Potential_Type),
                         level_data[lev  ].get_new_data(Elec_Potential_Type), ratio);
        }
        else if (is_new == 0)
        {
            average_down(level_data[lev-1].get_old_data(Elec_Potential_Type),
                         level_data[lev  ].get_old_data(Elec_Potential_Type), ratio);
        }
    }

    // Average grad_phi from fine to coarse level
    for (int lev = finest_level; lev > level; lev--)
        average_fine_ec_onto_crse_ec(lev-1,is_new);

//  *****************************************************************************
//  FOR TIMINGS
//  *****************************************************************************
    if (show_timings)
    {
        Real      end    = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
        {
            *gmsg << "Electrostatic:: all                  : time = " << end << '\n'
                  << "Electrostatic:: all but solve        : time = " << end - end_solve << '\n';
        }
    }
//  *****************************************************************************
//  END TIMINGS
//  *****************************************************************************

    if (ParallelDescriptor::IOProcessor())
        *gmsg << "Leaving multilevel solver now ... " << endl;
}

void
Electrostatic::get_old_e_field (int       level,
                                MultiFab& e_field,
                                Real      time)
{
    // Fill grow cells in grad_phi, will need to compute grad_phi_cc in 1 grow cell
    const Geometry& geom = parent->Geom(level);
    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        {
            grad_phi_prev[level][i].setBndry(0);
            grad_phi_prev[level][i].FillBoundary();
//             geom.FillPeriodicBoundary(grad_phi_prev[level][i]);
        }
    }
    else
    {
        PArray<MultiFab> crse_grad_phi(BL_SPACEDIM, PArrayManage);
        get_crse_grad_phi(level, crse_grad_phi, time);
        fill_ec_grow(level, grad_phi_prev[level], crse_grad_phi);
    }

    int lo_bc[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        lo_bc[dir] = phys_bc->lo(dir);

    int         symmetry_type  = Symmetry;
    const Real* problo         = parent->Geom(level).ProbLo();

    // Average edge-centered gradients to cell centers, including grow cells
    //   Grow cells are filled either by physical bc's in AVG_EC_TO_CC or
    //   by FillBoundary call to e_field afterwards.
    for (MFIter mfi(e_field); mfi.isValid(); ++mfi)
    {
        const int  i  = mfi.index();
        const Box& bx = grids[level][i];

        BL_FORT_PROC_CALL(FORT_AVG_EC_TO_CC, fort_avg_ec_to_cc)
            (bx.loVect(), bx.hiVect(), lo_bc, &symmetry_type,
             BL_TO_FORTRAN(e_field[i]),
             D_DECL(BL_TO_FORTRAN(grad_phi_prev[level][0][i]),
                    BL_TO_FORTRAN(grad_phi_prev[level][1][i]),
                    BL_TO_FORTRAN(grad_phi_prev[level][2][i])),
             problo);
    }

    e_field.FillBoundary();
//     geom.FillPeriodicBoundary(e_field, 0, BL_SPACEDIM);

    MultiFab& G_old = level_data[level].get_old_data(Elec_Field_Type);
    MultiFab::Copy(G_old, e_field, 0, 0, BL_SPACEDIM, 0);

    // This is a hack-y way to fill the ghost cell values of e_field
    //   before returning it
    AmrLevel* amrlev = &parent->getLevel(level);
 
    for (FillPatchIterator fpi(*amrlev,G_old,e_field.nGrow(),time,Elec_Field_Type,0,BL_SPACEDIM);
         fpi.isValid(); ++fpi)
      {
         int i = fpi.index();
         e_field[i].copy(fpi());
      }
}

void
Electrostatic::get_new_e_field (int       level,
                                MultiFab& e_field,
                                Real      time)
{
    /// Fill grow cells in `grad_phi`, will need to compute `grad_phi_cc` in
    // 1 grow cell
    const Geometry& geom = parent->Geom(level);
    if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM ; i++)
        {
            grad_phi_curr[level][i].setBndry(0);
            grad_phi_curr[level][i].FillBoundary();
//             geom.FillPeriodicBoundary(grad_phi_curr[level][i]);
        }
    }
    else
    {
        PArray<MultiFab> crse_grad_phi(BL_SPACEDIM, PArrayManage);
        get_crse_grad_phi(level, crse_grad_phi, time);
        fill_ec_grow(level, grad_phi_curr[level], crse_grad_phi);
    }

    int lo_bc[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
        lo_bc[dir] = phys_bc->lo(dir);

    int         symmetry_type = Symmetry;
    const Real* problo        = parent->Geom(level).ProbLo();

    // Average edge-centered gradients to cell centers, including grow cells
    //   Grow cells are filled either by physical bc's in AVG_EC_TO_CC or
    //   by FillBoundary call to e_field afterwards.
    for (MFIter mfi(e_field); mfi.isValid(); ++mfi)
    {
        const int  i  = mfi.index();
        const Box& bx = grids[level][i];

        BL_FORT_PROC_CALL(FORT_AVG_EC_TO_CC, fort_avg_ec_to_cc)
            (bx.loVect(), bx.hiVect(), lo_bc, &symmetry_type,
             BL_TO_FORTRAN(e_field[i]),
             D_DECL(BL_TO_FORTRAN(grad_phi_curr[level][0][i]),
                    BL_TO_FORTRAN(grad_phi_curr[level][1][i]), BL_TO_FORTRAN(grad_phi_curr[level][2][i])),
             problo);
    }

    e_field.FillBoundary();
//     geom.FillPeriodicBoundary(e_field, 0, BL_SPACEDIM);

    MultiFab& E_new = level_data[level].get_new_data(Elec_Field_Type);
    MultiFab::Copy(E_new, e_field, 0, 0, BL_SPACEDIM, 0);

#ifdef DBG_SCALARFIELD
    /* Writes electric field vector to a file.
     * Format: x, y, z, ex, ey, ez, processor
     */
    for (MFIter mfi(E_new); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox& fab = E_new[mfi];
        
        for (int proc = 0; proc < ParallelDescriptor::NProcs(); ++proc) {
            
            if ( proc == ParallelDescriptor::MyProc() ) {
                std::string outfile = "data/amr-e_field-level-" + std::to_string(level);
                
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                   out.open(outfile, std::ios_base::app);
                
                if ( !out.is_open() )
                    throw OpalException("Error in Electrostatic::get_new_e_field",
                                        "Couldn't open the file: " + outfile);
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            out << i + 1 << " " << j + 1 << " " << k + 1
                                << " " << fab(ivec, 0) << " " << fab(ivec, 1) << " " << fab(ivec, 2) << " "
                                << proc << std::endl;
                        }
                    }
                }
                out.close();
            }
            ParallelDescriptor::Barrier();
        }
    }
#endif

    // This is a hack-y way to fill the ghost cell values of e_field
    //   before returning it
    AmrLevel* amrlev = &parent->getLevel(level) ;

    for (FillPatchIterator fpi(*amrlev,E_new,e_field.nGrow(),time,Elec_Field_Type,0,BL_SPACEDIM);
         fpi.isValid(); ++fpi)
      {
         int i = fpi.index();
         e_field[i].copy(fpi());
      }
}

void
Electrostatic::add_to_fluxes(int level, int iteration, int ncycle)
{
    const int       finest_level      = parent->finestLevel();
    FluxRegister*   phi_fine          = (level<finest_level ? &phi_flux_reg[level+1] : 0);
    FluxRegister*   phi_current       = (level>0 ? &phi_flux_reg[level] : 0);
    const Geometry& geom              = parent->Geom(level);
    const Real*     dx                = geom.CellSize();
    const Real      area[BL_SPACEDIM] = { dx[1]*dx[2], dx[0]*dx[2], dx[0]*dx[1] };

    if (phi_fine)
    {
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            BoxArray ba = grids[level];
            ba.surroundingNodes(n);
            MultiFab fluxes(ba, 1, 0);
            for (MFIter mfi(fluxes); mfi.isValid(); ++mfi)
            {
                FArrayBox& gphi_flux = fluxes[mfi];
                gphi_flux.copy(grad_phi_curr[level][n][mfi]);
                gphi_flux.mult(area[n]);
            }

            phi_fine->CrseInit(fluxes, n, 0, 0, 1, -1);
        }
    }

    if (phi_current && (iteration == ncycle))
    {
        for (int n = 0; n < BL_SPACEDIM; ++n)
        {
            phi_current->FineAdd(grad_phi_curr[level][n], n, 0, 0, 1, area[n]);
        }
    }
}

void
Electrostatic::average_fine_ec_onto_crse_ec(int level, int is_new)
{
    //
    // NOTE: this is called with level == the coarser of the two levels involved.
    //
    if (level == parent->finestLevel()) return;
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_gphi_fine_BA(grids[level+1].size());

    IntVect  fine_ratio = parent->refRatio(level);

    for (int i = 0; i < crse_gphi_fine_BA.size(); ++i)
        crse_gphi_fine_BA.set(i, BoxLib::coarsen(grids[level+1][i],
                                                 fine_ratio));

    PArray<MultiFab> crse_gphi_fine(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        const BoxArray eba = BoxArray(crse_gphi_fine_BA).surroundingNodes(n);
        crse_gphi_fine.set(n, new MultiFab(eba, 1, 0));
    }

    if (is_new == 1)
    {
        for (MFIter mfi(grad_phi_curr[level+1][0]); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            const Box& overlap = crse_gphi_fine_BA[i];

            BL_FORT_PROC_CALL(FORT_AVERAGE_EC, fort_average_ec)
                (D_DECL(BL_TO_FORTRAN(grad_phi_curr[level+1][0][mfi]),
                        BL_TO_FORTRAN(grad_phi_curr[level+1][1][mfi]),
                        BL_TO_FORTRAN(grad_phi_curr[level+1][2][mfi])),
                 D_DECL(BL_TO_FORTRAN(crse_gphi_fine[0][mfi]),
                        BL_TO_FORTRAN(crse_gphi_fine[1][mfi]),
                        BL_TO_FORTRAN(crse_gphi_fine[2][mfi])),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());
        }

        for (int n = 0; n < BL_SPACEDIM; ++n)
            grad_phi_curr[level][n].copy(crse_gphi_fine[n]);
    }
    else if (is_new == 0)
    {
        for (MFIter mfi(grad_phi_prev[level+1][0]); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            const Box& overlap = crse_gphi_fine_BA[i];

            BL_FORT_PROC_CALL(FORT_AVERAGE_EC, fort_average_ec)
                (D_DECL(BL_TO_FORTRAN(grad_phi_prev[level+1][0][mfi]),
                        BL_TO_FORTRAN(grad_phi_prev[level+1][1][mfi]),
                        BL_TO_FORTRAN(grad_phi_prev[level+1][2][mfi])),
                 D_DECL(BL_TO_FORTRAN(crse_gphi_fine[0][mfi]),
                        BL_TO_FORTRAN(crse_gphi_fine[1][mfi]),
                        BL_TO_FORTRAN(crse_gphi_fine[2][mfi])),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());
        }

        for (int n = 0; n < BL_SPACEDIM; ++n)
            grad_phi_prev[level][n].copy(crse_gphi_fine[n]);
    }
}

void
Electrostatic::average_down (MultiFab&       crse,
                       const MultiFab& fine,
                       const IntVect&  ratio)
{
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_fine_BA(fine.boxArray().size());

    for (int i = 0; i < fine.boxArray().size(); ++i)
    {
        crse_fine_BA.set(i, BoxLib::coarsen(fine.boxArray()[i], ratio));
    }

    MultiFab crse_fine(crse_fine_BA, 1, 0);

    for (MFIter mfi(fine); mfi.isValid(); ++mfi)
    {
        const int        i        = mfi.index();
        const Box&       overlap  = crse_fine_BA[i];
        FArrayBox&       crse_fab = crse_fine[i];
        const FArrayBox& fine_fab = fine[i];

        BL_FORT_PROC_CALL(FORT_AVGDOWN_PHI, fort_avgdown_phi)
            (BL_TO_FORTRAN(crse_fab), BL_TO_FORTRAN(fine_fab),
             overlap.loVect(), overlap.hiVect(), ratio.getVect());
    }

    crse.copy(crse_fine);
}

void
Electrostatic::reflux_phi (int       level,
                     MultiFab& dphi)
{
    const Geometry& geom = parent->Geom(level);
    dphi.setVal(0);
    phi_flux_reg[level+1].Reflux(dphi, 1.0, 0, 0, 1, geom);
}

void
Electrostatic::fill_ec_grow (int                     level,
                       PArray<MultiFab>&       ecF,
                       const PArray<MultiFab>& ecC) const
{
    //
    // Fill grow cells of the edge-centered mfs.  Assume
    // ecF built on edges of grids at this amr level, and ecC
    // is build on edges of the grids at amr level-1
    //
    BL_ASSERT(ecF.size() == BL_SPACEDIM);

    const int nGrow = ecF[0].nGrow();

    if (nGrow == 0 || level == 0) return;

    BL_ASSERT(nGrow == ecF[1].nGrow());
    BL_ASSERT(nGrow == ecF[2].nGrow());

    const BoxArray& fgrids = grids[level];
    const Geometry& fgeom  = parent->Geom(level);

    BoxList bl = BoxLib::GetBndryCells(fgrids, 1);

    BoxArray f_bnd_ba(bl);

    bl.clear();

    BoxArray c_bnd_ba   = BoxArray(f_bnd_ba.size());
    IntVect  crse_ratio = parent->refRatio(level-1);

    for (int i = 0; i < f_bnd_ba.size(); ++i)
    {
        c_bnd_ba.set(i, Box(f_bnd_ba[i]).coarsen(crse_ratio));
        f_bnd_ba.set(i, Box(c_bnd_ba[i]).refine(crse_ratio));
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        //
        // crse_src & fine_src must have same parallel distribution.
        // We'll use the KnapSack distribution for the fine_src_ba.
        // Since fine_src_ba should contain more points, this'll lead
        // to a better distribution.
        //
        BoxArray crse_src_ba(c_bnd_ba);
        BoxArray fine_src_ba(f_bnd_ba);

        crse_src_ba.surroundingNodes(n);
        fine_src_ba.surroundingNodes(n);

        std::vector<long> wgts(fine_src_ba.size());

        for (unsigned int i = 0; i < wgts.size(); i++)
        {
            wgts[i] = fine_src_ba[i].numPts();
        }
        DistributionMapping dm;
        //
        // This call doesn't invoke the MinimizeCommCosts() stuff.
        // There's very little to gain with these types of coverings
        // of trying to use SFC or anything else.
        // This also guarantees that these DistributionMaps won't be put into the
        // cache, as it's not representative of that used for more
        // usual MultiFabs.
        //
        dm.KnapSackProcessorMap(wgts, ParallelDescriptor::NProcs());

        MultiFab crse_src; crse_src.define(crse_src_ba, 1, 0, dm, Fab_allocate);
        MultiFab fine_src; fine_src.define(fine_src_ba, 1, 0, dm, Fab_allocate);

        crse_src.setVal(1.e200);
        fine_src.setVal(1.e200);

        //
        // We want to fill crse_src from ecC[n].
        // Gotta do it in steps since parallel copy only does valid region.
        //
        {
            BoxArray edge_grids = ecC[n].boxArray();
            edge_grids.grow(ecC[n].nGrow());
            MultiFab ecCG(edge_grids, 1, 0);

            for (MFIter mfi(ecC[n]); mfi.isValid(); ++mfi)
                ecCG[mfi].copy(ecC[n][mfi]);

            crse_src.copy(ecCG);
        }

        for (MFIter mfi(crse_src); mfi.isValid(); ++mfi)
        {
            const int nComp = 1;
            const Box box = crse_src[mfi].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_PC_EDGE_INTERP, fort_pc_edge_interp)
                (box.loVect(), box.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(crse_src[mfi]), BL_TO_FORTRAN(fine_src[mfi]));
        }

        crse_src.clear();
        //
        // Replace pc-interpd fine data with preferred u_mac data at
        // this level u_mac valid only on surrounding faces of valid
        // region - this op will not fill grow region.
        //
        fine_src.copy(ecF[n]); // parallel copy

        for (MFIter mfi(fine_src); mfi.isValid(); ++mfi)
        {
            //
            // Interpolate unfilled grow cells using best data from
            // surrounding faces of valid region, and pc-interpd data
            // on fine edges overlaying coarse edges.
            //
            const int nComp = 1;
            const Box& fbox = fine_src[mfi.index()].box();
            const int* rat = crse_ratio.getVect();
            BL_FORT_PROC_CALL(FORT_EDGE_INTERP, fort_edge_interp)
                (fbox.loVect(), fbox.hiVect(), &nComp, rat, &n,
                 BL_TO_FORTRAN(fine_src[mfi]));
        }
        //
        // Build a mf with no grow cells on ecF grown boxes, do parallel copy into.
        //
        BoxArray fgridsG = ecF[n].boxArray();
        fgridsG.grow(ecF[n].nGrow());

        MultiFab ecFG(fgridsG, 1, 0);

        ecFG.copy(fine_src); // Parallel copy

        fine_src.clear();

        ecFG.copy(ecF[n]);   // Parallel copy

        for (MFIter mfi(ecF[n]); mfi.isValid(); ++mfi)
            ecF[n][mfi].copy(ecFG[mfi]);
    }

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        ecF[n].FillBoundary();
//         fgeom.FillPeriodicBoundary(ecF[n], true);
    }
}

/*#if 0
void
Electrostatic::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);
    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        if (geom.isPeriodic(dir))
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
        else
        {
            if (phys_bc->lo(dir) == Symmetry)
            {
                mg_bc[2*dir + 0] = MGT_BC_NEU;
            }
            else if (phys_bc->lo(dir) == Outflow)
            {
                mg_bc[2*dir + 0] = MGT_BC_DIR;
            }
            else
            {
                BoxLib::Abort("Unknown lo bc in make_mg_bc");
            }

            if (phys_bc->hi(dir) == Symmetry)
            {
                mg_bc[2*dir + 1] = MGT_BC_NEU;
            }
            else if (phys_bc->hi(dir) == Outflow)
            {
                mg_bc[2*dir + 1] = MGT_BC_DIR;
            }
            else
            {
                BoxLib::Abort("Unknown hi bc in make_mg_bc");
            }
        }
    }
}
#endif
*/

void
Electrostatic::AddParticlesToRhs (int               level,
                            MultiFab&         Rhs,
                            int               ngrow)
{
    // Use the same multifab for all particle types
    MultiFab particle_mf(grids[level], 1, ngrow);

    particle_mf.setVal(0.);
    Accel::thePAPC()->AssignDensitySingleLevel(0, particle_mf, level);
    MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
}

void
Electrostatic::AddParticlesToRhs(int base_level, int finest_level, PArray<MultiFab>& Rhs_particles)
{
    const int num_levels = finest_level - base_level + 1;

    PArray<MultiFab> PartMF;
    //BEGIN MATTHIAS
    //PartMF.resize(num_levels, PArrayManage);
    //PartMF.set(0, new MultiFab(grids[0] , 1, 1));
    //PartMF[0].setVal(0.0);
    //END MATTHIAS
    
    // call: BoxLib ParticleContainer::AssignDensity
//     Accel::thePAPC()->AssignDensity(0, PartMF, base_level, 1, finest_level);
    
    //BEGIN MATTHIAS
    //Accel::thePAPC()->AssignDensitySingleLevel(PartMF[0], 0, 1, 0);
    
    //for (int level = finest_level - 1 - base_level; level >= 0; level--)
    //Electrostatic::average_down(PartMF[level], PartMF[level+1],  parent->refRatio(level+base_level));

    //for (int level = 0; level < num_levels; ++level)
    //	MultiFab::Add(Rhs_particles[base_level + level], PartMF[level], 0, 0, 1, 0);
    
    //std::cerr << Accel::thePAPC()->TotalNumberOfParticles() << " "
	//      << Accel::thePAPC()->sumParticleMass(0) << std::endl;
    //Array<Real> charge;
    //Accel::thePAPC()->GetParticleData(charge, 0, 1);

    //double sum = 0;
    //for (int i = 0; i < charge.size(); ++i)
	//sum += charge[i];
    //std::cerr << sum << std::endl;
    //END MATTHIAS


  // TULIN

    for (int lev = 0; lev < num_levels; lev++)
    {
        if (PartMF[lev].contains_nan())
        {
            *gmsg << "Testing particle density at level " << base_level+lev << endl;
            BoxLib::Abort("...PartMF has NaNs in Electrostatic::actual_multilevel_solve()");
        }
    }
    *gmsg << "finest_level=" << finest_level <<  " base_level=" << base_level << endl;

    for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
    {
        const IntVect ratio = parent->refRatio(lev+base_level);
        *gmsg << "lev=" << lev << " ratio=" << endl; 
        average_down(PartMF[lev], PartMF[lev+1], ratio);
    }

    for (int lev = 0; lev < num_levels; lev++)
    {
        MultiFab::Add(Rhs_particles[lev], PartMF[lev], 0, 0, 1, 0);
    }
}

void
Electrostatic::AddVirtualParticlesToRhs (int               level,
                                         MultiFab&         Rhs,
                                         int               ngrow)
{
    if (level <  parent->finestLevel() && Accel::theVirtPC() != 0)
    {
        // If we have virtual particles, add their density to the single level solve
        MultiFab particle_mf(grids[level], 1, ngrow);
        particle_mf.setVal(0.);

        Accel::theVirtPC()->AssignDensitySingleLevel(0, particle_mf, level, 1, 1);
        MultiFab::Add(Rhs, particle_mf, 0, 0, 1, 0);
    }
}

void
Electrostatic::AddVirtualParticlesToRhs(int finest_level, PArray<MultiFab>& Rhs_particles)
{
    if (finest_level < parent->finestLevel() && Accel::theVirtPC() != 0)
    {
        // Should only need ghost cells for virtual particles if they're near 
        // the simulation boundary and even then only maybe
        MultiFab VirtPartMF(grids[finest_level], 1, 1);
        VirtPartMF.setVal(0.0);

        Accel::theVirtPC()->AssignDensitySingleLevel(0, VirtPartMF, finest_level, 1, 1);
        MultiFab::Add(Rhs_particles[finest_level], VirtPartMF, 0, 0, 1, 0);
    }
}

void
Electrostatic::AddGhostParticlesToRhs (int               level,
                                       MultiFab&         Rhs,
                                       int               ngrow)
{
    if (level > 0 && Accel::theGhostPC() != 0)
    {
        // If we have ghost particles, add their density to the single level solve
        MultiFab ghost_mf(grids[level], 1, ngrow);

        ghost_mf.setVal(0.);
        Accel::theGhostPC()->AssignDensitySingleLevel(0, ghost_mf, level, 1, -1);
        MultiFab::Add(Rhs, ghost_mf, 0, 0, 1, 0);
    }
}

void
Electrostatic::AddGhostParticlesToRhs(int level, PArray<MultiFab>& Rhs_particles)
{
    if (level > 0 && Accel::theGhostPC() != 0)
    {
        // No grow cells for Ghost particles since they're duplicated for
        // all relevant regions.
        MultiFab GhostPartMF(grids[level], 1, 0);
        GhostPartMF.setVal(0.0);

        // Get the Ghost particle mass function. Note that Ghost particles should
        // only affect the coarsest level so we use a single level solve
        Accel::theGhostPC()->AssignDensitySingleLevel(0, GhostPartMF, level, 1, -1);
        MultiFab::Add(Rhs_particles[0], GhostPartMF, 0, 0, 1, 0);
    }
}

void 
Electrostatic::compute_face_fluxes(BoxArray& these_grids, PArray<MultiFab>& grad_phi, 
                                   MultiFab& phi, const Real* dx)
{
    // Initialize the gradients to zero to deal with ghost cell fluxes.
    grad_phi[0].setVal(0.,grad_phi[0].nGrow());
    grad_phi[1].setVal(0.,grad_phi[1].nGrow());
    grad_phi[2].setVal(0.,grad_phi[2].nGrow());

    // Compute the face-based fluxes by simple centered differencing.
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const int  i  = mfi.index();
        const Box& bx = these_grids[i];

        BL_FORT_PROC_CALL(FORT_COMPUTE_EC, fort_compute_ec)
            (bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(phi[i]),
             D_DECL(BL_TO_FORTRAN(grad_phi[0][i]),
                    BL_TO_FORTRAN(grad_phi[1][i]),
                    BL_TO_FORTRAN(grad_phi[2][i])),
             dx);
    }
}
