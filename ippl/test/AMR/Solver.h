#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MacBndry.H>
#include <AMReX_MGT_Solver.H>
#include <mg_cpp_f.h>
#include <AMReX_stencil_types.H>
#include <AMReX_VisMF.H>
#include <AMReX_FMultiGrid.H>

using namespace amrex;

#include <memory>
#include <vector>

// #define USEHYPRE

#ifdef USEHYPRE
#include "HypreABecLap.H"
#endif

/*!
 * @file Solver.h
 * @author Matthias Frey
 *         Ann Almgren
 * @date October 2016, LBNL
 * @details The functions defined in this class
 * are copied from the BoxLib library
 * (i.e. BoxLib/Tutorials/PIC_C). It solves
 * the Poisson equation using a multigrid
 * solver (Gauss-Seidel, V-cycle).
 * @brief V-cycle multi grid solver
 */

/// Do a MultiGrid solve
class Solver {

public:
    typedef Array<std::unique_ptr<MultiFab> > container_t;
    typedef Array<MultiFab*> container_pt;

    /*!
     * Prepares the solver and calls the solve_with_f90 function.
     * @param rhs is the density at each level (cell-centered)
     * @param phi is the potential at each level (cell-centered)
     * @param grad_phi is the electric field at each level (cell-centered)
     * @param geom is the geometry at each level
     * @param base_level from which the solve starts
     * @param finest_level up to which solver goes
     * @param offset is zero in case of Dirichlet boundary conditions.
     * @param timing of solver parts
     * @param doGradient compute the gradient (true) or not (false)
     */
    void solve_for_accel(const container_t& rhs,
                         const container_t& phi,
                         const container_t& grad_phi,
                         const Array<Geometry>& geom,
                         int base_level,
                         int finest_level,
                         Real offset,
                         bool timing=true,
                         bool doGradient=true);
    /*!
     * Actual solve.
     * @param rhs is the density at each level (cell-centered)
     * @param phi is the potential at each level (cell-centered)
     * @param grad_phi_edge is the electric field at each level (at cell-faces)
     * @param geom is the geometry at each level.
     * @param base_level from which the solve starts
     * @param finest_level up to which solver goes
     * @param tol is \f$ 10^{-10}\f$ (specified in solve_for_accel)
     * @param abs_tol is \f$ 10^{-14}\f$ (specified in solve_for_accel)
     * @param timing of solver parts
     * @param doGradient  compute the gradient (true) or not (false)
     */
    void solve_with_f90(const container_pt& rhs,
                        const container_pt& phi, const Array<container_pt>& grad_phi_edge, 
                        const Array<Geometry>& geom, int base_level, int finest_level, Real tol, Real abs_tol,
                        bool timing, bool doGradient);
    
#ifdef USEHYPRE
    void solve_with_hypre(MultiFab& soln, MultiFab& rhs, const BoxArray& bs,
                          const Geometry& geom);

private:
    void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);

#endif
};


#endif