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


#include <memory>
#include <vector>

#ifdef UNIQUE_PTR
    #include <PArray.H>
#endif

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
#ifdef UNIQUE_PTR
    typedef Array<std::unique_ptr<MultiFab> > container_t;
#else
    typedef PArray<MultiFab> container_t;
#endif
    typedef Array<MultiFab*> container_pt;      // instead of PArray<MultiFab>

    /*!
     * Prepares the solver and calls the solve_with_f90 function.
     * @param rhs is the density at each level (cell-centered)
     * @param phi is the potential at each level (cell-centered)
     * @param grad_phi is the electric field at each level (cell-centered)
     * @param geom is the geometry at each level
     * @param base_level from which the solve starts
     * @param finest_level up to which solver goes
     * @param offset is zero in case of Dirichlet boundary conditions.
     */
    void solve_for_accel(container_t& rhs,
                         container_t& phi,
                         container_t& grad_phi,
                         const Array<Geometry>& geom,
                         int base_level,
                         int finest_level,
                         Real offset);


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
     */
    void solve_with_f90(container_t& rhs,
                        container_t& phi, Array< container_t >& grad_phi_edge,
                        const Array<Geometry>& geom, int base_level, int finest_level, Real tol, Real abs_tol);

#ifdef USEHYPRE
    void solve_with_hypre(MultiFab& soln, MultiFab& rhs, const BoxArray& bs,
                          const Geometry& geom);

private:
    void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);

#endif
};


#endif