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


#endif
