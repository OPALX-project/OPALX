#ifndef MGTSOLVER_H
#define MGTSOLVER_H

#include <AMReX_MGT_Solver.H>
#include <memory>

/// Based on https://github.com/AMReX-Codes/Nyx.git Gravity.cpp
class MGTSolver {
    
public:
    typedef amrex::Array<std::unique_ptr<amrex::MultiFab> > container_t;
    
    void solve(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const amrex::Array<amrex::Geometry>& geom);
    
    
//     void set_dirichlet_bcs (int level, MultiFab* phi);
    
};

#endif
