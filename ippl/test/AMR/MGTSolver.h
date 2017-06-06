#ifndef MGTSOLVER_H
#define MGTSOLVER_H

#include <AMReX_MGT_Solver.H>
#include <memory>

using namespace amrex;

class MGTSolver {
    
public:
    typedef amrex::Array<std::unique_ptr<amrex::MultiFab> > container_t;
    
    void solve(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const amrex::Array<amrex::Geometry>& geom);
    
};

#endif
