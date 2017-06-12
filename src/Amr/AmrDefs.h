#ifndef AMR_DEFS_H
#define AMR_DEFS_H

#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <memory>

namespace amr {
    typedef amrex::MultiFab                                 AmrField_t;
    typedef amrex::DistributionMapping                      AmrProcMap_t;
    typedef amrex::Geometry                                 AmrGeometry_t;
    typedef amrex::BoxArray                                 AmrGrid_t;
    typedef amrex::Array< std::unique_ptr<AmrField_t>  >    AmrFieldContainer_t;
    typedef amrex::Array< AmrGeometry_t >                   AmrGeomContainer_t;
    typedef amrex::Array< AmrGrid_t >                       AmrGridContainer_t;
    typedef amrex::Array< AmrProcMap_t >                    AmrProcMapContainer_t;
    typedef amrex::RealBox                                  AmrDomain_t;
    typedef amrex::Array<int>                               AmrIntArray_t;
    typedef amrex::IntVect                                  AmrIntVect_t;
    typedef amrex::Array< AmrIntVect_t >                    AmrIntVectContainer_t;
    typedef amrex::Box                                      AmrBox_t;
    typedef amrex::Real                                     AmrReal_t;
};

#endif
