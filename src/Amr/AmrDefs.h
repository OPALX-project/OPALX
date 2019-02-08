#ifndef AMR_DEFS_H
#define AMR_DEFS_H

#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <memory>

/// Some AMR types used a lot
namespace amr {
    typedef amrex::MultiFab                                 AmrField_t;
    typedef std::array< std::unique_ptr<AmrField_t>,
                        AMREX_SPACEDIM
                       >                                    AmrVectorField_t;
    typedef amrex::DistributionMapping                      AmrProcMap_t;
    typedef amrex::Geometry                                 AmrGeometry_t;
    typedef amrex::BoxArray                                 AmrGrid_t;
    typedef amrex::Vector< std::unique_ptr<AmrField_t> >    AmrScalarFieldContainer_t;
    typedef amrex::Vector< AmrVectorField_t >               AmrVectorFieldContainer_t;
    typedef amrex::Vector< AmrGeometry_t >                  AmrGeomContainer_t;
    typedef amrex::Vector< AmrGrid_t >                      AmrGridContainer_t;
    typedef amrex::Vector< AmrProcMap_t >                   AmrProcMapContainer_t;
    typedef amrex::RealBox                                  AmrDomain_t;
    typedef amrex::Vector<int>                              AmrIntArray_t;
    typedef amrex::IntVect                                  AmrIntVect_t;
    typedef amrex::Vector< AmrIntVect_t >                   AmrIntVectContainer_t;
    typedef amrex::Box                                      AmrBox_t;
    typedef amrex::Real                                     AmrReal_t;
};

#endif
