/*!
 * @file helper_functions.h
 * @author Matthias Frey
 * @date 23. Dec. 2016
 * @brief Define several functions to avoid
 *        code duplication.
 * @version 1.0
 */

#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <fstream>

using amrex::Array;
using amrex::MultiFab;
using amrex::DistributionMapping;
using amrex::BoxArray;
using amrex::Geometry;
using amrex::MFIter;
using amrex::Box;
using amrex::FArrayBox;
using amrex::IntVect;
using amrex::RealBox;
using amrex::Real;

typedef Array<std::unique_ptr<MultiFab> > container_t;

/*!
 * @param scalfield to write
 * @param filename to be written to
 * @param step of the file
 */
inline void writeScalarField(const container_t& scalfield, std::string filename, int step)
{
    for (MFIter mfi(*scalfield[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& lhs = (*scalfield[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                std::string outfile = filename + "_step_" + std::to_string(step) + ".grid";
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
                    for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
                        for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " "
                                << lhs(ivec, 0)  << std::endl;
                        }
                    }
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
}

/*!
 * Write the grid data along the horizontal direction
 * (centered in y and z) to a file. Just single level supported.
 * @param scalfield is the scalar field
 * @param dx is the mesh spacing
 * @param dlow is the lower boundary of the domain (not grid point)
 * @param filename where to write
 */
inline void writeScalarField(const container_t& scalfield,
                             double dx,
                             double dlow,
                             std::string filename)
{
    for (MFIter mfi(*scalfield[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& lhs = (*scalfield[0])[mfi];
        
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                std::string outfile = filename + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0)  << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
}

/*!
 * Write the grid data along the horizontal direction
 * (centered in y and z) to a file. Just single level supported.
 * @param vecfield is the scalar field
 * @param dx is the mesh spacing
 * @param dlow is the lower boundary of the domain (not grid point)
 */
inline void writeVectorField(const container_t& vecfield,
                             double dx,
                             double dlow)
{
    for (MFIter mfi(*vecfield[0 /*level*/]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const FArrayBox& lhs = (*vecfield[0])[mfi];
        /* Write the potential of all levels to file in the format: x, y, z, ex, ey, ez
         */
        for (int proc = 0; proc < amrex::ParallelDescriptor::NProcs(); ++proc) {
            if ( proc == amrex::ParallelDescriptor::MyProc() ) {
                std::string outfile = "amr-efield_scalar-level-" + std::to_string(0);
                std::ofstream out;
                
                if ( proc == 0 )
                    out.open(outfile);
                else
                    out.open(outfile, std::ios_base::app);
                
                int j = 0.5 * (bx.hiVect()[1] - bx.loVect()[1]);
                int k = 0.5 * (bx.hiVect()[2] - bx.loVect()[2]);
                
                for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                     for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                         for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k) {
                            IntVect ivec(i, j, k);
                            // add one in order to have same convention as PartBunch::computeSelfField()
                            out << i + 1 << " " << j + 1 << " " << k + 1 << " " << dlow + dx * i << " "
                                << lhs(ivec, 0) << " " << lhs(ivec, 1) << " " << lhs(ivec, 2) << std::endl;
//                         }
//                     }
                }
                out.close();
            }
            amrex::ParallelDescriptor::Barrier();
        }
    }
}

/*!
 * Compute the electric field energy in \f$ [CV\epsilon_0] \f$ units.
 * The values at finer grids is averaged down to the coarsest grid.
 * @param efield specifies the electric field
 * @param rr are the refinement ratios.
 */
inline double totalFieldEnergy(container_t& efield, const Array<int>& rr) {
    
    for (int lev = efield.size() - 2; lev >= 0; lev--)
        amrex::average_down(*(efield[lev+1].get()), *(efield[lev].get()), 0, 3, rr[lev]);
    
    double fieldEnergy = 0.0;
//     long volume = 0;
//     for (int l = 0; l < efield.size(); ++l) {
        long nPoints = 0;

        int l = 0;
        fieldEnergy += MultiFab::Dot(*efield[l], 0, *efield[l], 0, 3, 0);
        for (unsigned int b = 0; b < efield[l]->boxArray().size(); ++b) {
//             volume += efield[l]->box(b).volume(); // cell-centered --> volume() == numPts
            nPoints += efield[l]->box(b).numPts();
        }
        fieldEnergy /= double(nPoints);
//     }
    
    return 0.5 * /*Physics::epsilon_0 * */fieldEnergy;
}

/*!
 * @param domain is the physical domain
 * @param ba are the boxes at the coarsest level
 * @param dmap is the mapping to the processors
 * @param geom are the geometries at the different levels
 * @param rr are the refinement ratios among levels
 * @param nr are the number of grid points in each dimension
 * @param nLevels is 1 for single level
 * @param maxBoxSize is the maximum size per grid
 * @param lower is the physical lower bound of the domain
 * @param upper is the physical upper bound of the domain
 */
inline void init(RealBox& domain,
                 Array<BoxArray>& ba,
                 Array<DistributionMapping>& dmap,
                 Array<Geometry>& geom,
                 Array<int>& rr,
                 const Vektor<size_t, 3>& nr,
                 int nLevels,
                 size_t maxBoxSize,
                 const std::array<double, BL_SPACEDIM>& lower,
                 const std::array<double, BL_SPACEDIM>& upper)
{
    /*
     * nLevels is the number of levels allowed, i.e if nLevels = 1
     * we just run single-level
     */
    
    /*
     * set up the geometry
     */
    IntVect low(0, 0, 0);
    IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
    Box bx(low, high);
    
    // box
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, lower[i]); // m
        domain.setHi(i, upper[i]); // m
    }
    
    // Dirichlet boundary conditions in all directions
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    geom.resize(nLevels);
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    // Container for boxes at all levels
    ba.resize(nLevels);
    
    // box at level 0
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);
    
    /*
     * distribution mapping
     */
    dmap.resize(nLevels);
    dmap[0].define(ba[0], amrex::ParallelDescriptor::NProcs() /*nprocs*/);
    
    // refinement ratio
    rr.resize(nLevels - 1);
    for (std::size_t lev = 0; lev < rr.size(); ++lev)
        rr[lev] = 2;
    
    // geometries of refined levels
    for (int lev = 1; lev < nLevels; ++lev) {
        geom[lev].define(amrex::refine(geom[lev - 1].Domain(),
                                        rr[lev - 1]),
                         &domain, 0, bc);
    }
}

/*!
 * @param domain is the physical domain
 * @param nr are the number of grid points in each dimension
 * @param lower is the physical lower bound of the domain
 * @param upper is the physical upper bound of the domain
 */
inline void init(RealBox& domain,
                 const Vektor<size_t, 3>& nr,
                 const std::array<double, BL_SPACEDIM>& lower,
                 const std::array<double, BL_SPACEDIM>& upper)
{
    /*
     * nLevels is the number of levels allowed, i.e if nLevels = 1
     * we just run single-level
     */
    
    /*
     * set up the geometry
     */
//     IntVect low(0, 0, 0);
//     IntVect high(nr[0] - 1, nr[1] - 1, nr[2] - 1);    
//     Box bx(low, high);
    
    // box
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, lower[i]); // m
        domain.setHi(i, upper[i]); // m
    }
    
//     // Dirichlet boundary conditions in all directions
//     int bc[BL_SPACEDIM] = {0, 0, 0};
}

/*!
 * Allocate memory for the solver and initialize
 * the grid data to zero.
 * @param rhs is the right-hand side of the equation
 * @param phi is th scalar potential
 * @param is the electric field
 * @param ba is the box array per level
 * @param level for which we initialize the data
 */
inline void initGridData(container_t& rhs,
                         container_t& phi,
                         container_t& grad_phi,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         int level)
{
    //                                                               # component # ghost cells                                                                                                                                          
    rhs[level] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm,      1          , 0));
    phi[level] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm,      1          , 1));
    grad_phi[level] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm, BL_SPACEDIM, 1));

    rhs[level]->setVal(0.0);
    phi[level]->setVal(0.0, 1);
    grad_phi[level]->setVal(0.0, 1);
}

/*!
 * Allocate memory for the solver and initialize
 * the grid data to zero.
 * @param rhs is the right-hand side of the equation
 * @param is the electric field
 * @param ba is the box array per level
 * @param level for which we initialize the data
 */
inline void initGridData(container_t& rhs,
                         container_t& grad_phi,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         int level)
{
    //                                                               # component # ghost cells                                                                                                                                          
    rhs[level] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm,      1          , 0));
    grad_phi[level] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm, BL_SPACEDIM, 1));

    rhs[level]->setVal(0.0);
    grad_phi[level]->setVal(0.0);
}

/*!
 * Compute the total charge from the grid data
 * @param rhs is the charge density on the grid
 * @param finest_level of AMR
 * @param geom are the geometries of each level
 * @param scale with the cell volume
 */
inline double totalCharge(const container_t& rhs,
                          int finest_level,
                          const Array<Geometry>& geom,
                          bool scale = true)
{
    
//     Real totCharge = 0.0;
//     for (int i = 0; i <= finest_level; ++i) {
//         Real vol = (*(geom[i].CellSize()) * *(geom[i].CellSize()) * *(geom[i].CellSize()) );
//         Real sum = (scale) ? rhs[i]->sum(0) * vol : rhs[i]->sum(0);
//         totCharge += sum;
//     }
    
    //AssignDensityFort averages down
    Real totCharge = 0.0;
    Real vol = (*(geom[0].CellSize()) * *(geom[0].CellSize()) * *(geom[0].CellSize()) );
    totCharge = (scale) ? rhs[0]->sum(0) * vol : rhs[0]->sum(0);
    return totCharge;
}

#endif
