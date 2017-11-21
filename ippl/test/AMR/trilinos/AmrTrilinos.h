#ifndef AMR_TRILINOS_H
#define AMR_TRILINOS_H

#include <AMReX.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>


#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <map>
#include <memory>

#include "EllipticDomain.h"

#include "electrostatic_pic_F.H"


#include <AMReX_BCUtil.H>
#include <AMReX_MacBndry.H>

#include <AMReX_FluxRegister.H>

using namespace amrex;

class AmrTrilinos {
    
public:
    typedef std::unique_ptr<MultiFab>   AmrField_t;
    typedef Array<AmrField_t>           container_t;
    
public:
    
    AmrTrilinos();
    
    void solveWithGeometry(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const Array<Geometry>& geom,
               int lbase,
               int lfine);
    
    void solve(const container_t& rho,
               container_t& phi,
               container_t& efield,
               const Array<Geometry>& geom,
               int lbase,
               int lfine,
               double scale);
    
    
private:
    
    void solveLevel_m(AmrField_t& phi,
                          const AmrField_t& rho,
                          const Geometry& geom, int lev);
    
    // single level solve
    void linSysSolve_m(AmrField_t& phi, const Geometry& geom);
    
//     void solveGrid_m(
    
    // create map + right-hand side
    void build_m(const AmrField_t& rho, const AmrField_t& phi, const Geometry& geom, int lev);
    
    void fillMatrix_m(const Geometry& geom, const AmrField_t& phi, int lev);
    
    int serialize_m(const IntVect& iv) const;
    
//     IntVect deserialize_m(int idx) const;
    
    bool isInside_m(const IntVect& iv) const;
    
//     bool isBoundary_m(const IntVect& iv) const;
    
//     void findBoundaryBoxes_m(const BoxArray& ba);
    
    void copyBack_m(AmrField_t& phi,
                    const Teuchos::RCP<Epetra_MultiVector>& sol,
                    const Geometry& geom);
    
    
//     void grid_m(AmrField_t& rhs, const container_t& rho,
//                 const Geometry& geom, int lev);
    
    void interpFromCoarseLevel_m(container_t& phi,
                                 const Array<Geometry>& geom,
                                 int lev);
    
    // efield = -grad phi
    void getGradient(AmrField_t& phi,
                     AmrField_t& efield,
                     const Geometry& geom, int lev);
    
    void getGradient_1(AmrField_t& phi, AmrField_t& efield,
                       const Geometry& geom, int lev);
    
    
    void buildLevelMask_m(const BoxArray& grids,
                          const DistributionMapping& dmap,
                          Geometry& geom,
                          int lev)
    {
        
        masks_m[lev].reset(new FabArray<BaseFab<int> >(grids, dmap, 1, 1));
        
        // covered   : ghost cells covered by valid cells of this FabArray
        //             (including periodically shifted valid cells)
        // notcovered: ghost cells not covered by valid cells
        //             (including ghost cells outside periodic boundaries)
        // physbnd   : boundary cells outside the domain (excluding periodic boundaries)
        // interior  : interior cells (i.e., valid cells)
        int covered = -1;
        int interior = 0;
        int notcovered = 1; // -> boundary
        int physbnd = 2;
        
        masks_m[lev]->BuildMask(geom.Domain(), geom.periodicity(),
                                covered, notcovered, physbnd, interior);
    }
    
    
//     void buildGeometry_m(const AmrField_t& rho, const AmrField_t& phi,
//                          const Geometry& geom, int lev, EllipticDomain& bp);
    
//     void solveLevelGeometry_m(AmrField_t& phi,
//                               const AmrField_t& rho,
//                               const Geometry& geom, int lev, EllipticDomain& bp);
    
//     void fillMatrixGeometry_m(const Geometry& geom, const AmrField_t& phi, int lev, EllipticDomain& bp);
    
    void setBoundaryValue_m(AmrField_t& phi, const AmrField_t& crse_phi = 0,
                            IntVect crse_ratio = IntVect::TheZeroVector())
    {
        // The values of phys_bc & ref_ratio do not matter
        // because we are not going to use those parts of MacBndry.
        IntVect ref_ratio = IntVect::TheZeroVector();
        Array<int> lo_bc(AMREX_SPACEDIM, 0);
        Array<int> hi_bc(AMREX_SPACEDIM, 0);
        BCRec phys_bc(lo_bc.dataPtr(), hi_bc.dataPtr());
    
//         if (crse_phi == 0 && phi == 0) 
//         {
//             bndry_m.setHomogValues(phys_bc, ref_ratio);
//         }
        if (crse_phi == 0 && phi != 0)
        {
            bndry_m->setBndryValues(*phi, 0, 0, phi->nComp(), phys_bc); 
        }
        else if (crse_phi != 0 && phi != 0) 
        {
            BL_ASSERT(crse_ratio != IntVect::TheZeroVector());
    
            const int ncomp      = phi->nComp();
            const int in_rad     = 0;
            const int out_rad    = 1;
            const int extent_rad = 2;
    
            BoxArray crse_boxes(phi->boxArray());
            crse_boxes.coarsen(crse_ratio);
    
            BndryRegister crse_br(crse_boxes, phi->DistributionMap(),
                                in_rad, out_rad, extent_rad, ncomp);
            crse_br.copyFrom(*crse_phi, crse_phi->nGrow(), 0, 0, ncomp);
            
            bndry_m->setBndryValues(crse_br, 0, *phi, 0, 0, ncomp, crse_ratio, phys_bc, 3);
            
//             print_m(phi);
            
        }
        else
        {
            amrex::Abort("FMultiGrid::Boundary::build_bndry: How did we get here?");
        }
    }
    
    void print_m(const AmrField_t& phi) {
        for (MFIter mfi(*phi, false); mfi.isValid(); ++mfi) {
            const Box&          bx  = mfi.validbox();
            const FArrayBox&    fab = (*phi)[mfi];
            
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();
        
            for (int i = lo[0]-1; i <= hi[0]+1; ++i) {
                for (int j = lo[1]-1; j <= hi[1]+1; ++j) {
                    for (int k = lo[2]-1; k <= hi[2]+1; ++k) {
                        IntVect iv(i, j, k);
                        std::cout << iv << " " << fab(iv) << std::endl;// std::cin.get();
                    }
                }
            }
        }
    }
    
    void boundaryprint_m() {
        
        bc_m.clear();
        for (OrientationIter oitr; oitr; ++oitr) {
            const FabSet& fs = bndry_m->bndryValues(oitr());
            
            for (FabSetIter fsi(fs); fsi.isValid(); ++fsi) {
                const Box&          bx  = fsi.validbox();
                const FArrayBox&    fab = fs[fsi];
                const int* lo = bx.loVect();
                const int* hi = bx.hiVect();
            
                for (int ii = lo[0]; ii <= hi[0]; ++ii) {
                    for (int j = lo[1]; j <= hi[1]; ++j) {
                        for (int k = lo[2]; k <= hi[2]; ++k) {
                            IntVect iv(ii, j, k);
                            bc_m.insert(iv);
//                             std::cout << iv << " " << fab(iv) << std::endl; std::cin.get();
                        }
                    }
                }
            }
        }
    }
    
    void initGhosts_m() {
        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face  = fi();
            FabSet& fs = (*bndry_m.get())[face];
            
            for (FabSetIter fsi(fs); fsi.isValid(); ++fsi)
            {
                fs[fsi].setVal(0);
            }
        }
        
    }
    
    void fillDomainBoundary_m(MultiFab& phi, const Geometry& geom) {
        
        Array<BCRec> bc(phi.nComp());
        for (int n = 0; n < phi.nComp(); ++n) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if ( Geometry::isPeriodic(idim) ) {
                    bc[n].setLo(idim, BCType::int_dir); // interior
                    bc[n].setHi(idim, BCType::int_dir);
                } else {
                    bc[n].setLo(idim, BCType::foextrap); // first oder extrapolation from last cell in interior
                    bc[n].setHi(idim, BCType::foextrap);
                }
            }
        }
        
        // fill physical boundary + ghosts at physical boundary
        amrex::FillDomainBoundary(phi, geom, bc);
    }
    
    void computeFlux_m(const AmrField_t& phi, const Geometry& geom, int lev) {
        
        fluxes_m.resize(lev + 1);
        
        fluxes_m[lev].resize(AMREX_SPACEDIM);
        
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            BoxArray ba = phi->boxArray();
            const DistributionMapping& dmap = phi->DistributionMap();
            fluxes_m[lev][n].reset(new MultiFab(ba.surroundingNodes(n), dmap, 1, 1));
            fluxes_m[lev][n]->setVal(0.0, 1);
        }
        
        const double* dx = geom.CellSize();
        
        for (MFIter mfi(*masks_m[lev], false); mfi.isValid(); ++mfi) {
            const Box&          bx   = mfi.validbox();
            const FArrayBox&    pfab = (*phi)[mfi];
            FArrayBox&          xface = (*fluxes_m[lev][0])[mfi];
            FArrayBox&          yface = (*fluxes_m[lev][1])[mfi];
            FArrayBox&          zface = (*fluxes_m[lev][2])[mfi];
            
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();
            
            for (int i = lo[0]; i <= hi[0]+1; ++i) {
                for (int j = lo[1]; j <= hi[1]+1; ++j) {
                    for (int k = lo[2]; k <= hi[2]+1; ++k) {
                        
                        IntVect iv(i, j, k);
                        
                        
                        // x direction
                        IntVect left(i-1, j, k);
                        xface(iv) = (pfab(left) - pfab(iv)) / dx[0];
                        
                        // y direction
                        IntVect down(i, j-1, k);
                        yface(iv) = (pfab(down) - pfab(iv)) / dx[1];
                        
                        // z direction
                        IntVect front(i, j, k-1);
                        zface(iv) = (pfab(front) - pfab(iv)) / dx[2];
                    }
                }
            }
        }
    }
    
    void correctFlux_m(container_t& phi, const Array<Geometry>& geom, int lev) {
        
        const BoxArray& grids = phi[lev+1]->boxArray();
        const DistributionMapping& dmap = phi[lev+1]->DistributionMap();
        IntVect crse_ratio(2, 2, 2);
        std::unique_ptr<FluxRegister> fine = std::unique_ptr<FluxRegister>(new FluxRegister(grids,
                                                                                            dmap,
                                                                                            crse_ratio,
                                                                                            lev+1,
                                                                                            1));
        computeFlux_m(phi[lev+1], geom[lev+1], lev+1);
        
        
        for (int i = 0; i < AMREX_SPACEDIM ; i++) {
            Array<MultiFab*> p = amrex::GetArrOfPtrs(fluxes_m[lev+1]);
            fine->FineAdd((*p[i]), i, 0, 0, 1, 1.);
        }
        
        
        fine->Reflux(*phi[lev], 1.0, 0, 0, 1, geom[lev]);
        
    }
    
    
    void fixRHSForSolve(Array<std::unique_ptr<MultiFab> >& rhs,
                        const Array<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
                        const Array<Geometry>& geom, const IntVect& ratio)
    {
        int num_levels = rhs.size();
        for (int lev = 0; lev < num_levels; ++lev) {
            MultiFab& fine_rhs = *rhs[lev];
            const FabArray<BaseFab<int> >& mask = *masks[lev];        
            const BoxArray& fine_ba = fine_rhs.boxArray();
            const DistributionMapping& fine_dm = fine_rhs.DistributionMap();
            MultiFab fine_bndry_data(fine_ba, fine_dm, 1, 1);
            zeroOutBoundary(fine_rhs, fine_bndry_data, mask);
        }
    }
    
    void zeroOutBoundary(MultiFab& input_data,
                     MultiFab& bndry_data,
                     const FabArray<BaseFab<int> >& mask)
    {
        bndry_data.setVal(0.0, 1);
        for (MFIter mfi(input_data); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            zero_out_bndry(bx.loVect(), bx.hiVect(),
                        input_data[mfi].dataPtr(),
                        bndry_data[mfi].dataPtr(),
                        mask[mfi].dataPtr());
        }
        bndry_data.FillBoundary();
    }
    
    
private:
    
    AmrField_t fab_set_m;
    
    std::unique_ptr<MacBndry> bndry_m;
    
    Array<std::unique_ptr<FabArray<BaseFab<int> > > > masks_m;
    
    Array<container_t> fluxes_m;
    
    std::set<IntVect> bc_m;
    
    Epetra_MpiComm epetra_comm_m;
    
    Teuchos::RCP<Epetra_Map> map_mp;
    
    int nLevel_m;
    
    // info for a level
    IntVect nPoints_m;
    
    // (x, y, z) --> idx
    std::map<int, IntVect> indexMap_m;
    
    
    Teuchos::RCP<Epetra_CrsMatrix> A_mp;
    Teuchos::RCP<Epetra_Vector> rho_mp;
    Teuchos::RCP<Epetra_Vector> x_mp;
    
};

#endif
