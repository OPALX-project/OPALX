#include "AmrMultiGrid.h"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Array.H>

#include <AMReX_MacBndry.H>


#define WRITE 0

#if WRITE
    #include <fstream>
#endif

AmrMultiGrid::AmrMultiGrid(Boundary bc,
                           Interpolater interp,
                           Interpolater interface,
                           LinSolver solver)
    : nIter_m(0),
      nsmooth_m(12),
      lbase_m(0),
      lfine_m(0),
      bc_m(bc)
{
    comm_mp = Teuchos::rcp( new comm_t( Teuchos::opaqueWrapper(Ippl::getComm()) ) );
    node_mp = KokkosClassic::DefaultNode::getDefaultNode();
    
    buildTimer_m        = IpplTimings::getTimer("build");
    restrictTimer_m     = IpplTimings::getTimer("restrict");
    smoothTimer_m       = IpplTimings::getTimer("smooth");
    interpTimer_m       = IpplTimings::getTimer("prolongate");
    residnofineTimer_m  = IpplTimings::getTimer("resid-no-fine");
    
    this->initInterpolater_m(interp);
    
    // interpolater for crse-fine-interface
    this->initCrseFineInterp_m(interface);
    
    // base level solver
    this->initBaseSolver_m(solver);
}

void AmrMultiGrid::solve(const amrex::Array<AmrField_u>& rho,
                         amrex::Array<AmrField_u>& phi,
                         amrex::Array<AmrField_u>& efield,
                         const amrex::Array<AmrGeometry_t>& geom,
                         int lbase, int lfine, bool previous)
{
    nIter_m = 0;
    lbase_m = lbase;
    lfine_m = lfine;
    int nLevel = lfine - lbase + 1;
    
    AmrIntVect_t rr = AmrIntVect_t(D_DECL(2, 2, 2));
    
    // initialize all levels
    mglevel_m.resize(nLevel);
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        switch ( bc_m ) {
            
            case Boundary::DIRICHLET:
            {
                mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev]->boxArray(),
                                                             rho[ilev]->DistributionMap(),
                                                             geom[ilev],
                                                             rr,
                                                             new AmrDirichletBoundary<AmrMultiGridLevel_t>(), 
                                                             comm_mp,
                                                             node_mp));
                break;
            }
            case Boundary::OPEN:
            {
                mglevel_m[lev].reset(new AmrMultiGridLevel_t(rho[ilev]->boxArray(),
                                                             rho[ilev]->DistributionMap(),
                                                             geom[ilev],
                                                             rr,
                                                             new AmrOpenBoundary<AmrMultiGridLevel_t>(),
                                                             comm_mp,
                                                             node_mp));
                break;
            }
            default:
                throw std::runtime_error("Error: This type of boundary is not supported");
        }
    }
    
    if ( !previous ) {
        // reset
        for (int lev = 0; lev < nLevel; ++lev) {
            int ilev = lbase + lev;
            phi[ilev]->setVal(0.0, phi[ilev]->nGrow());
        }
    }
    
    // build all necessary matrices and vectors
    IpplTimings::startTimer(buildTimer_m);
    for (int lev = 0; lev < nLevel; ++lev) {
        
        this->buildRestrictionMatrix_m(lev);
        
        this->buildInterpolationMatrix_m(lev);
        
        this->buildCrseBoundaryMatrix_m(lev);
        
        this->buildFineBoundaryMatrix_m(lev);
        
        this->buildSmootherMatrix_m(lev);
        
        this->buildNoFinePoissonMatrix_m(lev);
        
        this->buildCompositePoissonMatrix_m(lev);
        
        int ilev = lbase + lev;
        
        this->buildDensityVector_m(*rho[ilev], lev);
        
        this->buildPotentialVector_m(*phi[ilev], lev);
        
        this->buildGradientMatrix_m(lev);
        
        nSweeps_m.push_back(nsmooth_m >> lev);
    }
    IpplTimings::stopTimer(buildTimer_m);
    
    mglevel_m[lfine_m]->error_p->putScalar(0.0);
    
    double err = 1.0e7;
    double eps = 1.0e-8;
    
    // initial error
    double residualsum = 0.0;
    double rhosum = 0.0;
    
    for (int lev = 0; lev < nLevel; ++lev) {
        this->residual_m(mglevel_m[lev]->residual_p,
                         mglevel_m[lev]->rho_p,
                         mglevel_m[lev]->phi_p, lev);
        
        double tmp = mglevel_m[lev]->residual_p->norm2();
        residualsum += tmp;
        
        tmp = mglevel_m[lev]->rho_p->norm2();
        rhosum += tmp;
    }
    
#if WRITE
    std::ofstream out("residual.dat");
#endif
    
    
    while ( residualsum > eps * rhosum) {
        
        std::cout << residualsum << " " << eps * rhosum << std::endl;
        
        relax_m(lfine);
        
//         // update residual
//         for (int lev = 0; lev < nLevel; ++lev) {
//             
//             this->residual_m(mglevel_m[lev]->residual_p,
//                              mglevel_m[lev]->rho_p,
//                              mglevel_m[lev]->phi_p, lev);
//         }
        
        residualsum = l2error_m();
        
#if WRITE
        if ( Ippl::myNode() == 0 ) {
            out << residualsum << std::endl;
            std::cerr << residualsum << std::endl;
        }
#endif
        
        ++nIter_m;
    }
    
#if WRITE
    out.close();
#endif
    
    // evaluate the electric field
    Teuchos::RCP<vector_t> efield_p = Teuchos::null;
    for (int lev = nLevel - 1; lev > -1; --lev) {
        int ilev = lbase + lev;
        
        efield_p = Teuchos::rcp( new vector_t(mglevel_m[lev]->map_p, false) );
        
        averageDown_m(lev);
        
        for (int d = 0; d < BL_SPACEDIM; ++d) {
            mglevel_m[lev]->G_p[d]->apply(*mglevel_m[lev]->phi_p, *efield_p);
            this->trilinos2amrex_m(*efield[ilev], d, efield_p);
        }
    }
    
    
    // copy solution back
    for (int lev = 0; lev < nLevel; ++lev) {
        int ilev = lbase + lev;
        
        this->trilinos2amrex_m(*phi[ilev], 0, mglevel_m[lev]->phi_p);
    }
}


void AmrMultiGrid::residual_m(Teuchos::RCP<vector_t>& r,
                              const Teuchos::RCP<vector_t>& b,
                              const Teuchos::RCP<vector_t>& x,
                              int level)
{
    /*
     * r = b - A*x
     */
    if ( level < lfine_m ) {
        vector_t fine2crse(mglevel_m[level]->Anf_p->getDomainMap());
        
        // get boundary for 
        if ( mglevel_m[level]->Bfine_p != Teuchos::null ) {
            mglevel_m[level]->Bfine_p->apply(*mglevel_m[level+1]->phi_p, fine2crse);
        }
        
        vector_t tmp2(mglevel_m[level]->Awf_p->getDomainMap());
        mglevel_m[level]->Awf_p->apply(*mglevel_m[level]->phi_p, tmp2);
        
        vector_t crse2fine(mglevel_m[level]->Anf_p->getDomainMap());
        
        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p, crse2fine);
        }
        
        tmp2.update(1.0, fine2crse, 1.0, crse2fine, 1.0);
        
        Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    
        mglevel_m[level]->UnCovered_p->apply(*mglevel_m[level]->rho_p, *tmp3);
    
        // ONLY subtract coarse rho
        mglevel_m[level]->residual_p->update(1.0, *tmp3, -1.0, tmp2, 0.0);
        
    } else {
        this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                 mglevel_m[level]->phi_p,
                                 mglevel_m[level-1]->phi_p,
                                 mglevel_m[level]->rho_p, level);
        
    }
}


void AmrMultiGrid::relax_m(int level) {
    
    if ( level == lfine_m ) {
        
        if ( level == lbase_m ) {
            Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
            this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     tmp,
                                     mglevel_m[level]->rho_p, level);
        } else {
            //TODO residual is already computed at beginning
            
            this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     mglevel_m[level-1]->phi_p,
                                     mglevel_m[level]->rho_p, level);
        }
    }
    
    if ( level > 0 ) {
        //TODO
        
        // phi^(l, save) = phi^(l)
        vector_t phi_save = *mglevel_m[level]->phi_p;
        
        mglevel_m[level-1]->error_p->putScalar(0.0);
        
        // smoothing
        IpplTimings::startTimer(smoothTimer_m);
        for (std::size_t iii = 0; iii < nSweeps_m[level]; ++iii)
            this->gsrb_level_m(mglevel_m[level]->error_p,
                               mglevel_m[level]->residual_p, level);
        IpplTimings::stopTimer(smoothTimer_m);
        
        // phi = phi + e
        mglevel_m[level]->phi_p->update(1.0, *mglevel_m[level]->error_p, 1.0);
        
        /*
         * restrict
         */
        IpplTimings::startTimer(restrictTimer_m);
        this->restrict_m(level-1);
        IpplTimings::stopTimer(restrictTimer_m);
        
        this->relax_m(level - 1);
        
        /*
         * prolongate / interpolate
         */
        
        IpplTimings::startTimer(interpTimer_m);
        // interpolate error from l-1 to l
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        mglevel_m[level-1]->I_p->apply(*mglevel_m[level-1]->error_p, *tmp);
        IpplTimings::stopTimer(interpTimer_m);
        
        // e^(l) += tmp
        mglevel_m[level]->error_p->update(1.0, *tmp, 1.0);
        
        // residual update
        tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        this->residual_no_fine_m(tmp,
                                 mglevel_m[level]->error_p,
                                 mglevel_m[level-1]->error_p,
                                 mglevel_m[level]->residual_p,
                                 level);
        
        *mglevel_m[level]->residual_p = *tmp;
        
        // delta error
        Teuchos::RCP<vector_t> derror = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        
        // smoothing
        IpplTimings::startTimer(smoothTimer_m);
        for (std::size_t iii = 0; iii < nSweeps_m[level]; ++iii)
            this->gsrb_level_m(derror, mglevel_m[level]->residual_p, level);
        IpplTimings::stopTimer(smoothTimer_m);
        
        // e^(l) += de^(l)
        mglevel_m[level]->error_p->update(1.0, *derror, 1.0);
        
        // phi^(l) = phi^(l, save) + e^(l)
        mglevel_m[level]->phi_p->update(1.0, phi_save, 1.0, *mglevel_m[level]->error_p, 0.0);
        
    } else {
        // e = A^(-1)r
        
        mglevel_m[level]->Anf_p->resumeFill();
        mglevel_m[level]->Anf_p->scale(-1.0);
        mglevel_m[level]->Anf_p->fillComplete();
        
        mglevel_m[level]->residual_p->scale(-1.0);
        
        
        double tol = mglevel_m[lfine_m]->residual_p->normInf();
        
        tol = std::abs(tol) * 1.0e-1;
        
        tol = std::min(tol, 1.0e-1);
        
        std::cout << tol << std::endl; //std::cin.get();
        
        solver_mp->solve(mglevel_m[level]->Anf_p,
                         mglevel_m[level]->error_p,
                         mglevel_m[level]->residual_p, tol);
        
        mglevel_m[level]->Anf_p->resumeFill();
        mglevel_m[level]->Anf_p->scale(-1.0);
        mglevel_m[level]->Anf_p->fillComplete();
        
        mglevel_m[level]->residual_p->scale(-1.0);
        
        // phi = phi + e
        
        mglevel_m[level]->phi_p->update(1.0, *mglevel_m[level]->error_p, 1.0);
    }
}


void AmrMultiGrid::residual_no_fine_m(Teuchos::RCP<vector_t>& result,
                                      const Teuchos::RCP<vector_t>& rhs,
                                      const Teuchos::RCP<vector_t>& crs_rhs,
                                      const Teuchos::RCP<vector_t>& b,
                                      int level)
{
    IpplTimings::startTimer(residnofineTimer_m);
    vector_t crse2fine(mglevel_m[level]->Anf_p->getDomainMap());
    
    // get boundary for 
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->apply(*crs_rhs, crse2fine);
    }
    
    // tmp = A * x
    vector_t tmp(mglevel_m[level]->Anf_p->getDomainMap());
    mglevel_m[level]->Anf_p->apply(*rhs, tmp);
    
    
    if ( level > lbase_m ) {
        vector_t t2mp(mglevel_m[level]->Anf_p->getDomainMap());
        mglevel_m[level]->B_p->apply(*rhs, t2mp);
    
        tmp.update(1.0, t2mp, 1.0);
    }
    
    tmp.update(1.0, crse2fine, 1.0);
    
    result->update(1.0, *b, -1.0, tmp, 0.0);
    IpplTimings::stopTimer(residnofineTimer_m);
}


double AmrMultiGrid::l2error_m() {
    int nLevel = lfine_m - lbase_m + 1;
    
    double err = 0.0;
    
    for (int lev = 0; lev < nLevel; ++lev) {
        
        double tmp = mglevel_m[lev]->residual_p->norm2();
        err += tmp;
    }
    
    return err;
}


void AmrMultiGrid::buildNoFinePoissonMatrix_m(int level) {
    /*
     * Laplacian of "no fine"
     */
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nPhysBoundary = 2 * BL_SPACEDIM * mglevel_m[level]->getBCStencilNum();
    
    int nEntries = (BL_SPACEDIM << 1) + 1 /* plus boundaries */ + nPhysBoundary;
    
    // number of internal stencil points
    int nIntBoundary = BL_SPACEDIM * interface_mp->getNumberOfPoints();
    
    amrex::BoxArray empty;
    
    mglevel_m[level]->Anf_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p, nEntries, Tpetra::StaticProfile) );
    
    if ( level > lbase_m ) {
        mglevel_m[level]->B_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p, nIntBoundary, Tpetra::StaticProfile) );
    }
    
    indices_t indices, bindices;
    coefficients_t values, bvalues;
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    indices.clear();
                    values.clear();
                    
                    bindices.clear();
                    bvalues.clear();
                    
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    /*
                     * check neighbours in all directions (Laplacian stencil --> cross)
                     */
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case AmrMultiGridLevel_t::Mask::COVERED:
                                    // covered --> interior cell
                                case AmrMultiGridLevel_t::Mask::INTERIOR:
                                {
                                    indices.push_back( mglevel_m[level]->serialize(biv) );
                                    values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                    break;
                                }
                                case AmrMultiGridLevel_t::Mask::BNDRY:
                                {
                                    // boundary cell
                                    // only level > 0 have this kind of boundary
                                    
                                    /* Dirichlet boundary conditions from coarser level.
                                     */
                                    std::size_t nn = bindices.size();
                                    
                                    interface_mp->fine(biv, bindices, bvalues, d, -shift, empty,
                                                       mglevel_m[level].get());
                                    
                                    double value = 1.0 / ( dx[d] * dx[d] );
                                    for (std::size_t iter = nn; iter < bindices.size(); ++iter)
                                        bvalues[iter] *= value;
                                    
                                    break;
                                }
                                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                                {
                                    // physical boundary cell
                                    double value = 1.0 / ( dx[d] * dx[d] );
                                    mglevel_m[level]->applyBoundary(biv,
                                                                    indices,
                                                                    values,
                                                                    value /*matrix coefficient*/);
                                    break;
                                }
                                default:
                                    throw std::runtime_error("Error in mask for level "
                                                             + std::to_string(level) + "!");
                                    break;
                            }
                        }
                    }
                    
                    // check center
                    indices.push_back( globalidx );
                    values.push_back( -2.0 / ( dx[0] * dx[0] ) +
                                      -2.0 / ( dx[1] * dx[1] )
#if BL_SPACEDIM == 3
                                      -2.0 / ( dx[2] * dx[2] )
#endif
                        );
                    
                    this->unique_m(indices, values);
                    
                    mglevel_m[level]->Anf_p->insertGlobalValues(globalidx,
                                                                indices.size(),
                                                                &values[0],
                                                                &indices[0]);
                    
                    if ( level > lbase_m && bindices.size() > 0 )
                        mglevel_m[level]->B_p->insertGlobalValues(globalidx,
                                                                  bindices.size(),
                                                                  &bvalues[0],
                                                                  &bindices[0]);
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    mglevel_m[level]->Anf_p->fillComplete();
    
    if ( level > lbase_m )
        mglevel_m[level]->B_p->fillComplete();
}


void AmrMultiGrid::buildCompositePoissonMatrix_m(int level) {
    /*
     * Laplacian of "with fine"
     * 
     * For the finest level: Awf == Anf
     */
    
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level]->refinement();
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    const double* fdx = 0;
    
    amrex::BoxArray crse_fine_ba; // empty box array
    
    if ( level < lfine_m ) {
        mglevel_m[level+1]->geom.CellSize();
        crse_fine_ba = mglevel_m[level]->grids;
        crse_fine_ba.refine(rr);
        crse_fine_ba = amrex::intersect(mglevel_m[level+1]->grids, crse_fine_ba);
        crse_fine_ba.coarsen(rr);
    }
    
    /*
     * 1D not supported by AmrMultiGrid
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    int nPhysBoundary = 2 * BL_SPACEDIM * mglevel_m[level]->getBCStencilNum();
    
    // number of internal stencil points
    int nIntBoundary = BL_SPACEDIM * interface_mp->getNumberOfPoints();
    
    int nEntries = (BL_SPACEDIM << 1) + 5 /* plus boundaries */ +
                   nPhysBoundary + nIntBoundary;
    
    mglevel_m[level]->Awf_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                         nEntries,
                                                         Tpetra::StaticProfile)
                                          );
    
    mglevel_m[level]->UnCovered_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                               1,
                                                               Tpetra::StaticProfile)
                                                );
    
    indices_t indices;
    coefficients_t values;
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) ) {
                        /*
                         * Only cells that are not refined
                         */
                        indices.clear();
                        values.clear();
                        int globalidx = mglevel_m[level]->serialize(iv);
                        
                        /*
                         * check neighbours in all directions (Laplacian stencil --> cross)
                         */
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                AmrIntVect_t biv = iv;                        
                                biv[d] += shift;
                                
                                if ( !crse_fine_ba.contains(biv) )
                                {
                                    /*
                                     * It can't be a refined cell!
                                     */
                                    switch ( mfab(biv) )
                                    {
                                        case AmrMultiGridLevel_t::Mask::COVERED:
                                            // covered --> interior cell
                                        case AmrMultiGridLevel_t::Mask::INTERIOR:
                                        {
                                            indices.push_back( mglevel_m[level]->serialize(biv) );
                                            values.push_back( 1.0 / ( dx[d] * dx[d] ) );
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        case AmrMultiGridLevel_t::Mask::BNDRY:
                                        {
                                            // boundary cell
                                            // only level > 0 have this kind of boundary
                                            
                                            /* We are on the fine side of the crse-fine interface
                                             * --> normal stencil --> no flux matching required
                                             * --> interpolation of fine ghost cell required
                                             * (used together with Bcrse)
                                             */
                                            
                                            
                                            /* Dirichlet boundary conditions from coarser level.
                                             */
                                            std::size_t nn = indices.size();
                                            
                                            interface_mp->fine(biv, indices, values, d, -shift, crse_fine_ba,
                                                               mglevel_m[level].get());
                                            
                                            double value = 1.0 / ( dx[d] * dx[d] );
                                            for (std::size_t iter = nn; iter < indices.size(); ++iter)
                                                values[iter] *= value;
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                                        {
                                            // physical boundary cell
                                            double value = 1.0 / ( dx[d] * dx[d] );
                                            
                                            mglevel_m[level]->applyBoundary(biv,
                                                                    indices,
                                                                    values,
                                                                    value /*matrix coefficient*/);
                                            
                                            // add center once
                                            indices.push_back( globalidx );
                                            values.push_back( -1.0 / ( dx[d] * dx[d] ) );
                                            
                                            break;
                                        }
                                        default:
                                            break;
                                    }
                                } else {
                                    /*
                                     * If neighbour cell is refined, we are on the coarse
                                     * side of the crse-fine interface --> flux matching
                                     * required --> interpolation of fine ghost cell
                                     * (used together with Bfine)
                                     */
                                    
                                    
                                    // flux matching, coarse part
                                    
                                    /* 2D --> 2 fine cells to compute flux per coarse-fine-interace --> avg = 2
                                     * 3D --> 4 fin cells to compute flux per coarse-fine-interace --> avg = 4
                                     * 
                                     * @precondition: refinement of 2
                                     */
                                    std::size_t nn = indices.size();
                                    
                                    // top and bottom for all directions
                                    for (int d1 = 0; d1 < 2; ++d1) {
#if BL_SPACEDIM == 3
                                        for (int d2 = 0; d2 < 2; ++d2) {
#endif
                                            
                                            /* in order to get a top iv --> needs to be odd value in "d"
                                             * in order to get a bottom iv --> needs to be even value in "d"
                                             */
                                            AmrIntVect_t fake(D_DECL(0, 0, 0));
                                    
                                            fake[(d+1)%BL_SPACEDIM] = d1;
#if BL_SPACEDIM == 3
                                            fake[(d+2)%BL_SPACEDIM] = d2;
#endif
                                            interface_mp->coarse(iv, indices, values, d, shift, crse_fine_ba,
                                                                 fake, mglevel_m[level].get());
                                    
#if BL_SPACEDIM == 3
                                        }
#endif
                                    }                                    
#if BL_SPACEDIM == 2
                                    double avg = 2;
#elif BL_SPACEDIM == 3
                                    double avg = 4;
#endif
                                    double value = -1.0 / ( avg * cdx[d] * fdx[d] );
                                    
                                    for (std::size_t iter = nn; iter < indices.size(); ++iter) {
                                        values[iter] *= value;
                                    }
                                }
                            }
                        }
                        
                        this->unique_m(indices, values);
                        
                        mglevel_m[level]->Awf_p->insertGlobalValues(globalidx,
                                                                    indices.size(),
                                                                    &values[0],
                                                                    &indices[0]);
                        
                        double vv = 1.0;
                        mglevel_m[level]->UnCovered_p->insertGlobalValues(globalidx,
                                                                          1,
                                                                          &vv,
                                                                          &globalidx);
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    mglevel_m[level]->Awf_p->fillComplete();
    
    mglevel_m[level]->UnCovered_p->fillComplete();
}


void AmrMultiGrid::buildRestrictionMatrix_m(int level) {
    /*
     * x^(l-1) = R * x^(l)
     */
    
    // base level does not need to have a restriction matrix
    if ( level == lbase_m )
        return;
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    
    
    // find all coarse cells that are covered by fine cells
    AmrIntVect_t rr = mglevel_m[level-1]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level-1]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
        
#if BL_SPACEDIM == 3
    int nNeighbours = rr[0] * rr[1] * rr[2];
#else
    int nNeighbours = rr[0] * rr[1];
#endif
    
    indices_t indices(nNeighbours);
    coefficients_t values(nNeighbours);
    
    double val = 1.0 / double(nNeighbours);
    
    mglevel_m[level]->R_p = Teuchos::rcp( new matrix_t(mglevel_m[level-1]->map_p, nNeighbours, Tpetra::StaticProfile) );
    
    for (amrex::MFIter mfi(mglevel_m[level-1]->grids,
                           mglevel_m[level-1]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            int ii = i * rr[0];
            for (int j = lo[1]; j <= hi[1]; ++j) {
                int jj = j * rr[1];
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    int kk = k * rr[2];
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( crse_fine_ba.contains(iv) ) {
                        int crse_globalidx = mglevel_m[level-1]->serialize(iv);
                    
                        int numEntries = 0;
                    
                        // neighbours
                        for (int iref = 0; iref < rr[0]; ++iref) {
                            for (int jref = 0; jref < rr[1]; ++jref) {
#if BL_SPACEDIM == 3
                                for (int kref = 0; kref < rr[2]; ++kref) {
#endif
                                    AmrIntVect_t riv(D_DECL(ii + iref, jj + jref, kk + kref));
                                    
                                    int fine_globalidx = mglevel_m[level]->serialize(riv);
                                    
                                    indices[numEntries] = fine_globalidx;
                                    values[numEntries]  = val;
                                    ++numEntries;
#if BL_SPACEDIM == 3
                                }
#endif
                            }
                        }
                    
                        mglevel_m[level]->R_p->insertGlobalValues(crse_globalidx,
                                                                  numEntries,
                                                                  &values[0],
                                                                  &indices[0]);
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
     mglevel_m[level]->R_p->fillComplete(mglevel_m[level]->map_p,
                                         mglevel_m[level-1]->map_p);
}


void AmrMultiGrid::buildInterpolationMatrix_m(int level) {
    /* crse (this): level
     * fine: level + 1
     */
    
    /*
     * This does not include ghost cells of *this* level
     * --> no boundaries
     * 
     * x^(l+1) = I * x^(l)
     */
    
    if ( level == lfine_m )
        return;
    
    int nNeighbours = (mglevel_m[level]->getBCStencilNum() + 1) * interp_mp->getNumberOfPoints();
    
    indices_t indices;
    coefficients_t values;
    
    mglevel_m[level]->I_p = Teuchos::rcp( new matrix_t(mglevel_m[level+1]->map_p, nNeighbours, Tpetra::StaticProfile) );
    
    for (amrex::MFIter mfi(mglevel_m[level+1]->grids,
                           mglevel_m[level+1]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx  = mfi.validbox();
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level+1]->serialize(iv);
                    
                    indices.clear();
                    values.clear();
                    
                    /*
                     * we need boundary + indices from coarser level
                     */
                    interp_mp->stencil(iv, indices, values, mglevel_m[level].get());
                    
                    this->unique_m(indices, values);
                    
                    mglevel_m[level]->I_p->insertGlobalValues(globalidx,
                                                              indices.size(),
                                                              &values[0],
                                                              &indices[0]);
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    mglevel_m[level]->I_p->fillComplete(mglevel_m[level]->map_p,
                                        mglevel_m[level+1]->map_p);
}


void AmrMultiGrid::buildCrseBoundaryMatrix_m(int level) {
    
    /*
     * fine (this): level
     * coarse:      level - 1
     */
    
    // the base level has only physical boundaries
    if ( level == lbase_m )
        return;
    
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level-1]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level-1]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    
    /* Find all fine cells that are at the crse-fine interace
     * 
     * Put them into a map (to avoid duplicates, e.g. due to corners).
     * Finally, iterate through the list of cells
     */
    // std::list<std::pair<int, int> > contains the shift and direction to come to the fine cell
    // not part of the box array
    map_t cells;
    
    for (amrex::MFIter mfi(mglevel_m[level]->grids, mglevel_m[level]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box&          bx  = mfi.validbox();
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    // iv is a fine cell
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    // check its neighbours to see if at crse-fine interface
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            // neighbour
                            AmrIntVect_t niv = iv;
                            niv[d] += shift;
                                
                            if ( !mglevel_m[level]->grids.contains(niv) ) {
                                // neighbour is does not belong to fine grids
                                
                                // insert iv if not yet exists
                                cells[iv].push_back(std::make_pair(shift, d));
                            }
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    int nNeighbours = 2 * BL_SPACEDIM * mglevel_m[level]->getBCStencilNum() *
                      2 * BL_SPACEDIM * interface_mp->getNumberOfPoints();
    
    mglevel_m[level]->Bcrse_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                           nNeighbours,
                                                           Tpetra::StaticProfile)
                                            );
    
    indices_t indices;
    coefficients_t values;
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    for (map_t::iterator it = cells.begin(); it != cells.end(); ++it) {
        // fine cell at crse-fine interface
        AmrIntVect_t iv = it->first;
        int globalidx = mglevel_m[level]->serialize(iv);
        
        indices.clear();
        values.clear();
        
        while ( !it->second.empty() ) {
            /*
             * "shift" is the amount of cells to fine cell that does not belong to fine grids
             * "d" is the direction to shift
             */
            int shift = it->second.front().first;
            int d     = it->second.front().second;
            it->second.pop_front();
            
            // coarse cell that is not refined
            AmrIntVect_t civ = iv;
            civ[d] += shift;
            civ.coarsen(mglevel_m[level]->refinement());
            
            std::size_t numEntries = indices.size();
            // we need boundary + indices from coarser level
            interface_mp->coarse(civ, indices, values, d, shift, crse_fine_ba,
                                 iv, mglevel_m[level-1].get());
            
            // we need normalization by mesh size squared
            for (std::size_t n = numEntries; n < indices.size(); ++n)
                values[n] *= 1.0 / ( dx[d] * dx[d] );
        }
        
        /* In some cases indices occur several times, e.g. for corner points
         * at the physical boundary --> sum them up
         */
        this->unique_m(indices, values);
        
        mglevel_m[level]->Bcrse_p->insertGlobalValues(globalidx,
                                                      indices.size(),
                                                      &values[0],
                                                      &indices[0]);
    }
    
    mglevel_m[level]->Bcrse_p->fillComplete(mglevel_m[level-1]->map_p,  // col map
                                            mglevel_m[level]->map_p);   // row map
}


void AmrMultiGrid::buildFineBoundaryMatrix_m(int level)
{
    /* fine: level + 1
     * coarse (this): level
     */
    
    // the finest level does not need data from a finer level
    if ( level == lfine_m )
        return;
    
    // find all cells with refinement
    AmrIntVect_t rr = mglevel_m[level]->refinement();
    
    amrex::BoxArray crse_fine_ba = mglevel_m[level]->grids;
    crse_fine_ba.refine(rr);
    crse_fine_ba = amrex::intersect(mglevel_m[level+1]->grids, crse_fine_ba);
    crse_fine_ba.coarsen(rr);
    
    // 2 interfaces per direction --> 2 * rr
    int nNeighbours = 4 * rr[0] * rr[1];
#if BL_SPACEDIM == 3
    nNeighbours = 6 * rr[0] * rr[1] * rr[2];
#endif
    
    mglevel_m[level]->Bfine_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                           nNeighbours, Tpetra::StaticProfile) );
    
    indices_t indices;
    coefficients_t values;
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    const double* fdx = mglevel_m[level+1]->geom.CellSize();
    
    /* Find all coarse cells that are at the crse-fine interace but are
     * not refined.
     * Put them into a map (to avoid duplicates, e.g. due to corners).
     * Finally, iterate through the list of cells
     */
    
    // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
    map_t cells;
    
    for (amrex::MFIter mfi(mglevel_m[level]->grids, mglevel_m[level]->dmap, false);
         mfi.isValid(); ++mfi)
    {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
            
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    if ( !crse_fine_ba.contains(iv) && mfab(iv) != 2 /*not physical bc*/ ) {
                        /*
                         * iv is a coarse cell that got not refined
                         * 
                         * --> check all neighbours to see if at crse-fine
                         * interface
                         */
                        for (int d = 0; d < BL_SPACEDIM; ++d) {
                            for (int shift = -1; shift <= 1; shift += 2) {
                                // neighbour
                                AmrIntVect_t covered = iv;
                                covered[d] += shift;
                                
                                if ( crse_fine_ba.contains(covered) &&
                                     !mglevel_m[level]->isBoundary(covered) /*not physical bc*/ )
                                {
                                    // neighbour is covered by fine cells
                                    
                                    // insert if not yet exists
                                    cells[iv].push_back(std::make_pair(shift, d));
                                }
                            }
                        }
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    auto fill = [&](indices_t& indices,
                    coefficients_t& values,
                    D_DECL(int ii, int jj, int kk),
                    int* begin, int* end, int d,
                    const AmrIntVect_t& iv, int shift,
                    int sign)
    {
        // number of fine cell gradients
        double avg = 1.0;
        switch ( d ) {
            case 0:
                // horizontal
                avg = rr[1];
#if BL_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
            case 1:
                // vertical
                avg = rr[0];
#if BL_SPACEDIM == 3
                avg *= rr[2];
#endif
                break;
#if BL_SPACEDIM == 3
            case 2:
                avg = rr[0] * rr[1];
                break;
#endif
            default:
                throw std::runtime_error("This dimension does not exist!");
                break;
        }
        
        for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
            for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
#if BL_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     */
                    AmrIntVect_t riv(D_DECL(iref, jref, kref));
                    
                    if ( riv[d] / rr[d] == iv[d] ) {
                        /* the fine cell is on the coarse side --> fine
                         * ghost cell --> we need to interpolate
                         */
                        double value = -1.0 / ( avg * cdx[d] * fdx[d] );
                        
                        std::size_t nn = indices.size();
                        
                        interface_mp->fine(riv, indices, values, d, shift, crse_fine_ba,
                                           mglevel_m[level+1].get());
                        
                        for (std::size_t iter = nn; iter < indices.size(); ++iter)
                            values[iter] *= value;
                        
                    } else {
                        double value = 1.0 / ( avg * cdx[d] * fdx[d] );
                        
                        indices.push_back( mglevel_m[level+1]->serialize(riv) );
                        values.push_back( value );
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    };
    
    for (map_t::iterator it = cells.begin(); it != cells.end(); ++it) {
        // not covered coarse cell at crse-fine interface
        AmrIntVect_t iv = it->first;
        int globalidx = mglevel_m[level]->serialize(iv);
        
        indices.clear();
        values.clear();
        
        while ( !it->second.empty() ) {
            /*
             * "shift" is the amount of is a coarse cell that got refined
             * "d" is the direction to shift
             * 
             * --> check all covered neighbour cells
             */
            int shift = it->second.front().first;
            int d     = it->second.front().second;
            it->second.pop_front();
            
            
            /* we need to iterate over correct fine cells. It depends
             * on the orientation of the interface
             */
            int begin[BL_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
            int end[BL_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
            
            // neighbour
            AmrIntVect_t covered = iv;
            covered[d] += shift;
                    
            /*
             * neighbour cell got refined but is not on physical boundary
             * --> we are at a crse-fine interface
             * 
             * we need now to find out which fine cells
             * are required to satisfy the flux matching
             * condition
             */
            
            switch ( shift ) {
                case -1:
                {
                    // --> interface is on the lower face
                    int ii = iv[0] * rr[0];
                    int jj = iv[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = iv[2] * rr[2];
#endif
                    // iterate over all fine cells at the interface
                    // start with lower cells --> cover coarse neighbour
                    // cell
                    fill(indices, values, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
                case 1:
                default:
                {
                    // --> interface is on the upper face
                    int ii = covered[0] * rr[0];
                    int jj = covered[1] * rr[1];
#if BL_SPACEDIM == 3
                    int kk = covered[2] * rr[2];
#endif
                    fill(indices, values, D_DECL(ii, jj, kk), &begin[0], &end[0], d, iv, shift, 1.0);
                    break;
                }
            }
        }
        
        this->unique_m(indices, values);
        
        mglevel_m[level]->Bfine_p->insertGlobalValues(globalidx,
                                                      indices.size(),
                                                      &values[0],
                                                      &indices[0]);
    }
    
    mglevel_m[level]->Bfine_p->fillComplete(mglevel_m[level+1]->map_p,
                                            mglevel_m[level]->map_p);
}


inline void AmrMultiGrid::buildDensityVector_m(const AmrField_t& rho, int level) {
    this->amrex2trilinos_m(rho, 0, mglevel_m[level]->rho_p, level);
}


inline void AmrMultiGrid::buildPotentialVector_m(const AmrField_t& phi, int level) {
    this->amrex2trilinos_m(phi, 0, mglevel_m[level]->phi_p, level);
}


void AmrMultiGrid::buildSmootherMatrix_m(int level) {
    // the base level does not require a smoother matrix
    if ( level == lbase_m )
        return;
    
    mglevel_m[level]->S_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p, 1, Tpetra::StaticProfile) );
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    double h = std::max(dx[0], dx[1]);
#if BL_SPACEDIM == 3
    h = std::max(h, dx[2]);
#endif
    h = h * h;
    
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    /*
                     * check all directions (Laplacian stencil --> cross)
                     */
                    bool interior = true;
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t biv = iv;                        
                            biv[d] += shift;
                            
                            switch ( mfab(biv) )
                            {
                                case AmrMultiGridLevel_t::Mask::COVERED:
                                    // covered --> interior cell
                                case AmrMultiGridLevel_t::Mask::INTERIOR:
                                    // interior cells
                                    interior *= true;
                                    break;
                                case AmrMultiGridLevel_t::Mask::BNDRY:
                                    // boundary cells
                                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                                    // boundary cells
                                    interior *= false;
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                    
#if BL_SPACEDIM == 2
                    double value = h * ( (interior) ? 0.25 : 0.1875 );
#elif BL_SPACEDIM == 3
                    double value = h * ( (interior) ? 1.0 / 6.0 : 1.0 / 24.0 );
#endif
                    mglevel_m[level]->S_p->insertGlobalValues(globalidx,
                                                              1,
                                                              &value,
                                                              &globalidx);
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    mglevel_m[level]->S_p->fillComplete();
}


void AmrMultiGrid::buildGradientMatrix_m(int level) {
    
    int nEntries = 5;
    
    for (int d = 0; d < BL_SPACEDIM; ++d) {
        mglevel_m[level]->G_p[d] = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p, nEntries, Tpetra::StaticProfile) );
    }
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    indices_t indices;
    coefficients_t values;
    
    auto check = [&](const AmrIntVect_t& iv,
                     const amrex::BaseFab<int>& mfab,
                     int dir,
                     double shift)
    {
        switch ( mfab(iv) )
        {
            case AmrMultiGridLevel_t::Mask::COVERED:
                // covered --> interior cell
            case AmrMultiGridLevel_t::Mask::INTERIOR:
                // interior cells
                
                indices.push_back( mglevel_m[level]->serialize(iv) );
                values.push_back(  - shift * 0.5 / dx[dir] );
                
                break;
            case AmrMultiGridLevel_t::Mask::BNDRY:
            {
                // interior boundary cells --> only level > 0
                
                double value = - shift * 0.5 / dx[dir];
                
                // use 1st order Lagrange --> only cells of this level required
                AmrIntVect_t tmp = iv;
                // first fine cell on refined coarse cell (closer to interface)
                tmp[dir] -= shift;
                indices.push_back( mglevel_m[level]->serialize(tmp) );
                values.push_back( 2.0 * value );
                
                // second fine cell on refined coarse cell (further away from interface)
                tmp[dir] -= shift;
                indices.push_back( mglevel_m[level]->serialize(tmp) );
                values.push_back( -1.0 * value );
                
                break;
            }
            case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
            {
                // physical boundary cells
                
                double value = - shift * 0.5 / dx[dir];
                
                mglevel_m[level]->applyBoundary(iv,
                                                indices,
                                                values,
                                                value);
                break;
            }
            default:
                break;
        }
    };
    
    for (amrex::MFIter mfi(*mglevel_m[level]->mask, false); mfi.isValid(); ++mfi) {
        const amrex::Box& bx  = mfi.validbox();
        const amrex::BaseFab<int>& mfab = (*mglevel_m[level]->mask)[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    for (int d = 0; d < BL_SPACEDIM; ++d) {
                        
                        indices.clear();
                        values.clear();
                        
                        for (int shift = -1; shift <= 1; shift += 2) {
                            AmrIntVect_t niv = iv;                        
                            niv[d] += shift;
                            check(niv, mfab, d, shift);
                        }
                        
                        mglevel_m[level]->G_p[d]->insertGlobalValues(globalidx,
                                                                     indices.size(),
                                                                     &values[0],
                                                                     &indices[0]);
                    }
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
    
    
    for (int d = 0; d < BL_SPACEDIM; ++d)
        mglevel_m[level]->G_p[d]->fillComplete();
}


void AmrMultiGrid::amrex2trilinos_m(const AmrField_t& mf, int comp,
                                    Teuchos::RCP<vector_t>& mv, int level)
{
    if ( mv.is_null() )
        mv = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        const amrex::FArrayBox&    fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    mv->replaceGlobalValue(globalidx, fab(iv, comp));
#if BL_SPACEDIM == 3
                }
#endif
            }
        }
    }
}


void AmrMultiGrid::trilinos2amrex_m(AmrField_t& mf, int comp,
                                    const Teuchos::RCP<vector_t>& mv)
{
    Teuchos::ArrayRCP<const amr::scalar_t> data =  mv->get1dView();
    
    int localidx = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box&          bx  = mfi.validbox();
        amrex::FArrayBox&          fab = mf[mfi];
        
        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();
        
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
#if BL_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    fab(iv, comp) = data[localidx++];
                }
#if BL_SPACEDIM == 3
            }
#endif
        }
    }
}


void AmrMultiGrid::unique_m(indices_t& indices,
                            coefficients_t& values)
{
    indices_t uindices;
    coefficients_t uvalues;
    
    // we need to sort for std::unique_copy
    // 20. Sept. 2017,
    // https://stackoverflow.com/questions/34878329/how-to-sort-two-vectors-simultaneously-in-c-without-using-boost-or-creating-te
    std::vector< std::pair<int, double> > pair;
    for (uint i = 0; i < indices.size(); ++i)
        pair.push_back( { indices[i], values[i] } );
    
    std::sort(pair.begin(), pair.end(),
              [](const std::pair<int, double>& a,
                  const std::pair<int, double>& b) {
                  return a.first < b.first;
              });
    
    for (uint i = 0; i < indices.size(); ++i) {
        indices[i] = pair[i].first;
        values[i]  = pair[i].second;
    }
    
    // requirement: duplicates are consecutive
    std::unique_copy(indices.begin(), indices.end(), std::back_inserter(uindices));
    
    uvalues.resize(uindices.size());
    
    for (std::size_t i = 0; i < uvalues.size(); ++i) {
        for (std::size_t j = 0; j < values.size(); ++j) {
            if ( uindices[i] == indices[j] )
                uvalues[i] += values[j];
        }
    }
    
    std::swap(values, uvalues);
    std::swap(indices, uindices);
}


void AmrMultiGrid::gsrb_level_m(Teuchos::RCP<vector_t>& e,
                                Teuchos::RCP<vector_t>& r,
                                int level)
{
    // apply "no fine" Laplacian
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
    
    // tmp = A * x
    mglevel_m[level]->Anf_p->apply(/**mglevel_m[level]->error_p*/*e, *tmp);
    
    // tmp = tmp - r
    tmp->update(-1.0, *r, 1.0);
    
    // tmp2 = S * tmp
    Teuchos::RCP<vector_t> tmp2 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    mglevel_m[level]->S_p->apply(*tmp, *tmp2);
    
    e->update(1.0, *tmp2, 1.0);
}


void AmrMultiGrid::restrict_m(int level) {
    
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level+1]->map_p) );
    
    this->residual_no_fine_m(tmp,
                             mglevel_m[level+1]->error_p,
                             mglevel_m[level]->error_p,
                             mglevel_m[level+1]->residual_p, level+1);
    
//     std::cout << *tmp << std::endl;
    
    // average down: residual^(l-1) = R^(l) * tmp
    mglevel_m[level+1]->R_p->apply(*tmp, *mglevel_m[level]->residual_p);
    
    
//     std::cout << *mglevel_m[level]->residual_p << std::endl;
    
    // special matrix, i.e. matrix without covered cells
    // r^(l-1) = rho^(l-1) - A * phi^(l-1)
    
    vector_t fine2crse(mglevel_m[level]->Anf_p->getDomainMap());
    
    // get boundary for 
    mglevel_m[level]->Bfine_p->apply(*mglevel_m[level+1]->phi_p, fine2crse);
    
    vector_t tmp2(mglevel_m[level]->Awf_p->getDomainMap());
    mglevel_m[level]->Awf_p->apply(*mglevel_m[level]->phi_p, tmp2);
    
    vector_t crse2fine(mglevel_m[level]->Anf_p->getDomainMap());
    
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p, crse2fine);
    }
    
    tmp2.update(1.0, fine2crse, 1.0, crse2fine, 1.0);
    
    Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
    
    mglevel_m[level]->UnCovered_p->apply(*mglevel_m[level]->rho_p, *tmp3);
    
    // ONLY subtract coarse rho
    mglevel_m[level]->residual_p->update(1.0, *tmp3, -1.0, tmp2, 1.0);
}


void AmrMultiGrid::averageDown_m(int level) {
    
    if (level == lfine_m )
        return;
    
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    
    mglevel_m[level+1]->R_p->apply(*mglevel_m[level+1]->phi_p, *tmp);
    
    Teuchos::RCP<vector_t> tmp2 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    
    mglevel_m[level]->UnCovered_p->apply(*mglevel_m[level]->phi_p, *tmp2);
    
    mglevel_m[level]->phi_p->update(1.0, *tmp, 1.0, *tmp2, 0.0);
}


void AmrMultiGrid::initInterpolater_m(const Interpolater& interp) {
    switch ( interp ) {
        case Interpolater::TRILINEAR:
            interp_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            std::runtime_error("Not yet implemented.");
        case Interpolater::PIECEWISE_CONST:
            interp_mp.reset( new AmrPCInterpolater<AmrMultiGridLevel_t>() );
            break;
        default:
            std::runtime_error("No such interpolater available.");
    }
}


void AmrMultiGrid::initCrseFineInterp_m(const Interpolater& interface) {
    switch ( interface ) {
        case Interpolater::TRILINEAR:
            interface_mp.reset( new AmrTrilinearInterpolater<AmrMultiGridLevel_t>() );
            break;
        case Interpolater::LAGRANGE:
            interface_mp.reset( new AmrLagrangeInterpolater<AmrMultiGridLevel_t>(
                AmrLagrangeInterpolater<AmrMultiGridLevel_t>::Order::QUADRATIC) );
            break;
        case Interpolater::PIECEWISE_CONST:
            interface_mp.reset( new AmrPCInterpolater<AmrMultiGridLevel_t>() );
            break;
        default:
            std::runtime_error("No such interpolater for the coarse-fine interface available.");
    }
}


void AmrMultiGrid::initBaseSolver_m(const LinSolver& solver) {
    switch ( solver ) {
        case LinSolver::BLOCK_CG:
            solver_mp.reset( new BlockCGSolMgr() );
            break;
        default:
            std::runtime_error("No such solver available.");
    }
}
