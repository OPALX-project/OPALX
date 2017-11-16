#include "AmrMultiGrid.h"

#define AMR_MG_WRITE 1

#define DEBUG 1

#if AMR_MG_WRITE
    #include <iomanip>
    #include <fstream>
#endif

AmrMultiGrid::AmrMultiGrid(Boundary bc,
                           Interpolater interp,
                           Interpolater interface,
                           BaseSolver solver,
                           Preconditioner precond,
                           Smoother smoother,
                           Norm norm
                          )
    : nIter_m(0),
      nSweeps_m(12),
      smootherType_m(smoother),
      lbase_m(0),
      lfine_m(0),
      bc_m(bc),
      norm_m(norm)
{
    comm_mp = Teuchos::rcp( new comm_t( Teuchos::opaqueWrapper(Ippl::getComm()) ) );
    node_mp = KokkosClassic::DefaultNode::getDefaultNode();
    
#if AMR_MG_TIMER
    buildTimer_m        = IpplTimings::getTimer("build");
    restrictTimer_m     = IpplTimings::getTimer("restrict");
    smoothTimer_m       = IpplTimings::getTimer("smooth");
    interpTimer_m       = IpplTimings::getTimer("prolongate");
    residnofineTimer_m  = IpplTimings::getTimer("resid-no-fine");
    bottomTimer_m       = IpplTimings::getTimer("bottom-solver");
#endif
    
    this->initInterpolater_m(interp);
    
    // interpolater for crse-fine-interface
    this->initCrseFineInterp_m(interface);
    
    // base level solver
    this->initBaseSolver_m(solver, precond);
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
    nlevel_m = lfine - lbase + 1;
    
    this->initLevels_m(rho, geom);
    
    this->initGuess_m(phi, previous);
    
    // build all necessary matrices and vectors
    this->setup_m(rho, phi);
    
    // actual solve
    this->iterate_m();
    
    // write efield to AMReX
    this->computeEfield_m(efield);    
    
    // copy solution back
    for (int lev = 0; lev < nlevel_m; ++lev) {
        int ilev = lbase + lev;
        
        this->trilinos2amrex_m(*phi[ilev], 0, mglevel_m[lev]->phi_p);
    }
}


void AmrMultiGrid::initLevels_m(const amrex::Array<AmrField_u>& rho,
                                const amrex::Array<AmrGeometry_t>& geom)
{
    mglevel_m.resize(nlevel_m);
    
    AmrIntVect_t rr = AmrIntVect_t(D_DECL(2, 2, 2));
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        int ilev = lbase_m + lev;
        
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
}


void AmrMultiGrid::initGuess_m(amrex::Array<AmrField_u>& phi, bool previous) {
    if ( !previous ) {
        // reset
        for (int lev = 0; lev < nlevel_m; ++lev) {
            int ilev = lbase_m + lev;
            phi[ilev]->setVal(0.0, phi[ilev]->nGrow());
        }
    }
}


void AmrMultiGrid::iterate_m() {
    
    double eps = 1.0e-8;
    
    // initial error
    double max_residual = 0.0;
    double max_rho = 0.0;
    
    this->initResidual_m(max_residual, max_rho);
    
    while ( max_residual > eps * max_rho) {
        
        relax_m(lfine_m);
        
        // update residual
        for (int lev = 0; lev < nlevel_m; ++lev) {
            
            this->residual_m(mglevel_m[lev]->residual_p,
                             mglevel_m[lev]->rho_p,
                             mglevel_m[lev]->phi_p, lev);
        }
        
        max_residual = residualNorm_m();
        
        ++nIter_m;
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
        
        vector_t fine2crse(mglevel_m[level]->Awf_p->getDomainMap());
        
        // get boundary for 
        if ( mglevel_m[level]->Bfine_p != Teuchos::null ) {
            mglevel_m[level]->Bfine_p->apply(*mglevel_m[level+1]->phi_p, fine2crse);
        }
        
        vector_t tmp2(mglevel_m[level]->Awf_p->getDomainMap());
        mglevel_m[level]->Awf_p->apply(*x, tmp2);
        
        vector_t crse2fine(mglevel_m[level]->Awf_p->getDomainMap());
        
        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p, crse2fine);
        }
        
        tmp2.update(1.0, fine2crse, 1.0, crse2fine, 1.0);
        
        Teuchos::RCP<vector_t> tmp4 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
        mglevel_m[level]->UnCovered_p->apply(tmp2, *tmp4);
        
        Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
    
        mglevel_m[level]->UnCovered_p->apply(*b, *tmp3);
    
        // ONLY subtract coarse rho
//         mglevel_m[level]->residual_p->putScalar(0.0);
        
        r->update(1.0, *tmp3, -1.0, *tmp4, 0.0);
        
    } else {
        /* finest level: Awf_p == Anf_p
         * 
         * In this case we use Awf_p instead of Anf_p since Anf_p might be
         * made positive definite for the bottom solver.
         */
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
        mglevel_m[level]->Awf_p->apply(*x, *tmp);
        
        if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
            Teuchos::RCP<vector_t> crse2fine = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
            mglevel_m[level]->Bcrse_p->apply(*mglevel_m[level-1]->phi_p, *crse2fine);
            tmp->update(1.0, *crse2fine, 1.0);
        }
        r->update(1.0, *b, -1.0, *tmp, 0.0);        
    }
}


void AmrMultiGrid::relax_m(int level) {
    
    if ( level == lfine_m ) {
        
        if ( level == lbase_m ) {
            /* Anf_p == Awf_p
             * 
             * In this case we use Awf_p instead of Anf_p since Anf_p might be
             * made positive definite for the bottom solver.
             */
            Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, true) );
            mglevel_m[level]->Awf_p->apply(*mglevel_m[level]->phi_p, *tmp);
            mglevel_m[level]->residual_p->update(1.0, *mglevel_m[level]->rho_p, -1.0, *tmp, 0.0);
        } else {
            this->residual_no_fine_m(mglevel_m[level]->residual_p,
                                     mglevel_m[level]->phi_p,
                                     mglevel_m[level-1]->phi_p,
                                     mglevel_m[level]->rho_p, level);
        }
    }
    
    if ( level > 0 ) {
        // phi^(l, save) = phi^(l)        
        Teuchos::RCP<vector_t> phi_save = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        Tpetra::deep_copy(*phi_save, *mglevel_m[level]->phi_p);
        
        mglevel_m[level-1]->error_p->putScalar(0.0);
        
        // smoothing
        this->smooth_m(mglevel_m[level]->error_p,
                       mglevel_m[level]->residual_p, level);
        
        
        // phi = phi + e
        mglevel_m[level]->phi_p->update(1.0, *mglevel_m[level]->error_p, 1.0);
        
        /*
         * restrict
         */
        this->restrict_m(level);
        
        this->relax_m(level - 1);
        
        /*
         * prolongate / interpolate
         */
        this->prolongate_m(level);
        
        // residual update
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        this->residual_no_fine_m(tmp,
                                 mglevel_m[level]->error_p,
                                 mglevel_m[level-1]->error_p,
                                 mglevel_m[level]->residual_p,
                                 level);
        
        Tpetra::deep_copy(*mglevel_m[level]->residual_p, *tmp);
        
        // delta error
        Teuchos::RCP<vector_t> derror = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        
        // smoothing
        this->smooth_m(derror, mglevel_m[level]->residual_p, level);
        
        // e^(l) += de^(l)
        mglevel_m[level]->error_p->update(1.0, *derror, 1.0);
        
        // phi^(l) = phi^(l, save) + e^(l)
        mglevel_m[level]->phi_p->update(1.0, *phi_save, 1.0, *mglevel_m[level]->error_p, 0.0);
        
    } else {
        // e = A^(-1)r
#if AMR_MG_TIMER
    IpplTimings::startTimer(bottomTimer_m);
#endif
    
        solver_mp->solve(mglevel_m[level]->error_p,
                         mglevel_m[level]->residual_p);
        
#if AMR_MG_TIMER
    IpplTimings::stopTimer(bottomTimer_m);
#endif
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
#if AMR_MG_TIMER
    IpplTimings::startTimer(residnofineTimer_m);
#endif
    vector_t crse2fine(mglevel_m[level]->Anf_p->getDomainMap());
    
    // get boundary for 
    if ( mglevel_m[level]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level]->Bcrse_p->apply(*crs_rhs, crse2fine);
    }
    
    
    // tmp = A * x
    vector_t tmp(mglevel_m[level]->Anf_p->getDomainMap());
    mglevel_m[level]->Anf_p->apply(*rhs, tmp);
    
    tmp.update(1.0, crse2fine, 1.0);
    
    result->update(1.0, *b, -1.0, tmp, 0.0);
#if AMR_MG_TIMER
    IpplTimings::stopTimer(residnofineTimer_m);
#endif
}

double AmrMultiGrid::residualNorm_m() {
    double err = 0.0;
    
#if AMR_MG_WRITE
    std::ofstream out;
    if ( Ippl::myNode() == 0 )
        out.open("residual.dat", std::ios::app);
#endif
    
    double tmp = 0.0;
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        
        switch ( norm_m ) {
            case Norm::L1:
            {
                tmp = mglevel_m[lev]->residual_p->norm1();
                break;
            }
            case Norm::L2:
            {
                tmp = mglevel_m[lev]->residual_p->norm2();
                break;
            }
            case Norm::LINF:
            {
                tmp = mglevel_m[lev]->residual_p->normInf();
                break;
            }
        }
        
        err = std::max(err, tmp);
        
#if AMR_MG_WRITE
        if ( Ippl::myNode() == 0 )
            out << std::setw(15) << std::right << tmp;
#endif
    }
    
#if AMR_MG_WRITE
    if ( Ippl::myNode() == 0 ) {
        out << std::setw(15) << err << std::endl;
        out.close();
    }
#endif
    
    
    return err;
}


void AmrMultiGrid::initResidual_m(double& maxResidual, double& maxRho) {
#if AMR_MG_WRITE
    std::ofstream out;
    
    if ( Ippl::myNode() == 0) {
        out.open("residual.dat", std::ios::out);
    
        for (int lev = 0; lev < nlevel_m; ++lev)
            out << std::setw(14) << std::right << "level" << lev;
        out << std::setw(15) << std::right << "max";
        out << std::endl;
    }
#endif
    
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
        this->residual_m(mglevel_m[lev]->residual_p,
                         mglevel_m[lev]->rho_p,
                         mglevel_m[lev]->phi_p, lev);
        
        double tmp = mglevel_m[lev]->residual_p->normInf();
        maxResidual = std::max(maxResidual, tmp);
        
#if AMR_MG_WRITE
        if ( Ippl::myNode() == 0 )
            out << std::setw(15) << std::right << tmp;
#endif
        
        tmp = mglevel_m[lev]->rho_p->normInf();
        maxRho = std::max(maxRho, tmp);
    }
    
#if AMR_MG_WRITE
    if ( Ippl::myNode() == 0 ) {
        out << std::setw(15) << maxResidual << std::endl;
        out.close();
    }
#endif
}


void AmrMultiGrid::computeEfield_m(amrex::Array<AmrField_u>& efield) {
    Teuchos::RCP<vector_t> efield_p = Teuchos::null;
    for (int lev = nlevel_m - 1; lev > -1; --lev) {
        int ilev = lbase_m + lev;
        
        efield_p = Teuchos::rcp( new vector_t(mglevel_m[lev]->map_p, false) );
        
        averageDown_m(lev);
        
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            mglevel_m[lev]->G_p[d]->apply(*mglevel_m[lev]->phi_p, *efield_p);
            this->trilinos2amrex_m(*efield[ilev], d, efield_p);
        }
    }
}


void AmrMultiGrid::setup_m(const amrex::Array<AmrField_u>& rho,
                           const amrex::Array<AmrField_u>& phi,
                           bool matrices)
{
#if AMR_MG_TIMER
    IpplTimings::startTimer(buildTimer_m);
#endif
        
    // the base level has no smoother --> nlevel_m - 1
    smoother_m.resize(nlevel_m-1);
    
    for (int lev = 0; lev < nlevel_m; ++lev) {
    
        this->open_m(lev, matrices);
        
        int ilev = lbase_m + lev;
        
        // find all coarse cells that are covered by fine cells
        AmrIntVect_t rr = mglevel_m[lev]->refinement();
        
        amrex::BoxArray this_fine_ba;   // empty
        
        if ( lev < lfine_m ) {
            this_fine_ba = mglevel_m[lev]->grids;
            this_fine_ba.refine(rr);
            this_fine_ba = amrex::intersect(mglevel_m[lev+1]->grids, this_fine_ba);
            this_fine_ba.coarsen(rr);
        }
        
        
        /* Find all fine cells that are at the crse-fine interace
         * 
         * Put them into a map (to avoid duplicates, e.g. due to corners).
         * Finally, iterate through the list of cells
         */
        // std::list<std::pair<int, int> > contains the shift and direction to come to the fine cell
        // not part of the box array
        map_t cells_fine;
        
        /* Find all coarse cells that are at the crse-fine interace but are
         * not refined.
         * Put them into a map (to avoid duplicates, e.g. due to corners).
         * Finally, iterate through the list of cells
         */
        // std::list<std::pair<int, int> > contains the shift and direction to come to the covered cell
        map_t cells_crse;
        
        for (amrex::MFIter mfi(*mglevel_m[lev]->mask, false);
             mfi.isValid(); ++mfi)
        {
            const box_t&       bx  = mfi.validbox();
            const basefab_t&  mfab = (*mglevel_m[lev]->mask)[mfi];
            const farraybox_t& rfab = (*rho[ilev])[mfi];
            const farraybox_t& pfab = (*phi[ilev])[mfi];
            
            const int* lo = bx.loVect();
            const int* hi = bx.hiVect();
        
            for (int i = lo[0]; i <= hi[0]; ++i) {
                int ii = i * rr[0];
                for (int j = lo[1]; j <= hi[1]; ++j) {
                    int jj = j * rr[1];
#if AMREX_SPACEDIM == 3
                    for (int k = lo[2]; k <= hi[2]; ++k) {
                        int kk = k * rr[2];
#endif
                        AmrIntVect_t iv(D_DECL(i, j, k));
                        
                        this->buildRestrictionMatrix_m(iv, this_fine_ba,
                                                       D_DECL(ii, jj, kk),
                                                       lev);
                        
                        this->buildInterpolationMatrix_m(iv, lev);
                        
                        this->buildCrseBoundaryMatrix_m(iv, cells_fine, lev);
                        
                        this->buildFineBoundaryMatrix_m(iv, cells_crse, this_fine_ba, lev);
                        
                        this->buildNoFinePoissonMatrix_m(iv, mfab, lev);
                        
                        this->buildCompositePoissonMatrix_m(iv, mfab, this_fine_ba, lev);
                        
                        this->buildGradientMatrix_m(iv, mfab, lev);
                        
                        int globalidx = mglevel_m[lev]->serialize(iv);
                        mglevel_m[lev]->rho_p->replaceGlobalValue(globalidx, rfab(iv, 0));
                        mglevel_m[lev]->phi_p->replaceGlobalValue(globalidx, pfab(iv, 0));
                    }
                }
            }
        }
        
        
        if ( lev > lbase_m ) {
            amrex::BoxArray crse_this_ba = mglevel_m[lev-1]->grids;
            crse_this_ba.refine(rr);
            crse_this_ba = amrex::intersect(mglevel_m[lev]->grids, crse_this_ba);
            crse_this_ba.coarsen(rr);
            
            this->fillCrseBoundaryMatrix_m(cells_fine, crse_this_ba, lev);
        }
        
        this->fillFineBoundaryMatrix_m(cells_crse, this_fine_ba, lev);
        
        this->close_m(lev, matrices);
        
        if ( lev > lbase_m ) {
            smoother_m[lev-1].reset( new AmrSmoother(mglevel_m[lev]->Anf_p,
                                                     smootherType_m, nSweeps_m) );
        }
    }
    
    if ( matrices ) {
        // set the bottom solve operator
        solver_mp->setOperator(mglevel_m[lbase_m]->Anf_p);
    }
    
    mglevel_m[lfine_m]->error_p->putScalar(0.0);
    
#if AMR_MG_TIMER
    IpplTimings::stopTimer(buildTimer_m);
#endif
}


void AmrMultiGrid::open_m(int level, bool matrices) {
    
    if ( matrices ) {
    
        if ( level > lbase_m ) {
            
            /*
             * interpolation matrix
             */
            
            int nNeighbours = (mglevel_m[level]->getBCStencilNum() + 1) * interp_mp->getNumberOfPoints();
    
            mglevel_m[level]->I_p = Teuchos::rcp( new matrix_t(mglevel_m[level]->map_p,
                                                               nNeighbours,
                                                               Tpetra::StaticProfile) );
            
            /*
             * coarse boundary matrix
             */
            
            nNeighbours = 2 * AMREX_SPACEDIM * mglevel_m[level]->getBCStencilNum() *
                          2 * AMREX_SPACEDIM * interface_mp->getNumberOfPoints();
            
            mglevel_m[level]->Bcrse_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );
            
        }
        
        
        if ( level < lfine_m ) {
            
            /*
             * restriction matrix
             */
            
            // refinement 2
            int nNeighbours = AMREX_D_TERM(2, * 2, * 2);
            
            mglevel_m[level]->R_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );
            
            /*
             * fine boundary matrix
             */
            
            //                              refinement 2
            nNeighbours = 2 * AMREX_SPACEDIM * AMREX_D_TERM(2, * 2, * 2);
            
            mglevel_m[level]->Bfine_p = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nNeighbours,
                             Tpetra::StaticProfile) );
            
        }
        
        /*
         * no-fine Poisson matrix
         */
        
        int nPhysBoundary = 2 * AMREX_SPACEDIM * mglevel_m[level]->getBCStencilNum();
    
        // number of internal stencil points
        int nIntBoundary = AMREX_SPACEDIM * interface_mp->getNumberOfPoints();
    
        int nEntries = (AMREX_SPACEDIM << 1) + 1 /* plus boundaries */ + nPhysBoundary + nIntBoundary;
    
        mglevel_m[level]->Anf_p = Teuchos::rcp(
            new matrix_t(mglevel_m[level]->map_p, nEntries,
                         Tpetra::StaticProfile) );
        
        /*
         * with-fine / composite Poisson matrix
         */
        
        nEntries = (AMREX_SPACEDIM << 1) + 5 /* plus boundaries */ +
                   nPhysBoundary + nIntBoundary;
    
        mglevel_m[level]->Awf_p = Teuchos::rcp(
            new matrix_t(mglevel_m[level]->map_p, nEntries,
                         Tpetra::StaticProfile) );
        
        /*
         * uncovered cells matrix
         */
        mglevel_m[level]->UnCovered_p = Teuchos::rcp(
            new matrix_t(mglevel_m[level]->map_p, 1,
                         Tpetra::StaticProfile) );
        
        /*
         * gradient matrices
         */
        nEntries = 5;
    
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            mglevel_m[level]->G_p[d] = Teuchos::rcp(
                new matrix_t(mglevel_m[level]->map_p, nEntries,
                             Tpetra::StaticProfile) );
        }
    }
    
    
    mglevel_m[level]->rho_p = Teuchos::rcp(
        new vector_t(mglevel_m[level]->map_p, false) );
    
    mglevel_m[level]->phi_p = Teuchos::rcp(
        new vector_t(mglevel_m[level]->map_p, false) );
}


void AmrMultiGrid::close_m(int level, bool matrices) {
    
    if ( matrices ) {
        if ( level > lbase_m ) {
            
            mglevel_m[level]->I_p->fillComplete(mglevel_m[level-1]->map_p,  // col map (domain map)
                                                mglevel_m[level]->map_p);   // row map (range map)
            
            mglevel_m[level]->Bcrse_p->fillComplete(mglevel_m[level-1]->map_p,  // col map
                                                    mglevel_m[level]->map_p);   // row map
        }
        
        if ( level < lfine_m ) {
            
            mglevel_m[level]->R_p->fillComplete(mglevel_m[level+1]->map_p,
                                                  mglevel_m[level]->map_p);
            
            mglevel_m[level]->Bfine_p->fillComplete(mglevel_m[level+1]->map_p,
                                                    mglevel_m[level]->map_p);
        }
        
        mglevel_m[level]->Anf_p->fillComplete();
        
        mglevel_m[level]->Awf_p->fillComplete();
        
        mglevel_m[level]->UnCovered_p->fillComplete();
        
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            mglevel_m[level]->G_p[d]->fillComplete();
    }
}


void AmrMultiGrid::buildNoFinePoissonMatrix_m(const AmrIntVect_t& iv,
                                              const basefab_t& mfab, int level)
{
    /*
     * Laplacian of "no fine"
     */
    
    /*
     * 1D not supported
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    amrex::BoxArray empty;
    
    indices_t indices;
    coefficients_t values;
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    const double idx2[] = {
        D_DECL( 1.0 / ( dx[0] * dx[0] ),
                1.0 / ( dx[1] * dx[1] ),
                1.0 / ( dx[2] * dx[2] ) )
    };
                    
    int globalidx = mglevel_m[level]->serialize(iv);
    
    /*
     * check neighbours in all directions (Laplacian stencil --> cross)
     */
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
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
                    values.push_back( idx2[d] );
                    break;
                }
                case AmrMultiGridLevel_t::Mask::BNDRY:
                {
                    // boundary cell
                    // only level > 0 have this kind of boundary
#if DEBUG
                    if ( level == lbase_m )
                        throw std::runtime_error("Error in mask for level "
                                                 + std::to_string(level) + "!");
#endif
                    
                    
                    /* Dirichlet boundary conditions from coarser level.
                     */
                    std::size_t nn = indices.size();
                    
                    interface_mp->fine(biv, indices, values, d, -shift, empty,
                                       mglevel_m[level].get());
                    
                    double value = idx2[d];
                    for (std::size_t iter = nn; iter < indices.size(); ++iter)
                        values[iter] *= value;
                    
                    break;
                }
                case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                {
                    // physical boundary cell
                    double value = idx2[d];
                    mglevel_m[level]->applyBoundary(biv,
                                                    indices,
                                                    values,
                                                    value /*matrix coefficient*/);
                    break;
                }
                default:
                    throw std::runtime_error("Error in mask for level "
                                             + std::to_string(level) + "!");
            }
        }
    }
    
    // check center
    indices.push_back( globalidx );
    values.push_back( AMREX_D_TERM(- 2.0 * idx2[0],
                                   - 2.0 * idx2[1],
                                   - 2.0 * idx2[2]) );
    
    this->unique_m(indices, values);
    
    mglevel_m[level]->Anf_p->insertGlobalValues(globalidx,
                                                indices.size(),
                                                &values[0],
                                                &indices[0]);
}


void AmrMultiGrid::buildCompositePoissonMatrix_m(const AmrIntVect_t& iv,
                                                 const basefab_t& mfab,
                                                 const boxarray_t& crse_fine_ba,
                                                 int level)
{
    /*
     * Laplacian of "with fine"
     * 
     * For the finest level: Awf == Anf
     */
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    
    const double idx2[] = {
        D_DECL( 1.0 / ( cdx[0] * cdx[0] ),
                1.0 / ( cdx[1] * cdx[1] ),
                1.0 / ( cdx[2] * cdx[2] ) )
    };
    
    
    const double* fdx = 0;
    
    // finest level --> crse_fine_ba is an empty box array
    
    if ( level < lfine_m )
        fdx = mglevel_m[level+1]->geom.CellSize();
    
    /*
     * 1D not supported by AmrMultiGrid
     * 2D --> 5 elements per row
     * 3D --> 7 elements per row
     */
    
    indices_t indices;
    coefficients_t values;
    
    if ( !crse_fine_ba.contains(iv) ) {
        /*
         * Only cells that are not refined
         */
        int globalidx = mglevel_m[level]->serialize(iv);
        
        /*
         * check neighbours in all directions (Laplacian stencil --> cross)
         */
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
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
                            values.push_back( idx2[d] );
                            
                            // add center once
                            indices.push_back( globalidx );
                            values.push_back( -idx2[d] );
                            
                            break;
                        }
                        case AmrMultiGridLevel_t::Mask::BNDRY:
                        {
                            // boundary cell
                            // only level > 0 have this kind of boundary
#if DEBUG
                            if ( level == lbase_m )
                                throw std::runtime_error("Error in mask for level "
                                                         + std::to_string(level) + "!");
#endif
                                            
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
                            
                            double value = idx2[d];
                            for (std::size_t iter = nn; iter < indices.size(); ++iter)
                                values[iter] *= value;
                            
                            // add center once
                            indices.push_back( globalidx );
                            values.push_back( -idx2[d] );
                            
                            break;
                        }
                        case AmrMultiGridLevel_t::Mask::PHYSBNDRY:
                        {
                            // physical boundary cell
                            double value = idx2[d];
                            
                            mglevel_m[level]->applyBoundary(biv,
                                                            indices,
                                                            values,
                                                            value /*matrix coefficient*/);
                            
                            // add center once
                            indices.push_back( globalidx );
                            values.push_back( -idx2[d] );
                            
                            break;
                        }
                        default:
                            throw std::runtime_error("Error in mask for level "
                                                     + std::to_string(level) + "!");
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
#if AMREX_SPACEDIM == 3
                        for (int d2 = 0; d2 < 2; ++d2) {
#endif
                            
                            /* in order to get a top iv --> needs to be odd value in "d"
                             * in order to get a bottom iv --> needs to be even value in "d"
                             */
                            AmrIntVect_t fake(D_DECL(0, 0, 0));
                            
                            fake[(d+1)%AMREX_SPACEDIM] = d1;
#if AMREX_SPACEDIM == 3
                            fake[(d+2)%AMREX_SPACEDIM] = d2;
#endif
                            interface_mp->coarse(iv, indices, values, d, shift, crse_fine_ba,
                                                 fake, mglevel_m[level].get());
                            
#if AMREX_SPACEDIM == 3
                        }
#endif
                    }
                    
                    double avg = AMREX_D_PICK(1.0, 2.0, 4.0);
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
}


void AmrMultiGrid::buildRestrictionMatrix_m(const AmrIntVect_t& iv,
                                            const boxarray_t& crse_fine_ba,
                                            D_DECL(int& ii, int& jj, int& kk),
                                            int level)
{
    /*
     * x^(l) = R * x^(l+1)
     */
    
    // finest level does not need to have a restriction matrix
    if ( level == lfine_m )
        return;
    
    /* Difficulty:  If a fine cell belongs to another processor than the underlying
     *              coarse cell, we get an error when filling the matrix since the
     *              cell (--> global index) does not belong to the same processor.
     * Solution:    Find all coarse cells that are covered by fine cells, thus,
     *              the distributionmap is correct.
     * 
     * 
     */
    
    indices_t indices;
    coefficients_t values;
    
    if ( crse_fine_ba.contains(iv) ) {
        int crse_globalidx = mglevel_m[level]->serialize(iv);
        
        // neighbours
        for (int iref = 0; iref < 2; ++iref) {
            for (int jref = 0; jref < 2; ++jref) {
#if AMREX_SPACEDIM == 3
                for (int kref = 0; kref < 2; ++kref) {
#endif
                    AmrIntVect_t riv(D_DECL(ii + iref, jj + jref, kk + kref));
                                    
                    int fine_globalidx = mglevel_m[level+1]->serialize(riv);
                                    
                    indices.push_back( fine_globalidx );
                    values.push_back( AMREX_D_PICK(0.5, 0.25, 0.125) );
#if AMREX_SPACEDIM == 3
                }
#endif
            }
        }
        
        mglevel_m[level]->R_p->insertGlobalValues(crse_globalidx,
                                                  indices.size(),
                                                  &values[0],
                                                  &indices[0]);
    }
}


void AmrMultiGrid::buildInterpolationMatrix_m(const AmrIntVect_t& iv,
                                              int level)
{
    /* crse: level - 1
     * fine (this): level
     */
    
    /*
     * This does not include ghost cells
     * --> no boundaries
     * 
     * x^(l) = I * x^(l-1)
     */
    
    if ( level == lbase_m )
        return;
    
    indices_t indices;
    coefficients_t values;
                    
    int globalidx = mglevel_m[level]->serialize(iv);
                    
    /*
     * we need boundary + indices from coarser level
     */
    interp_mp->stencil(iv, indices, values, mglevel_m[level-1].get());
    
    this->unique_m(indices, values);
    
    mglevel_m[level]->I_p->insertGlobalValues(globalidx,
                                              indices.size(),
                                              &values[0],
                                              &indices[0]);
}


void AmrMultiGrid::buildCrseBoundaryMatrix_m(const AmrIntVect_t& iv,
                                             map_t& cells, int level)
{
    /*
     * fine (this): level
     * coarse:      level - 1
     */
    
    // the base level has only physical boundaries
    if ( level == lbase_m )
        return;
    
    // iv is a fine cell
    
    // check its neighbours to see if at crse-fine interface
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (int shift = -1; shift <= 1; shift += 2) {
            // neighbour
            AmrIntVect_t niv = iv;
            niv[d] += shift;
            
            if ( !mglevel_m[level]->grids.contains(niv) ) {
                // neighbour does not belong to fine grids
                
                // insert iv if not yet exists
                cells[iv].push_back(std::make_pair(shift, d));
            }
        }
    }
}


void AmrMultiGrid::fillCrseBoundaryMatrix_m(map_t& cells,
                                            const boxarray_t& crse_fine_ba,
                                            int level)
{
    // the base level has only physical boundaries
    if ( level == lbase_m )
        return;
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    const double idx2[] = {
        D_DECL( 1.0 / ( dx[0] * dx[0] ),
                1.0 / ( dx[1] * dx[1] ),
                1.0 / ( dx[2] * dx[2] ) )
    };
    
    indices_t indices;
    coefficients_t values;
    
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
                values[n] *= idx2[d];
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
}


void AmrMultiGrid::buildFineBoundaryMatrix_m(const AmrIntVect_t& iv,
                                             map_t& cells,
                                             const boxarray_t& crse_fine_ba,
                                             int level)
{
    /* fine: level + 1
     * coarse (this): level
     */
    
    // the finest level does not need data from a finer level
    if ( level == lfine_m )
        return;
    
    if ( !crse_fine_ba.contains(iv) /*&& mfab(iv) != 2 not physical bc*/ ) {
        /*
         * iv is a coarse cell that got not refined
         * 
         * --> check all neighbours to see if at crse-fine
         * interface
         */
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
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
}
    
    
void AmrMultiGrid::fillFineBoundaryMatrix_m(map_t& cells,
                                            const boxarray_t& crse_fine_ba,
                                            int level)
{
    // the finest level does not need data from a finer level
    if ( level == lfine_m )
        return;
    
    indices_t indices;
    coefficients_t values;
    
    const double* cdx = mglevel_m[level]->geom.CellSize();
    const double* fdx = mglevel_m[level+1]->geom.CellSize();
    
    auto fill = [&](indices_t& indices,
                    coefficients_t& values,
                    D_DECL(int ii, int jj, int kk),
                    int* begin, int* end, int d,
                    const AmrIntVect_t& iv, int shift,
                    int sign)
    {
        // number of fine cell gradients
        double avg = AMREX_D_PICK(1, 2, 4);
        
        for (int iref = ii - begin[0]; iref <= ii + end[0]; ++iref) {
            for (int jref = jj - begin[1]; jref <= jj + end[1]; ++jref) {
#if AMREX_SPACEDIM == 3
                for (int kref = kk - begin[2]; kref <= kk + end[2]; ++kref) {
#endif
                    /* Since all fine cells on the not-refined cell are
                     * outside of the "domain" --> we need to interpolate
                     */
                    AmrIntVect_t riv(D_DECL(iref, jref, kref));
                    
                    if ( riv[d] / 2 /*refinement*/ == iv[d] ) {
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
#if AMREX_SPACEDIM == 3
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
            int begin[AMREX_SPACEDIM] = { D_DECL( int(d == 0), int(d == 1), int(d == 2) ) };
            int end[AMREX_SPACEDIM]   = { D_DECL( int(d != 0), int(d != 1), int(d != 2) ) };
            
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
                    int ii = iv[0] * 2; // refinemet in x
                    int jj = iv[1] * 2; // refinemet in y
#if AMREX_SPACEDIM == 3
                    int kk = iv[2] * 2; // refinemet in z
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
                    int ii = covered[0] * 2; // refinemet in x
                    int jj = covered[1] * 2; // refinemet in y
#if AMREX_SPACEDIM == 3
                    int kk = covered[2] * 2; // refinemet in z
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
}


inline void AmrMultiGrid::buildDensityVector_m(const AmrField_t& rho, int level) {
    this->amrex2trilinos_m(rho, 0, mglevel_m[level]->rho_p, level);
}


inline void AmrMultiGrid::buildPotentialVector_m(const AmrField_t& phi, int level) {
    this->amrex2trilinos_m(phi, 0, mglevel_m[level]->phi_p, level);
}


void AmrMultiGrid::buildGradientMatrix_m(const AmrIntVect_t& iv,
                                         const basefab_t& mfab,
                                         int level)
{
    
    const double* dx = mglevel_m[level]->geom.CellSize();
    
    indices_t indices;
    coefficients_t values;
    
    auto check = [&](const AmrIntVect_t& iv,
                     const basefab_t& mfab,
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
#if DEBUG
                if ( level == lbase_m )
                    throw std::runtime_error("Error in mask for level "
                                             + std::to_string(level) + "!");
#endif
                
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
    
    
    int globalidx = mglevel_m[level]->serialize(iv);
    
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        
        indices.clear();
        values.clear();
        
        for (int shift = -1; shift <= 1; shift += 2) {
            AmrIntVect_t niv = iv;                        
            niv[d] += shift;
            check(niv, mfab, d, shift);
        }
        
        this->unique_m(indices, values);
        
        mglevel_m[level]->G_p[d]->insertGlobalValues(globalidx,
                                                     indices.size(),
                                                     &values[0],
                                                     &indices[0]);
    }
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
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    
                    int globalidx = mglevel_m[level]->serialize(iv);
                    
                    mv->replaceGlobalValue(globalidx, fab(iv, comp));
#if AMREX_SPACEDIM == 3
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
#if AMREX_SPACEDIM == 3
                for (int k = lo[2]; k <= hi[2]; ++k) {
#endif
                    AmrIntVect_t iv(D_DECL(i, j, k));
                    fab(iv, comp) = data[localidx++];
                }
#if AMREX_SPACEDIM == 3
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


void AmrMultiGrid::smooth_m(Teuchos::RCP<vector_t>& e,
                            Teuchos::RCP<vector_t>& r,
                            int level)
{
#if AMR_MG_TIMER
        IpplTimings::startTimer(smoothTimer_m);
#endif

    // apply "no fine" Laplacian
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
    
    // base level has no smoother --> l - 1
    smoother_m[level-1]->smooth(e, mglevel_m[level]->Anf_p, r);
    
#if AMR_MG_TIMER
        IpplTimings::stopTimer(smoothTimer_m);
#endif
}


void AmrMultiGrid::restrict_m(int level) {
    
#if AMR_MG_TIMER
        IpplTimings::startTimer(restrictTimer_m);
#endif
    
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
    
    this->residual_no_fine_m(tmp,
                             mglevel_m[level]->error_p,
                             mglevel_m[level-1]->error_p,
                             mglevel_m[level]->residual_p, level);
    
    mglevel_m[level-1]->residual_p->putScalar(0.0);
    
    // average down: residual^(l-1) = R^(l) * tmp
    mglevel_m[level-1]->R_p->apply(*tmp, *mglevel_m[level-1]->residual_p);
    
    
    // composite matrix, i.e. matrix without covered cells
    // r^(l-1) = rho^(l-1) - A * phi^(l-1)
    
    vector_t fine2crse(mglevel_m[level-1]->Awf_p->getDomainMap());
    
    // get boundary for 
    mglevel_m[level-1]->Bfine_p->apply(*mglevel_m[level]->phi_p, fine2crse);
    
    vector_t tmp2(mglevel_m[level-1]->Awf_p->getDomainMap());
    mglevel_m[level-1]->Awf_p->apply(*mglevel_m[level-1]->phi_p, tmp2);
    
    vector_t crse2fine(mglevel_m[level-1]->Awf_p->getDomainMap());
    
    if ( mglevel_m[level-1]->Bcrse_p != Teuchos::null ) {
        mglevel_m[level-1]->Bcrse_p->apply(*mglevel_m[level-2]->phi_p, crse2fine);
    }
    
    tmp2.update(1.0, fine2crse, 1.0, crse2fine, 1.0);
    
    Teuchos::RCP<vector_t> tmp3 = Teuchos::rcp( new vector_t(mglevel_m[level-1]->map_p, true) );
    
    mglevel_m[level-1]->UnCovered_p->apply(*mglevel_m[level-1]->rho_p, *tmp3);
    
    Teuchos::RCP<vector_t> tmp4 = Teuchos::rcp( new vector_t(mglevel_m[level-1]->map_p, true) );
    mglevel_m[level-1]->UnCovered_p->apply(tmp2, *tmp4);
    
    // ONLY subtract coarse rho
    mglevel_m[level-1]->residual_p->update(1.0, *tmp3, -1.0, *tmp4, 1.0);
    
#if AMR_MG_TIMER
        IpplTimings::stopTimer(restrictTimer_m);
#endif
}


void AmrMultiGrid::prolongate_m(int level) {
#if AMR_MG_TIMER
        IpplTimings::startTimer(interpTimer_m);
#endif
        // interpolate error from l-1 to l
        Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p) );
        mglevel_m[level]->I_p->apply(*mglevel_m[level-1]->error_p, *tmp);
        
        // e^(l) += tmp
        mglevel_m[level]->error_p->update(1.0, *tmp, 1.0);
        
#if AMR_MG_TIMER
        IpplTimings::stopTimer(interpTimer_m);
#endif
}


void AmrMultiGrid::averageDown_m(int level) {
    
    if (level == lfine_m )
        return;
    
    Teuchos::RCP<vector_t> tmp = Teuchos::rcp( new vector_t(mglevel_m[level]->map_p, false) );
    
    mglevel_m[level]->R_p->apply(*mglevel_m[level+1]->phi_p, *tmp);
    
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


void AmrMultiGrid::initBaseSolver_m(const BaseSolver& solver,
                                    const Preconditioner& precond)
{
    switch ( solver ) {
        case BaseSolver::BICGSTAB:
            solver_mp.reset( new BelosBottomSolver("BICGSTAB", precond) );
            break;
        case BaseSolver::MINRES:
            solver_mp.reset( new BelosBottomSolver("MINRES", precond) );
            break;
        case BaseSolver::PCPG:
            solver_mp.reset( new BelosBottomSolver("PCPG", precond) );
            break;
        case BaseSolver::CG:
            solver_mp.reset( new BelosBottomSolver("Pseudoblock CG", precond) );
            break;
        case BaseSolver::GMRES:
            solver_mp.reset( new BelosBottomSolver("Pseudoblock GMRES", precond) );
            break;
        case BaseSolver::STOCHASTIC_CG:
            solver_mp.reset( new BelosBottomSolver("Stochastic CG", precond) );
            break;
        case BaseSolver::RECYCLING_CG:
            solver_mp.reset( new BelosBottomSolver("RCG", precond) );
            break;
        case BaseSolver::RECYCLING_GMRES:
            solver_mp.reset( new BelosBottomSolver("GCRODR", precond) );
            break;
        default:
            std::runtime_error("No such bottom solver available.");
    }
}
