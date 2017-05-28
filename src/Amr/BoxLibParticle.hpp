#ifndef BOXLIB_PARTICLE_HPP
#define BOXLIB_PARTICLE_HPP

template<class PLayout>
BoxLibParticle<PLayout>::BoxLibParticle() : AmrParticleBase<PLayout>()
{
    AssignDensityTimer_m = IpplTimings::getTimer("AMR AssignDensity");
}


template<class PLayout>
BoxLibParticle<PLayout>::BoxLibParticle(PLayout *layout) : AmrParticleBase<PLayout>(layout)
{
    AssignDensityTimer_m = IpplTimings::getTimer("AMR AssignDensity");
}


template<class PLayout>
template<class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::scatter(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                                      ParticleAttrib<Vektor<PT, Dim> >& pp,
                                      int lbase, int lfine)
{
    const PLayout *layout_p = &this->getLayout();
    int nGrow = layout_p->refRatio(lbase)[0];
    
    AmrFieldContainer_t tmp(lfine+1);
    for (int lev = lbase; lev <= lfine; ++lev) {
        const AmrGrid_t& ba = f[lev]->boxArray();
        const AmrProcMap_t& dm = f[lev]->DistributionMap();
        tmp[lev].reset(new AmrField_t(ba, dm, 1, nGrow));
        tmp[lev]->setVal(0.0);
    }
    
    this->AssignDensityFort(attrib, tmp, lbase, 1, lfine);
    
    for (int lev = lbase; lev <= lfine; ++lev)
        AmrField_t::Copy(*f[lev], *tmp[lev], 0, 0, f[lev]->nComp(), f[lev]->nGrow());
}


template<class PLayout>
template <class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::scatter(ParticleAttrib<FT>& attrib, AmrField_t& f,
                                      ParticleAttrib<Vektor<PT, Dim> >& pp,
                                      int level)
{
    const AmrGrid_t& ba      = f.boxArray();
    const AmrProcMap_t& dmap = f.DistributionMap();
    
    AmrField_t tmp(ba, dmap, f.nComp(), 1);
    tmp.setVal(0.0);
    
    this->AssignCellDensitySingleLevelFort(attrib, tmp, level);
    
    AmrField_t::Copy(f, tmp, 0, 0, f.nComp(), f.nGrow());
}


template<class PLayout>
template<class FT, unsigned Dim, class PT>
void BoxLibParticle<PLayout>::gather(ParticleAttrib<FT>& attrib, AmrFieldContainer_t& f,
                                     ParticleAttrib<Vektor<PT, Dim> >& pp,
                                     int lbase, int lfine)
{
    this->InterpolateFort(attrib, f, lbase, lfine);
}


// Function from BoxLib adjusted to work with Ippl BoxLibParticle class
// Scatter the particle attribute pa on the grid
template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::AssignDensityFort(ParticleAttrib<AType> &pa,
                                                AmrFieldContainer_t& mf_to_be_filled, 
                                                int lev_min, int ncomp, int finest_level) const
{
//     BL_PROFILE("AssignDensityFort()");
    IpplTimings::startTimer(AssignDensityTimer_m);
    const PLayout *layout_p = &this->getLayout();
    
    // not done in amrex
    int rho_index = 0;
    
    amrex::PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    amrex::Array<amrex::BCRec> bcs(1, amrex::BCRec(lo_bc, hi_bc));
    amrex::CellConservativeLinear mapper;
    
    AmrFieldContainer_t tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const AmrGrid_t& ba = mf_to_be_filled[lev]->boxArray();
        const AmrProcMap_t& dm = mf_to_be_filled[lev]->DistributionMap();
        tmp[lev].reset(new AmrField_t(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }
    
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        AssignCellDensitySingleLevelFort(pa, *mf_to_be_filled[lev], lev, 1, 0);

        if (lev < finest_level) {
            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf_to_be_filled[lev],
                                          rho_index, rho_index, ncomp, 
                                          layout_p->Geom(lev), layout_p->Geom(lev+1),
                                          cphysbc, fphysbc,
                                          layout_p->refRatio(lev), &mapper, bcs);
        }

        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(*mf_to_be_filled[lev],
                                      *mf_to_be_filled[lev-1],
                                      rho_index, 1, layout_p->refRatio(lev-1),
                                      layout_p->Geom(lev-1), layout_p->Geom(lev));
        }
        
        mf_to_be_filled[lev]->plus(*tmp[lev], rho_index, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf_to_be_filled[lev+1], 
                             *mf_to_be_filled[lev], rho_index, ncomp, layout_p->refRatio(lev));
    }
    
    IpplTimings::stopTimer(AssignDensityTimer_m);
}


// This is the single-level version for cell-centered density
template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::AssignCellDensitySingleLevelFort(ParticleAttrib<AType> &pa,
                                                               AmrField_t& mf_to_be_filled,
                                                               int       level,
                                                               int       ncomp,
                                                               int       particle_lvl_offset) const
{
//     BL_PROFILE("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignCellDensitySingleLevelFort()");
    
    const PLayout *layout_p = &this->getLayout();
    
    int rho_index = 0;
    
//     if (rho_index != 0) amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0");

    AmrField_t* mf_pointer;

    if (layout_p->OnSameGrids(level, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new AmrField_t(layout_p->ParticleBoxArray(level), 
                                  layout_p->ParticleDistributionMap(level),
                                  ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread 
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignCellDensitySingleLevelFort");

    const AmrReal_t      strttime    = amrex::ParallelDescriptor::second();
    const AmrGeometry_t& gm          = layout_p->Geom(level);
    const AmrReal_t*     plo         = gm.ProbLo();
    const AmrReal_t*     dx_particle = layout_p->Geom(level + particle_lvl_offset).CellSize();
    const AmrReal_t*     dx          = gm.CellSize();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }
    
    for (amrex::MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }
    
    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    
    //while lev_min > this->Level[start_idx] we need to skip these particles since there level is
    //higher than the specified lev_min
    int start_idx = 0;
    while ((unsigned)level > this->Level[start_idx])
        start_idx++;
    
    AmrReal_t inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    for (size_t ip = start_idx; ip < LocalNum; ++ip) {
        //if particle doesn't belong on this level exit loop
        if (this->Level[ip] != (unsigned)level)
            break;
        
        const int grid = this->Grid[ip];
        FArrayBox_t& fab = (*mf_pointer)[grid];
        const AmrBox_t& box = fab.box();
        
        // not callable:
        // begin amrex_deposit_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        AmrIntVect_t i1(i-1, j-1, k-1);
        AmrIntVect_t i2(i-1, j-1, k);
        AmrIntVect_t i3(i-1, j,   k-1);
        AmrIntVect_t i4(i-1, j,   k);
        AmrIntVect_t i5(i,   j-1, k-1);
        AmrIntVect_t i6(i,   j-1, k);
        AmrIntVect_t i7(i,   j,   k-1);
        AmrIntVect_t i8(i,   j,   k);
        
        fab(i1, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i2, 0) += wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i3, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i4, 0) += wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        fab(i5, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*pa[ip];
        fab(i6, 0) += wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*pa[ip];
        fab(i7, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*pa[ip];
        fab(i8, 0) += wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*pa[ip];
        // end of amrex_deposit_cic
    }
    
    mf_pointer->SumBoundary(gm.periodicity());

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const AmrReal_t vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }
}


template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::InterpolateFort(ParticleAttrib<AType> &pa,
                                              AmrFieldContainer_t& mesh_data, 
                                              int lev_min, int lev_max)
{
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        InterpolateSingleLevelFort(pa, *mesh_data[lev], lev); 
    }
}


template<class PLayout>
template <class AType>
void BoxLibParticle<PLayout>::InterpolateSingleLevelFort(ParticleAttrib<AType> &pa,
                                                         AmrField_t& mesh_data, int lev)
{
    if (mesh_data.nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in InterpolateSingleLevelFort");
    
    PLayout *layout_p = &this->getLayout();
    
    const AmrGeometry_t& gm          = layout_p->Geom(lev);
    const AmrReal_t*     plo         = gm.ProbLo();
    const AmrReal_t*     dx          = gm.CellSize();
    
    //loop trough particles and distribute values on the grid
    size_t LocalNum = this->getLocalNum();
    
    //while lev_min > this->Level[start_idx] we need to skip these particles since there level is
    //higher than the specified lev_min
    int start_idx = 0;
    while ((unsigned)lev > this->Level[start_idx])
        start_idx++;
    
    AmrReal_t inv_dx[3] = { 1.0 / dx[0], 1.0 / dx[1], 1.0 / dx[2] };
    double lxyz[3] = { 0.0, 0.0, 0.0 };
    double wxyz_hi[3] = { 0.0, 0.0, 0.0 };
    double wxyz_lo[3] = { 0.0, 0.0, 0.0 };
    int ijk[3] = {0, 0, 0};
    for (size_t ip = start_idx; ip < LocalNum; ++ip) {
        //if particle doesn't belong on this level exit loop
        if (this->Level[ip] != (unsigned)lev)
            break;
        
        const int grid = this->Grid[ip];
        FArrayBox_t& fab = mesh_data[grid];
        const AmrBox_t& box = fab.box();
        int nComp = fab.nComp();
        
        // not callable
        // begin amrex_interpolate_cic(pbx.data(), nstride, N, fab.dataPtr(), box.loVect(), box.hiVect(), nComp, plo, dx);
        for (int i = 0; i < 3; ++i) {
            lxyz[i] = ( this->R[ip](i) - plo[i] ) * inv_dx[i] + 0.5;
            ijk[i] = lxyz[i];
            wxyz_hi[i] = lxyz[i] - ijk[i];
            wxyz_lo[i] = 1.0 - wxyz_hi[i];
        }
        
        int& i = ijk[0];
        int& j = ijk[1];
        int& k = ijk[2];
        
        AmrIntVect_t i1(i-1, j-1, k-1);
        AmrIntVect_t i2(i-1, j-1, k);
        AmrIntVect_t i3(i-1, j,   k-1);
        AmrIntVect_t i4(i-1, j,   k);
        AmrIntVect_t i5(i,   j-1, k-1);
        AmrIntVect_t i6(i,   j-1, k);
        AmrIntVect_t i7(i,   j,   k-1);
        AmrIntVect_t i8(i,   j,   k);
        
        for (int nc = 0; nc < nComp; ++nc) {
            pa[ip](nc) = wxyz_lo[0]*wxyz_lo[1]*wxyz_lo[2]*fab(i1, nc) +
                         wxyz_lo[0]*wxyz_lo[1]*wxyz_hi[2]*fab(i2, nc) +
                         wxyz_lo[0]*wxyz_hi[1]*wxyz_lo[2]*fab(i3, nc) +
                         wxyz_lo[0]*wxyz_hi[1]*wxyz_hi[2]*fab(i4, nc) +
                         wxyz_hi[0]*wxyz_lo[1]*wxyz_lo[2]*fab(i5, nc) +
                         wxyz_hi[0]*wxyz_lo[1]*wxyz_hi[2]*fab(i6, nc) +
                         wxyz_hi[0]*wxyz_hi[1]*wxyz_lo[2]*fab(i7, nc) +
                         wxyz_hi[0]*wxyz_hi[1]*wxyz_hi[2]*fab(i8, nc);
        }
        // end amrex_interpolate_cic
    }
}


#endif
