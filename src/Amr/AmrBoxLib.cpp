#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
// #include "Solvers/AmrPoissonSolver.h"
#include "Solvers/PoissonSolver.h"

#include "AmrBoxLib_F.h"
#include <MultiFabUtil.H>

#include <ParmParse.H> // used in initialize function

AmrBoxLib::AmrBoxLib() : AmrObject(),
                         AmrCore(),
                         nChargePerCell_m(0),
                         bunch_mp(nullptr),
                         rho_m(0),
                         phi_m(0),
                         eg_m(0)
{}


AmrBoxLib::AmrBoxLib(TaggingCriteria tagging,
                     double scaling,
                     double nCharge)
    : AmrObject(tagging, scaling, nCharge),
      AmrCore(),
      nChargePerCell_m(0),
      bunch_mp(nullptr),
      rho_m(0),
      phi_m(0),
      eg_m(0)
{}


// AmrBoxLib::AmrBoxLib(const DomainBoundary_t& realbox,
//               const NDIndex<3>& nGridPts,
//               short maxLevel,
//               const RefineRatios_t& refRatio)
//     : AmrObject(realbox, nGridPts, maxLevel, refRatio)
// {
//     
//     
// }



AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel,
                     const AmrIntArray_t& refRatio)
    : AmrObject(),
      AmrCore(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      nChargePerCell_m(maxLevel + 1),
      bunch_mp(nullptr),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      eg_m(maxLevel + 1)
{}


AmrBoxLib::AmrBoxLib(const AmrDomain_t& domain,
                     const AmrIntArray_t& nGridPts,
                     short maxLevel,
                     const AmrIntArray_t& refRatio,
                     AmrPartBunch* bunch)
    : AmrObject(),
      AmrCore(&domain, maxLevel, nGridPts, 0 /* cartesian */),
      nChargePerCell_m(maxLevel + 1),
      bunch_mp(bunch),
      rho_m(maxLevel + 1),
      phi_m(maxLevel + 1),
      eg_m(maxLevel + 1)
{
    initBaseLevel_m();
}


void AmrBoxLib::setBunch(AmrPartBunch* bunch) {
    bunch_mp = bunch;
}


// void AmrBoxLib::initialize(const DomainBoundary_t& lower,
//                            const DomainBoundary_t& upper,
//                            const NDIndex<3>& nGridPts,
//                            short maxLevel,
//                            const RefineRatios_t& refRatio)
// {
//     // setup the physical domain
//     RealBox domain;
//     for (int i = 0; i < 3; ++i) {
//         domain.setLo(i, lower[i]); // m
//         domain.setHi(i, upper[i]); // m
//     }
//     
//     Array<int> nCells(3);
//     for (int i = 0; i < 3; ++i)
//         nCells = nGridPts[i];
//     
//     AmrCore::Initialize();
//     Geometry::Setup(&domain, 0 /* cartesian coordinates */);
//     AmrCore::InitAmrCore(maxLevel, nCells);
// }


void AmrBoxLib::regrid(int lbase, int lfine, double time) {
    int new_finest = 0;
    AmrGridContainer_t new_grids(finest_level+2);
    
    MakeNewGrids(lbase, time, new_finest, new_grids);

    BL_ASSERT(new_finest <= finest_level+1);

    DistributionMapping::FlushCache();
    
    for (int lev = lbase+1; lev <= new_finest; ++lev)
    {
        if (lev <= finest_level) // an old level
        {
            if (new_grids[lev] != grids[lev]) // otherwise nothing
            {
                DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
                RemakeLevel(lev, time, new_grids[lev], new_dmap);
                
                /*
                 * particles need to know the BoxArray
                 * and DistributionMapping
                 */
                AmrLayout_t* layout_p = static_cast<AmrLayout_t*>(&bunch_mp->getLayout());
                layout_p->SetParticleBoxArray(lev, new_grids[lev]);
                layout_p->SetParticleDistributionMap(lev, new_dmap);
	    }
	}
	else  // a new level
	{
	    DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
	    MakeNewLevel(lev, time, new_grids[lev], new_dmap);
            
            /*
             * particles need to know the BoxArray
             * and DistributionMapping
             */
            AmrLayout_t* layout_p = static_cast<AmrLayout_t*>(&bunch_mp->getLayout());
            layout_p->SetParticleBoxArray(lev, new_grids[lev]);
            layout_p->SetParticleDistributionMap(lev, new_dmap);
        }
    }
    
    finest_level = new_finest;
    
//     // update to multilevel
//     bunch_m->update();
}


AmrBoxLib::VectorPair_t AmrBoxLib::getEExtrema() {
    Vector_t maxE = Vector_t(0, 0, 0), minE = Vector_t(0, 0, 0);
    /* check for maximum / minimum values over all levels
     * and all components, i.e. Ex, Ey, Ez
     */
    for (unsigned int lev = 0; lev < eg_m.size(); ++lev) {
        for (std::size_t i = 0; i < bunch_mp->Dimension; ++i) {
            /* calls:
             *  - max(comp, nghost (to search), local)
             *  - min(cmop, nghost, local)
             */
            double max = eg_m[lev]->max(i, 0, false);
            maxE[i] = (maxE[i] < max) ? max : maxE[i];
            
            double min = eg_m[lev]->min(i, 0, false);
            minE[i] = (minE[i] > min) ? min : minE[i];
        }
    }
    
    return VectorPair_t(maxE, minE);
}


double AmrBoxLib::getRho(int x, int y, int z) {
    //TODO
//     for (int lev = phi_m.size() - 2; lev >= 0; lev--)
//         BoxLib::average_down(phi_m[lev+1], phi_m[lev], 0, 1, this->refRatio(lev));
    throw OpalException("AmrBoxLib::getRho(x, y, z) ", "Not yet Implemented.");
    return 0.0;
}


void AmrBoxLib::computeSelfFields() {
    //FIXME Lorentz transformation
    //scatter charges onto grid
//     bunch_mp->Q *= bunch_mp->dt;
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, -1);
//     bunch_mp->Q /= bunch_mp->dt;
    
    int baseLevel = 0;
    int finestLevel = (&amrpbase_p->getAmrLayout())->finestLevel();
    
    int nLevel = finestLevel + 1;
    rho_m.resize(nLevel);
    phi_m.resize(nLevel);
    eg_m.resize(nLevel);
    
    double invDt = 1.0 / bunch_mp->getdT() * bunch_mp->getCouplingConstant();
    for (int i = 0; i <= finestLevel; ++i) {
        this->rho_m[i]->mult(invDt, 0, 1);
    }
    
    //calculating mesh-scale factor
    double gammaz = sum(bunch_mp->P)[2] / bunch_mp->getTotalNum();
    gammaz *= gammaz;
    gammaz = std::sqrt(gammaz + 1.0);
    
    // charge density is in rho_m
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, eg_m, baseLevel, finestLevel);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    amrpbase_p->gather(bunch_mp->Ef, this->eg_m, bunch_mp->R, 0, -1);
    
    /** Magnetic field in x and y direction induced by the eletric field
     *
     *  \f[ B_x = \gamma(B_x^{'} - \frac{beta}{c}E_y^{'}) = -\gamma \frac{beta}{c}E_y^{'} = -\frac{beta}{c}E_y \f]
     *  \f[ B_y = \gamma(B_y^{'} - \frac{beta}{c}E_x^{'}) = +\gamma \frac{beta}{c}E_x^{'} = +\frac{beta}{c}E_x \f]
     *  \f[ B_z = B_z^{'} = 0 \f]
     *
     */
    double betaC = sqrt(gammaz * gammaz - 1.0) / gammaz / Physics::c;

    bunch_mp->Bf(0) = bunch_mp->Bf(0) - betaC * bunch_mp->Ef(1);
    bunch_mp->Bf(1) = bunch_mp->Bf(1) + betaC * bunch_mp->Ef(0);
}


void AmrBoxLib::computeSelfFields(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::computeSelfFields_cycl(double gamma) {
    //FIXME Lorentz transformation
    //scatter charges onto grid
//     bunch_mp->Q *= bunch_mp->dt;
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, -1);
//     bunch_mp->Q /= bunch_mp->dt;
    int baseLevel = 0;
    int finestLevel = (&amrpbase_p->getAmrLayout())->finestLevel();
    int nLevel = finestLevel + 1;
    rho_m.resize(nLevel);
    phi_m.resize(nLevel);
    eg_m.resize(nLevel);
    double invDt = 1.0 / bunch_mp->getdT() * bunch_mp->getCouplingConstant();
    
    for (int i = 0; i <= finestLevel; ++i) {
        this->rho_m[i]->mult(invDt, 0, 1);
    }
    // charge density is in rho_m
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    IpplTimings::startTimer(bunch_mp->compPotenTimer_m);
    solver->solve(rho_m, phi_m, eg_m, baseLevel, finestLevel);
    IpplTimings::stopTimer(bunch_mp->compPotenTimer_m);
    
    amrpbase_p->gather(bunch_mp->Ef, this->eg_m, bunch_mp->R, 0, -1);
    /// calculate coefficient
    // Relativistic E&M says gamma*v/c^2 = gamma*beta/c = sqrt(gamma*gamma-1)/c
    // but because we already transformed E_trans into the moving frame we have to
    // add 1/gamma so we are using the E_trans from the rest frame -DW
    double betaC = sqrt(gamma * gamma - 1.0) / gamma / Physics::c;
    
    /// calculate B field from E field
    bunch_mp->Bf(0) =  betaC * bunch_mp->Ef(2);
    bunch_mp->Bf(2) = -betaC * bunch_mp->Ef(0);
    
//     std::cout << "self-field computed" << std::endl;
}


void AmrBoxLib::computeSelfFields_cycl(int bin) {
    //TODO
    throw OpalException("AmrBoxLib::computeSelfFields_cycl(int bin) ", "Not yet Implemented.");
}


void AmrBoxLib::updateMesh() {
    //FIXME What about resizing mesh, i.e. geometry?
    const AmrGeomContainer_t& geom = this->Geom();
    
    const Real* tmp = geom[0].CellSize();
    
    Vector_t hr;
    for (int i = 0; i < 3; ++i)
        hr[i] = tmp[i];
    
    bunch_mp->setBaseLevelMeshSpacing(hr);
}

Vektor<int, 3> AmrBoxLib::getBaseLevelGridPoints() {
    const AmrGeomContainer_t& geom = this->Geom();
    const Box& bx = geom[0].Domain();
    
    const IntVect& low = bx.smallEnd();
    const IntVect& high = bx.bigEnd();
    
    return Vektor<int, 3>(high[0] - low[0] + 1,
                          high[1] - low[1] + 1,
                          high[2] - low[2] + 1);
}


void AmrBoxLib::RemakeLevel (int lev, Real time,
                             const BoxArray& new_grids,
                             const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, 1, 1, new_dmap));
}


void AmrBoxLib::MakeNewLevel (int lev, Real time,
                              const BoxArray& new_grids,
                              const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m[lev].reset(new AmrField_t(new_grids, 1, 1, dmap[lev]));
}


void AmrBoxLib::ClearLevel(int lev) {
    nChargePerCell_m[lev].reset(nullptr);
    ClearBoxArray(lev);
    ClearDistributionMap(lev);
}


void AmrBoxLib::ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) {
    switch ( tagging_m ) {
        case CHARGE_DENSITY:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
        case POTENTIAL:
            tagForPotentialStrength_m(lev, tags, time, ngrow);
            break;
        case EFIELD:
            tagForEfield_m(lev, tags, time, ngrow);
            break;
        default:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
    }
}


void AmrBoxLib::tagForChargeDensity_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    for (int i = lev; i <= finest_level; ++i) {
        nChargePerCell_m[i]->setVal(0.0);
        // single level scatter
        amrpbase_p->scatter(bunch_mp->Q, *nChargePerCell_m[i], bunch_mp->R, i);
    }
    
    for (int i = finest_level-1; i >= lev; --i) {
        AmrField_t tmp(nChargePerCell_m[i]->boxArray(), 1, 0, nChargePerCell_m[i]->DistributionMap());
        tmp.setVal(0.0);
        BoxLib::average_down(*nChargePerCell_m[i+1], tmp, 0, 1, refRatio(i));
        AmrField_t::Add(*nChargePerCell_m[i], tmp, 0, 0, 1, 0);
    }
    
    
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        for (MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*nChargePerCell_m[lev])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &nCharge_m);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
}


void AmrBoxLib::tagForPotentialStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the potential is higher than 75 percent of the maximum potential
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[9] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev), 1, 1,
                                                  this->DistributionMap(lev)));
    nPartPerCell[baseLevel]->setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int lev = baseLevel; lev <= finestLevel; lev++)
            nPartPerCell[lev]->plus(+offset, 0, 1, 0);
    }
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi(1);
    AmrFieldContainer_t efield(1);
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, baseLevel, finestLevel);
    
    // 3. Tag cells for refinement
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    
    Real minp = 0.0;
    Real maxp = 0.0;
    maxp = scaling_m * phi[baseLevel]->max(0);
    minp = scaling_m * phi[baseLevel]->min(0);
    Real threshold = std::max( std::abs(minp), std::abs(maxp) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        for (MFIter mfi(*phi[baseLevel],false); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
            tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            tag_potential_strength(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                        BL_TO_FORTRAN_3D((*phi[baseLevel])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &threshold);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    phi[0].reset(nullptr);
    efield[0].reset(nullptr);
    nPartPerCell[0].reset(nullptr);
}


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int baseLevel   = 0;
    int finestLevel = 0;
    
    AmrFieldContainer_t nPartPerCell(1);
    nPartPerCell[baseLevel] = std::unique_ptr<AmrField_t>(new AmrField_t(this->boxArray(lev), 1, 1,
                                                          this->DistributionMap(lev)));
    nPartPerCell[baseLevel]->setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    // single level scatter
    amrpbase_p->scatter(bunch_mp->Q, *nPartPerCell[baseLevel], bunch_mp->R, lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int lev = baseLevel; lev <= finestLevel; lev++)
            nPartPerCell[lev]->plus(+offset, 0, 1, 0);
    }
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi(1);
    AmrFieldContainer_t efield(1);
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, baseLevel, finestLevel);
    
    // 3. Tag cells for refinement
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    
    Real min[3] = {0.0, 0.0, 0.0};
    Real max[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        max[i] = scaling_m * efield[baseLevel]->max(i);
        min[i] = scaling_m * efield[baseLevel]->min(i);
    }
    Real threshold[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        threshold[i] = std::max( std::abs(min[i]), std::abs(max[i]) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        // mfi(efield[baseLevel], true)
        for (MFIter mfi(*efield[baseLevel],false); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            FArrayBox&  fab     = (*efield[baseLevel])[mfi];
            // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
            // So we are going to get a temporary integer array.
//             tagfab.get_itags(itags, tilebx);
            
            // data pointer and index space
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tilebx.loVect();
            const int*  thi     = tilebx.hiVect();

            for (int i = tlo[0]; i <= thi[0]; ++i) {
                for (int j = tlo[1]; j <= thi[1]; ++j) {
                    for (int k = tlo[2]; k <= thi[2]; ++k) {
                        IntVect iv(i,j,k);
                        if (std::abs(fab(iv, 0)) >= threshold[0])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 1)) >= threshold[1])
                            tagfab(iv) = tagval;
                        else if (std::abs(fab(iv, 2)) >= threshold[2])
                            tagfab(iv) = tagval;
                        else
                            tagfab(iv) = clearval;
                    }
                }
            }
                            
                        
            //
            // Now update the tags in the TagBox.
            //
//             tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    phi[baseLevel].reset(nullptr);
    efield[baseLevel].reset(nullptr);
    nPartPerCell[baseLevel].reset(nullptr);
    
}


void AmrBoxLib::initBaseLevel_m() {
    AmrLayout_t* layout_p = static_cast<AmrLayout_t*>(&bunch_mp->getLayout());
    const BoxArray& ba = layout_p->ParticleBoxArray(0);
    const DistributionMapping& dm = layout_p->ParticleDistributionMap(0 /*level*/);
    
    SetBoxArray(0 /*level*/, ba);
    SetDistributionMap(0 /*level*/, dm);
    
    rho_m[0] = std::unique_ptr<AmrField_t>(new AmrField_t(ba, 1, 1, dm));
    rho_m[0]->setVal(0.0);
}

