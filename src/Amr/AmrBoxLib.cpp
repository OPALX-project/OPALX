#include "AmrBoxLib.h"

#include "Algorithms/AmrPartBunch.h"
// #include "Solvers/AmrPoissonSolver.h"
#include "Solvers/PoissonSolver.h"

#include "AmrBoxLib_F.h"
#include <MultiFabUtil.H>

#include <ParmParse.H> // used in initialize function

AmrBoxLib::AmrBoxLib() : AmrObject(),
                         AmrCore(),
                         nChargePerCell_m(PArrayManage),
                         bunch_mp(nullptr),
                         rho_m(PArrayManage),
                         phi_m(PArrayManage),
                         eg_m(PArrayManage)
{}


AmrBoxLib::AmrBoxLib(TaggingCriteria tagging,
                     double scaling,
                     double nCharge)
    : AmrObject(tagging, scaling, nCharge),
      AmrCore(),
      nChargePerCell_m(PArrayManage),
      bunch_mp(nullptr),
      rho_m(PArrayManage),
      phi_m(PArrayManage),
      eg_m(PArrayManage)
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
      nChargePerCell_m(PArrayManage),
      bunch_mp(nullptr),
      rho_m(PArrayManage),
      phi_m(PArrayManage),
      eg_m(PArrayManage)
{}


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
    for (int lev = 0; lev < eg_m.size(); ++lev) {
        for (std::size_t i = 0; i < bunch_mp->Dimension; ++i) {
            /* calls:
             *  - max(comp, nghost (to search), local)
             *  - min(cmop, nghost, local)
             */
            double max = eg_m[lev].max(i, 0, false);
            maxE[i] = (maxE[i] < max) ? max : maxE[i];
            
            double min = eg_m[lev].min(i, 0, false);
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
//     IpplTimings::startTimer(selfFieldTimer_m);
    
    //scatter charges onto grid
    bunch_mp->Q *= bunch_mp->dt;
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->scatter(bunch_mp->Q, this->rho_m, bunch_mp->R, 0, -1);
    bunch_mp->Q /= bunch_mp->dt;
    
    int baseLevel = 0;
    int finestLevel = (&amrpbase_p->getAmrLayout())->finestLevel();
    
    int nLevel = finestLevel + 1;
    rho_m.resize(nLevel);
    phi_m.resize(nLevel);
    eg_m.resize(nLevel);
    
    double invDt = 1.0 / bunch_mp->getdT() * bunch_mp->getCouplingConstant();
    for (int i = 0; i <= finestLevel; ++i) {
        this->rho_m[i].mult(invDt, 0, 1);
    }
    
    // charge density is in rho_m
//     IpplTimings::startTimer(compPotenTimer_m);
    
    //FIXME
//     fs_m->solver_m->solve(rho_m, eg_m, baseLevel, finestLevel);

//     IpplTimings::stopTimer(compPotenTimer_m);
    
//     IpplTimings::stopTimer(selfFieldTimer_m);
}


void AmrBoxLib::computeSelfFields(int b) {
    //TODO
}


void AmrBoxLib::computeSelfFields_cycl(double gamma) {
    //TODO
}


void AmrBoxLib::computeSelfFields_cycl(int b) {
    //TODO
}


void AmrBoxLib::RemakeLevel (int lev, Real time,
                             const BoxArray& new_grids,
                             const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m.clear(lev);
    nChargePerCell_m.set(lev, new AmrField_t(new_grids, 1, 1, new_dmap));
}


void AmrBoxLib::MakeNewLevel (int lev, Real time,
                              const BoxArray& new_grids,
                              const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m.set(lev, new AmrField_t(new_grids, 1, 1, dmap[lev]));
}


void AmrBoxLib::ClearLevel(int lev) {
    nChargePerCell_m.clear(lev);
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
        nChargePerCell_m[i].setVal(0.0);
        amrpbase_p->AssignDensitySingleLevel(bunch_mp->Q, nChargePerCell_m[i], i);
    }
    
    for (int i = finest_level-1; i >= lev; --i) {
        AmrField_t tmp(nChargePerCell_m[i].boxArray(), 1, 0, nChargePerCell_m[i].DistributionMap());
        tmp.setVal(0.0);
        BoxLib::average_down(nChargePerCell_m[i+1], tmp, 0, 1, refRatio(i));
        AmrField_t::Add(nChargePerCell_m[i], tmp, 0, 0, 1, 0);
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
        for (MFIter mfi(nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((nChargePerCell_m[lev])[mfi]),
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
    int base_level   = 0;
    int finest_level = 0;
    
    AmrFieldContainer_t nPartPerCell(1, PArrayManage);
    nPartPerCell.set(0, new AmrField_t(this->boxArray(lev), 1, 1,
                                       this->DistributionMap(lev)));
    nPartPerCell[base_level].setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->AssignDensitySingleLevel(bunch_mp->Q, nPartPerCell[base_level], lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int lev = base_level; lev <= finest_level; lev++)
            nPartPerCell[lev].plus(+offset, 0, 1, 0);
    }
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi(1, PArrayManage);
    AmrFieldContainer_t efield(1, PArrayManage);
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, base_level, finest_level);
    
    // 3. Tag cells for refinement
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    
    Real minp = 0.0;
    Real maxp = 0.0;
    maxp = scaling_m * phi[base_level].max(0);
    minp = scaling_m * phi[base_level].min(0);
    Real threshold = std::max( std::abs(minp), std::abs(maxp) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        for (MFIter mfi(phi[base_level],false); mfi.isValid(); ++mfi) {
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
                        BL_TO_FORTRAN_3D((phi[base_level])[mfi]),
                        &tagval, &clearval, 
                        ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                        ZFILL(dx), ZFILL(prob_lo), &time, &threshold);
            //
            // Now update the tags in the TagBox.
            //
            tagfab.tags_and_untags(itags, tilebx);
        }
    }
    
    phi.clear(base_level);
    efield.clear(base_level);
    nPartPerCell.clear(0);
}


void AmrBoxLib::tagForEfield_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    /* Perform a single level solve a level lev and tag all cells for refinement
     * where the value of the efield is higher than 75 percent of the maximum efield
     * value of this level.
     */
    
    // 1. Assign charge onto grid of level lev
    int base_level   = 0;
    int finest_level = 0;
    
    AmrFieldContainer_t nPartPerCell(PArrayManage);
    nPartPerCell.resize(1);
    nPartPerCell.set(0, new AmrField_t(this->boxArray(lev), 1, 1,
                                       this->DistributionMap(lev)));
    nPartPerCell[base_level].setVal(0.0);
    
    AmrPartBunch::pbase_t* amrpbase_p = bunch_mp->getAmrParticleBase();
    amrpbase_p->AssignDensitySingleLevel(bunch_mp->Q, nPartPerCell[base_level], lev);
    
    if ( geom[0].isAllPeriodic() ) {
        Real offset = 0.0;
        for (std::size_t i = 0; i < bunch_mp->getLocalNum(); ++i)
            offset += bunch_mp->Q[i];
        offset /= geom[0].ProbSize();
        for (int lev = base_level; lev <= finest_level; lev++)
            nPartPerCell[lev].plus(+offset, 0, 1, 0);
    }
    
    // 2. Solve Poisson's equation on level lev
    AmrFieldContainer_t phi(1, PArrayManage);
    AmrFieldContainer_t efield(1, PArrayManage);
    
//     AmrPoissonSolver<AmrBoxLib> *solver = bunch_mp->template getFieldSolver<AmrBoxLib>();
    PoissonSolver *solver = bunch_mp->getFieldSolver();
    
    solver->solve(nPartPerCell, phi, efield, base_level, finest_level);
    
    // 3. Tag cells for refinement
    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;
    
    Real min[3] = {0.0, 0.0, 0.0};
    Real max[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        max[i] = scaling_m * efield[base_level].max(i);
        min[i] = scaling_m * efield[base_level].min(i);
    }
    Real threshold[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i)
        threshold[i] = std::max( std::abs(min[i]), std::abs(max[i]) );
    
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;
        // mfi(efield[base_level], true)
        for (MFIter mfi(efield[base_level],false); mfi.isValid(); ++mfi) {
            const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
            
            TagBox&     tagfab  = tags[mfi];
            FArrayBox&  fab     = (efield[base_level])[mfi];
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
    
    phi.clear(base_level);
    efield.clear(base_level);
    nPartPerCell.clear(base_level);
    
}
