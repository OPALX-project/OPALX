#include "AmrBoxLib.h"

#include <ParmParse.H> // used in initialize function


void AmrBoxLib::initialize() {
    //TODO
}


void AmrBoxLib::regrid(int lbase, int lfine, double time) {
    int new_finest = 0;
    Array<BoxArray> new_grids(finest_level+2);
    
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
//                 amrplayout_t* PLayout = &bunch_m->getLayout();
//                 PLayout->SetParticleBoxArray(lev, new_grids[lev]);
//                 PLayout->SetParticleDistributionMap(lev, new_dmap);
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
//             amrplayout_t* PLayout = &bunch_m->getLayout();
//             PLayout->SetParticleBoxArray(lev, new_grids[lev]);
//             PLayout->SetParticleDistributionMap(lev, new_dmap);
        }
    }
    
    finest_level = new_finest;
    
//     // update to multilevel
//     bunch_m->update();
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
        case kChargeDensity:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
        case kPotentialStrength:
            tagForPotentialStrength_m(lev, tags, time, ngrow);
            break;
        case kEfieldGradient:
            tagForEfieldGradient_m(lev, tags, time, ngrow);
            break;
        default:
            tagForChargeDensity_m(lev, tags, time, ngrow);
            break;
    }
}


void AmrBoxLib::tagForChargeDensity_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    //TODO
}


void AmrBoxLib::tagForPotentialStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    //TODO
}


void AmrBoxLib::tagForEfieldGradient_m(int lev, TagBoxArray& tags, Real time, int ngrow) {
    //TODO
}
