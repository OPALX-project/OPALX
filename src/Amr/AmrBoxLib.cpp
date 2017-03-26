#include "AmrBoxLib.h"

#include <ParmParse.H> // used in initialize function


void AmrBoxLib::initialize() {
    //TODO
}


void AmrBoxLib::regrid() {
    //TODO
}


void AmrBoxLib::RemakeLevel (int lev, Real time,
                             const BoxArray& new_grids,
                             const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m.clear(lev);
    nChargePerCell_m.set(lev, new MultiFab(new_grids, 1, 1, new_dmap));
}


void AmrBoxLib::MakeNewLevel (int lev, Real time,
                              const BoxArray& new_grids,
                              const DistributionMapping& new_dmap)
{
    SetBoxArray(lev, new_grids);
    SetDistributionMap(lev, new_dmap);
    
    nChargePerCell_m.set(lev, new MultiFab(new_grids, 1, 1, dmap[lev]));
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
