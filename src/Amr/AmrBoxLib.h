#ifndef AMR_BOXLIB_H
#define AMR_BOXLIB_H

#include "Amr/AmrObject.h"

// BoxLib headers
#include <AmrCore.H>
#include <BoxLib.H>
#include <Geometry.H>

class AmrBoxLib : public AmrObject,
                  public AmrCore
{
    
public:
    typedef MultiFab AmrField_t;
    typedef PArray<AmrField_t> AmrFieldContainer_t;
    
public:
    
    /*!
     * Set all parameters for the AMR object, like #grid points per dimension etc.
     * (inherited from AmrObject)
     */
    void initialize();
    
    /*!
     * @param lbase start of regridding.
     * @param lfine end of regridding.
     * @param time of simulation (step).
     */
    void regrid(int lbase, int lfine, double time);
    
protected:
    /*!
     * Update the grids and the distributionmapping for a specific level (inherited from AmrCore)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void RemakeLevel (int lev, Real time,
                      const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    /*!
     * Create completeley new grids for a level (inherited from AmrCore)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    /*!
     * @param lev to free allocated memory.
     */
    void ClearLevel(int lev);
    
    /*!
     * Is called in the AmrCore function for performing tagging. (inherited from AmrCore)
     */
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) override;
    
private:
    void tagForChargeDensity_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    void tagForPotentialStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    void tagForEfieldGradient_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    
private:
    /// use in tagging tagForChargeDensity_m (needed when tracking)
    AmrFieldContainer_t nChargePerCell_m;
};

#endif
