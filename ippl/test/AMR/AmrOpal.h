#ifndef AMROPAL_H
#define AMROPAL_H

#include <AMReX_AmrCore.H>

#ifdef IPPL_AMR
    #include "ippl-amr/AmrParticleBase.h"
    #include "ippl-amr/ParticleAmrLayout.h"
    #include "ippl-amr/PartBunchAmr.h"
#else
    #include "boxlib-amr/AmrPartBunch.h"
#endif

#include <AMReX_MultiFabUtil.H>

#include <memory>

using namespace amrex;

/*!
 * @file AmrOpal.h
 * @authors Matthias Frey
 *          Ann Almgren
 *          Weiqun Zhang
 * @date October 2016, LBNL
 * @details This class implements the abstract base class
 * AmrCore of the BoxLib library. It defines functionalities
 * to do Amr calculations.
 * @brief Concrete AMR implementation
 */

/// Concrete AMR implementation
class AmrOpal : public AmrMesh {
    
private:
    typedef Array<std::unique_ptr<MultiFab> > mfs_mt;

public:
    /// Methods for tag cells for refinement
    enum TaggingCriteria {
        kChargeDensity = 0, // default
        kPotentialStrength,
        kEfieldStrength,
        kCenteredRegion     // only for boxlib-only/testDeposition comparison
    };
        
    
#ifdef IPPL_AMR
    typedef ParticleAmrLayout<double, BL_SPACEDIM> amrplayout_t;
    typedef AmrParticleBase<amrplayout_t> amrbase_t;
    typedef PartBunchAmr<amrplayout_t> amrbunch_t;
#endif
    
    
//     typedef PArray<MultiFab> mp_mt;
    
public:
    /*!
     * Create an AMR object.
     * @param rb is the physical domain
     * @param max_level_in is the max. number of allowed AMR levels
     * @param n_cell_in is the number of grid cells at the coarsest level
     * @param coord is the coordinate system (0: cartesian)
     * @param bunch is the particle bunch
     */
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord,
#ifdef IPPL_AMR
            PartBunchAmr<amrplayout_t>* bunch);
#else
            PartBunchBase* bunch);
#endif

    /*!
     * Create an AMR object.
     * @param rb is the physical domain
     * @param max_level_in is the max. number of allowed AMR levels
     * @param n_cell_in is the number of grid cells at the coarsest level
     * @param coord is the coordinate system (0: cartesian)
     * @param refratio
     */
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord,
            const std::vector<int>& refratio);
    
    /*!
     * Create an AMR object.
     * @param rb is the physical domain
     * @param max_level_in is the max. number of allowed AMR levels
     * @param n_cell_in is the number of grid cells at the coarsest level
     * @param coord is the coordinate system (0: cartesian)
     */
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);
    
    virtual ~AmrOpal();     ///< does nothing
    
    void initBaseLevel();   ///< defines BoxArray, DistributionMapping of level 0
    
    /*!
     * @param lbase is the level on top of which a new level is created
     * @param time not used
     */
    void regrid (int lbase, Real time);
    
    /*!
     * Update the grids and the distributionmapping for a specific level
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void RemakeLevel (int lev, Real time,
                      const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    /*!
     * Create completeley new grids for a level
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    
    void ClearLevel(int lev);
    
    void setBunch(
#ifdef IPPL_AMR
        PartBunchAmr<amrplayout_t>* bunch)
#else
        AmrPartBunch* bunch)
#endif
    {
        bunch_m = bunch;
    }
    
    /*!
     * Write a timestamp file for displaying with yt.
     */
    void writePlotFileYt(std::string filename, int step);
    
    /*!
     * Write a timestamp file for displaying with AmrVis.
     */
    void writePlotFile(std::string filename, int step);
    
    void setTagging(TaggingCriteria tagging) {
        tagging_m = tagging;
    }
    
    /*!
     * Scaling factor for tagging.
     * It is used in tagForPotentialStrength_m and tagForEfieldStrength_m
     * @param scaling factor in [0, 1]
     */
    void setScalingFactor(double scaling) {
        scaling_m = scaling;
    }
    
    /*!
     * Charge for tagging in tagForChargeDensity_m
     * @param charge >= 0.0 (e.g. 1e-14)
     */
    void setCharge(double charge) {
        nCharge_m = charge;
    }
    
protected:
    /*!
     * Is called in the AmrCore function for performing tagging.
     */
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) override;
    
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm);

    //! Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data.
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm);
    
private:
    
    void tagForChargeDensity_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    void tagForPotentialStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    void tagForEfieldStrength_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    void tagForCenteredRegion_m(int lev, TagBoxArray& tags, Real time, int ngrow);
    
    
    
#ifdef IPPL_AMR
    PartBunchAmr<amrplayout_t>* bunch_m;
#else
    AmrPartBunch* bunch_m;      ///< Particle bunch
#endif
    TaggingCriteria tagging_m;
    mfs_mt nChargePerCell_m;    ///< use in tagging tagForChargeDensity_m (needed when tracking)
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (tagForPotentialStrength_m, tagForEfieldStrength_m)
    Real   nCharge_m;           ///< Tagging value for tagForChargeDensity_m
    
};

#endif