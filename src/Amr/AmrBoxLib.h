#ifndef AMR_BOXLIB_H
#define AMR_BOXLIB_H

#include "Amr/AmrObject.h"

class AmrPartBunch;

// BoxLib headers
#include <AMReX_AmrMesh.H>
#include <AMReX.H>

class AmrBoxLib : public AmrObject,
                  public amrex::AmrMesh
{
    
public:
    typedef amr::AmrField_t             AmrField_t;
    typedef amr::AmrFieldContainer_t    AmrFieldContainer_t;
    typedef amr::AmrGeomContainer_t     AmrGeomContainer_t;
    typedef amr::AmrGridContainer_t     AmrGridContainer_t;
    typedef amr::AmrProcMapContainer_t  AmrProcMapContainer_t;
    typedef amr::AmrDomain_t            AmrDomain_t;
    typedef amr::AmrIntArray_t          AmrIntArray_t;
    typedef amr::AmrReal_t              AmrReal_t;
    typedef amr::AmrGrid_t              AmrGrid_t;
    typedef amr::AmrProcMap_t           AmrProcMap_t;
    typedef amr::AmrGeometry_t          AmrGeometry_t;
    typedef amr::AmrIntVect_t           AmrIntVect_t;
    
    typedef amrex::TagBoxArray          TagBoxArray_t;
    
    /**
     * This data structure is only used for creating an object
     * via the static member fucction AmrBoxLib::create()
     * that is called in FieldSolver::initAmrObject_m
     */
    struct AmrInitialInfo {
        int gridx;          ///< Number of grid points in x-direction
        int gridy;          ///< Number of grid points in y-direction
        int gridz;          ///< Number of grid points in z-direction
        int maxgrid;        ///< Maximum grid size allowed
        int maxlevel;       ///< Maximum level for AMR
        int refratx;        ///< Mesh refinement ratio in x-direction
        int refraty;        ///< Mesh refinement ratio in y-direction
        int refratz;        ///< Mesh refinement ratio in z-direction
    };
    
//     typedef typename AmrPartBunch::VectorPair_t VectorPair_t;
    
public:
    
//     AmrBoxLib(const DomainBoundary_t& realbox,
//               const NDIndex<3>& nGridPts,
//               short maxLevel,
//               const RefineRatios_t& refRatio);
    
    AmrBoxLib();
    
    AmrBoxLib(TaggingCriteria tagging,
              double scaling,
              double nCharge);
    
    AmrBoxLib(const AmrDomain_t& domain,
              const AmrIntArray_t& nGridPts,
              short maxLevel);
    
    AmrBoxLib(const AmrDomain_t& domain,
              const AmrIntArray_t& nGridPts,
              short maxLevel,
              AmrPartBunch* bunch);
    
    /**
     * Create a new object
     * @param fs having all data to create the object
     * @param bunch_p pointing to
     */
    static std::unique_ptr<AmrBoxLib> create(const AmrInitialInfo& fs, AmrPartBunch* bunch_p);
    
    inline void setBunch(AmrPartBunch* bunch);
    
//     /*!
//      * Set all parameters for the AMR object, like #grid points per dimension etc.
//      * (inherited from AmrObject)
//      */
//     void initialize(const DomainBoundary_t& lower,
//                     const DomainBoundary_t& upper,
//                     const NDIndex<3>& nGridPts,
//                     short maxLevel,
//                     const RefineRatios_t& refRatio);
    
    /*!
     * @param lbase start of regridding.
     * @param lfine end of regridding.
     * @param time of simulation (step).
     */
    void regrid(int lbase, int lfine, double time);
    
    
    VectorPair_t getEExtrema();
    
    double getRho(int x, int y, int z);
    
    void computeSelfFields();
    
    void computeSelfFields(int bin);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int bin);
    
    void updateMesh();
    
    Vektor<int, 3> getBaseLevelGridPoints();
    
    int maxLevel();
    int finestLevel();
    
//     void updateBunch();
    
protected:
    /*!
     * Update the grids and the distributionmapping for a specific level
     * (inherited from AmrMesh)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void RemakeLevel (int lev, AmrReal_t time,
                      const AmrGrid_t& new_grids, const AmrProcMap_t& new_dmap);
    
    /*!
     * Create completeley new grids for a level (inherited from AmrMesh)
     * @param lev is the current level
     * @param time not used
     * @param new_grids are the new created grids for this level
     * @param new_dmap is the new distribution to processors
     */
    void MakeNewLevel (int lev, AmrReal_t time,
                       const AmrGrid_t& new_grids, const AmrProcMap_t& new_dmap);
    
    /*!
     * @param lev to free allocated memory.
     */
    void ClearLevel(int lev);
    
    /*!
     * Is called in the AmrMesh function for performing tagging. (inherited from AmrMesh)
     */
    virtual void ErrorEst(int lev, TagBoxArray_t& tags,
                          AmrReal_t time, int ngrow) override;
    
    
    /*!
     * Make a new level from scratch using provided BoxArray and
     * DistributionMapping.
     * Only used during initialization.
     */
    void MakeNewLevelFromScratch (int lev, AmrReal_t time,
                                  const AmrGrid_t& ba,
                                  const AmrProcMap_t& dm);

    /*!
     * Make a new level using provided BoxArray and
     * DistributionMapping and fill with interpolated coarse level data.
     */
    void MakeNewLevelFromCoarse (int lev, AmrReal_t time,
                                 const AmrGrid_t& ba,
                                 const AmrProcMap_t& dm);
    
private:
    void tagForChargeDensity_m(int lev, TagBoxArray_t& tags,
                               AmrReal_t time, int ngrow);
    
    void tagForPotentialStrength_m(int lev, TagBoxArray_t& tags,
                                   AmrReal_t time, int ngrow);
    
    void tagForEfield_m(int lev, TagBoxArray_t& tags,
                        AmrReal_t time, int ngrow);
    
    /*!
     * Use particle BoxArray and DistributionMapping for AmrObject and
     * reset geometry for bunch
     */
    void initBaseLevel_m(const AmrIntArray_t& nGridPts);
    
//     void resizeBaseLevel_m(const AmrDomain_t& domain,
//                            const AmrIntArray_t& nGridPts);
    
private:
    /// use in tagging tagForChargeDensity_m (needed when tracking)
    AmrFieldContainer_t nChargePerCell_m;
    
    AmrPartBunch *bunch_mp;
    AmrLayout_t  *layout_mp;
    
    /// charge density on the grid for all levels
    AmrFieldContainer_t rho_m;
    
    /// scalar potential on the grid for all levels
    AmrFieldContainer_t phi_m;
    
    /// vector field on the grid for all levels
    AmrFieldContainer_t eg_m;
    
    int fieldDBGStep_m;
};

#endif
