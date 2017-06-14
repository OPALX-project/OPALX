#ifndef AMR_BOXLIB_H
#define AMR_BOXLIB_H

#include "Amr/AmrObject.h"

class AmrPartBunch;

// AMReX headers
#include <AMReX_AmrMesh.H>
#include <AMReX.H>

/*!
 * Concrete AMR object. It is based on the
 * <a href="https://ccse.lbl.gov/AMReX/">AMReX</a> library
 * developed at LBNL. This library is the successor of
 * <a href="https://ccse.lbl.gov/BoxLib/">BoxLib</a>.
 */
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
    
    /*!
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
    
public:
    
    /*!
     * Using this constructors leads to an invalid AmrBoxLib object.
     */
    AmrBoxLib();
    
    /*!
     * @param tagging different mesh refinement strategies
     * @param scaling is used in potential and electric field tagging where
     * all cells are marked for refinement if they have a value greater equal than
     * either
     * \f[
     *      \alpha\cdot\max_{i,j,k} |\phi_{i,j,k}|
     * \f]
     * in case of potential tagging or
     * \f[
     *      \alpha\cdot\max_{i,j,k} |\vec{E}_{i,j,k}|
     * \f]
     * in case of electric field tagging.\ The scalar \f$\alpha\f$ represents the
     * scaling value.\ In case of electric field tagging, each component is treated
     * independently.
     * @param nCharge is the amount of charge that a cell has to have in order to be
     * refined.\ The cell is marked for refinement if the value is greater equatl
     * to nCharge [C / m].
     */
    AmrBoxLib(TaggingCriteria tagging,
              double scaling,
              double nCharge);
    
    /*!
     * @param domain is the physical domain of the problem. In case of
     * AMReX the domain is specified by src/Amr/BoxLibLayout.h. The particles
     * are mapped to the domain \f$[-1, 1]^3\f$, thus the domain is a tiny bit
     * greater than that.
     * @param nGridPts per dimension (nx, ny, nz / nt)
     * @param maxLevel of mesh refinement
     */
    AmrBoxLib(const AmrDomain_t& domain,
              const AmrIntArray_t& nGridPts,
              short maxLevel);
    
    /*!
     * See other constructors documentation for further info.
     * @param domain is the physical domain of the problem
     * @param nGridPts per dimension (nx, ny, nz / nt)
     * @param maxLevel of mesh refinement
     * @param bunch is used when we tag for charges per cell
     */
    AmrBoxLib(const AmrDomain_t& domain,
              const AmrIntArray_t& nGridPts,
              short maxLevel,
              AmrPartBunch* bunch_p);
    
    /*!
     * Create a new object
     * @param info are the initial informations to construct an AmrBoxLib object.
     * It is set in src/Solvers/FieldSolver.cpp
     * @param bunch_p pointing to. Is used for getting the geometry (i.e. physical
     * domain of the problem).
     */
    static std::unique_ptr<AmrBoxLib> create(const AmrInitialInfo& info,
                                             AmrPartBunch* bunch_p);
    
    /*!
     * Inherited from AmrObject
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
    
protected:
    /*
     * AmrMesh functions
     */
    
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
     * Clean up a level.
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
     * 
     * Remark: Not used in OPAL
     * 
     * @param lev to create
     * @param time of simulation
     * @param ba the boxes
     * @param dm the grid distribution among cores
     */
    void MakeNewLevelFromScratch (int lev, AmrReal_t time,
                                  const AmrGrid_t& ba,
                                  const AmrProcMap_t& dm);

    /*!
     * Make a new level using provided BoxArray and
     * DistributionMapping and fill with interpolated coarse level data.
     * 
     * Remark: Not used in OPAL
     * 
     * @param lev to create
     * @param time of simulation
     * @param ba the boxes
     * @param dm the grid distribution among cores
     */
    void MakeNewLevelFromCoarse (int lev, AmrReal_t time,
                                 const AmrGrid_t& ba,
                                 const AmrProcMap_t& dm);
    
private:
    /*!
     * Mark a cell for refinement if the value is greater equal
     * than some amount of charge (AmrObject::nCharge_m).
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForChargeDensity_m(int lev, TagBoxArray_t& tags,
                               AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if the potential value is greater
     * equal than the maximum value of the potential on the grid
     * scaled by some factor [0, 1] (AmrObject::scaling_m)
     * It solves the Poisson equation on that level.
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForPotentialStrength_m(int lev, TagBoxArray_t& tags,
                                   AmrReal_t time, int ngrow);
    
    /*!
     * Mark a cell for refinement if one of the electric field components
     * is greater equal the maximum electric field value per direction scaled
     * by some factor [0, 1] (AmrObject::scaling_m). It solves the Poisson
     * equation on that level.
     * 
     * @param lev to check for refinement
     * @param tags is a special box array that marks cells for refinement
     * @param time of simulation (not used)
     * @param ngrow is the number of ghost cells (not used)
     */
    void tagForEfield_m(int lev, TagBoxArray_t& tags,
                        AmrReal_t time, int ngrow);
    
    /*!
     * Use particle BoxArray and DistributionMapping for AmrObject and
     * reset geometry for bunch
     * 
     * @param nGridPts per dimension (nx, ny, nz / nt)
     */
    void initBaseLevel_m(const AmrIntArray_t& nGridPts);
    
    
private:
    /// use in tagging tagForChargeDensity_m (needed when tracking)
    AmrFieldContainer_t nChargePerCell_m;
    
    /// bunch used for tagging strategies
    AmrPartBunch *bunch_mp;
    
    // the layout of the bunch
    AmrLayout_t  *layout_mp;
    
    /// charge density on the grid for all levels
    AmrFieldContainer_t rho_m;
    
    /// scalar potential on the grid for all levels
    AmrFieldContainer_t phi_m;
    
    /// vector field on the grid for all levels
    AmrFieldContainer_t eg_m;
    
    /*!
     * used for writing charge density, potential, electric field to a
     * file (compile with CMAKE option -DDBG_SCALARFIELD=1)
     */
    int fieldDBGStep_m;
};

#endif
