#ifndef AMR_OBECT_H
#define AMR_OBECT_H

#include "Index/NDIndex.h"

#include "Algorithms/PBunchDefs.h"
#include "Utilities/Util.h"

// #include "Algorithms/AmrPartBunch.h"

// class AmrPartBunch;

/**
 * The AMR interface to OPAL. A new AMR library needs
 * to inherit from this class in order to work properly
 * with OPAL. Among other things it specifies the refinement
 * strategies.
 */
class AmrObject {
    
public:
    // FIXME Why not using typedef of PartBunchBase
    typedef std::pair<Vector_t, Vector_t> VectorPair_t;
//     using VectorPair_t = typename AmrPartBunch::VectorPair_t;
    
    /// Methods for tagging cells for refinement
    enum TaggingCriteria {
        CHARGE_DENSITY = 0, // default
        POTENTIAL,
        EFIELD,
        MOMENTA,
        MAX_NUM_PARTICLES,      ///< max. #particles per cell
        MIN_NUM_PARTICLES       ///< min. #particles per cell
    };
    
    
    /*!
     * This data structure is only used for creating an object
     * via the static member function AmrBoxLib::create()
     * that is called in FieldSolver::initAmrObject_m
     */
    struct AmrInfo {
        int grid[3];        ///< Number of grid points in x-, y- and z-direction
        int maxgrid[3];     ///< Maximum grid size in x-, y- and z-direction
        int bf[3];          ///< Grid blocking factor in x-, y- and z-direction
        int maxlevel;       ///< Maximum level for AMR (0: single-level)
        int refratio[3];    ///< Mesh refinement ratio in x-, y- and z-direction
    };
    
public:
    
    AmrObject() : tagging_m(CHARGE_DENSITY),
                  scaling_m(0.75),
                  chargedensity_m(1.0e-15),
                  maxNumPart_m(1),
                  minNumPart_m(1),
                  refined_m(false)
    {}
    
    AmrObject(TaggingCriteria tagging,
              double scaling,
              double chargedensity_m) : tagging_m(tagging),
                                        scaling_m(scaling),
                                        chargedensity_m(chargedensity_m),
                                        maxNumPart_m(1),
                                        minNumPart_m(1),
                                        refined_m(false)
    {}
    
    virtual ~AmrObject() {}
    
    /*!
     * Setup all fine levels after object creation.
     */
    virtual void initFineLevels() = 0; 
    
    /*!
     * Update of mesh according to chosen refinement strategy.
     * @param lbase base level to regrid
     * @param lfine finest level to regrid
     * @param time of regrid
     */
    virtual void regrid(int lbase, int lfine, double time) = 0;
    
    /*!
     * Choose a new tagging strategy.
     * Is used in src/Structure/FieldSolver.cpp
     * @param tagging strategy
     */
    void setTagging(TaggingCriteria tagging) {
        tagging_m = tagging;
    }
    
    /*!
     * Choose a new tagging strategy (string version).
     * Is used in src/Structure/FieldSolver.cpp
     * @param tagging strategy
     */
    void setTagging(std::string tagging) {
        tagging = Util::toUpper(tagging);
        
        if ( tagging == "POTENTIAL" )
            tagging_m = TaggingCriteria::POTENTIAL;
        else if (tagging == "EFIELD" )
            tagging_m = TaggingCriteria::EFIELD;
        else if ( tagging == "MOMENTA" )
            tagging_m = TaggingCriteria::MOMENTA;
        else if ( tagging == "MAX_NUM_PARTICLES" )
            tagging_m = TaggingCriteria::MAX_NUM_PARTICLES;
        else if ( tagging == "MIN_NUM_PARTICLES" )
            tagging_m = TaggingCriteria::MIN_NUM_PARTICLES;
        else if ( tagging == "CHARGE_DENSITY" )
            tagging_m = TaggingCriteria::CHARGE_DENSITY;
        else
            throw OpalException("AmrObject::setTagging(std::string)",
                                "Not supported refinement criteria "
                                "[CHARGE_DENSITY | POTENTIAL | EFIELD | "
                                "MOMENTA | MAX_NUM_PARTICLES | MIN_NUM_PARTICLES].");
    }
    
    /*!
     * Scaling factor for tagging.
     * It is used with POTENTIAL and EFIELD
     * @param scaling factor in [0, 1]
     */
    void setScalingFactor(double scaling) {
        scaling_m = scaling;
    }
    
    /*!
     * Charge density for tagging with CHARGE_DENSITY
     * @param chargedensity >= 0.0 (e.g. 1e-14)
     */
    void setChargeDensity(double chargedensity) {
        chargedensity_m = chargedensity;
    }
    
    /*!
     * Maximum number of particles per cell for tagging
     * @param maxNumPart is upper bound for a cell to be marked
     * for refinement
     */
    void setMaxNumParticles(size_t maxNumPart) {
        maxNumPart_m = maxNumPart;
    }
    
    /*!
     * Minimum number of particles per cell for tagging
     * @param minNumPart is lower bound for a cell to be marked
     * for refinement
     */
    void setMinNumParticles(size_t minNumPart) {
        minNumPart_m = minNumPart;
    }
        
    /* Methods that are needed by the
     * bunch
     */
    virtual VectorPair_t getEExtrema() = 0;
    
    virtual double getRho(int x, int y, int z) = 0;
    
    virtual void computeSelfFields() = 0;
    
    virtual void computeSelfFields(int b) = 0;
    
    virtual void computeSelfFields_cycl(double gamma) = 0;
    
    virtual void computeSelfFields_cycl(int b) = 0;
    
    virtual void updateMesh() = 0;
    
    virtual Vektor<int, 3> getBaseLevelGridPoints() const = 0;
    
    virtual const int& maxLevel() const = 0;
    virtual const int& finestLevel() const = 0;
    
    /*!
     * Rebalance the grids among the
     * cores
     */
    virtual void redistributeGrids(int how) { }
    
    /*!
     * Used in AmrPartBunch to check if we need to refine
     * first.
     * @returns true fine grids are initialized
     */
    const bool& isRefined() const {
        return refined_m;
    }
    
protected:
    TaggingCriteria tagging_m;  ///< Tagging strategy
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (POTENTIAL, EFIELD)
    double chargedensity_m;     ///< Tagging value for CHARGE_DENSITY
    
    size_t maxNumPart_m;        ///< Tagging value for MAX_NUM_PARTICLES
    
    size_t minNumPart_m;        ///< Tagging value for MIN_NUM_PARTICLES
    
    bool refined_m;             ///< Only set to true in AmrObject::initFineLevels()
};

#endif
