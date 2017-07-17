#ifndef AMR_OBECT_H
#define AMR_OBECT_H

#include "Index/NDIndex.h"

#include "Algorithms/PBunchDefs.h"

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
        EFIELD
    };
    
public:
    
    AmrObject() : tagging_m(CHARGE_DENSITY),
                  scaling_m(0.75),
                  nCharge_m(1.0e-15),
                  regridFreq_m(1)
    {}
    
    AmrObject(TaggingCriteria tagging,
              double scaling,
              double nCharge) : tagging_m(tagging),
                                scaling_m(scaling),
                                nCharge_m(nCharge),
                                regridFreq_m(1)
    {}
    
    virtual ~AmrObject() {}
    
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
     * Scaling factor for tagging.
     * It is used with POTENTIAL and EFIELD
     * @param scaling factor in [0, 1]
     */
    void setScalingFactor(double scaling) {
        scaling_m = scaling;
    }
    
    /*!
     * Charge for tagging with CHARGE_DENSITY
     * @param charge >= 0.0 (e.g. 1e-14)
     */
    void setCharge(double charge) {
        nCharge_m = charge;
    }
    
    /*!
     * After how many steps a regrid
     * should be performed
     * @param regridFreq is set in FieldSolver.cpp
     */
    void setRegridFrequency(size_t regridFreq) {
        regridFreq_m = regridFreq;
    }
    
    /*!
     * @returns the frequency, used in AmrPartBunch::boundp()
     */
    const size_t& getRegridFrequency() {
        return regridFreq_m;
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
    
    virtual Vektor<int, 3> getBaseLevelGridPoints() = 0;
    
    virtual int maxLevel() = 0;
    virtual int finestLevel() = 0;
    
    /*!
     * Rebalance the grids among the
     * cores
     */
    virtual void redistributeGrids(int how) { }
    
protected:
    TaggingCriteria tagging_m;  ///< Tagging strategy
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (POTENTIAL, EFIELD)
    double nCharge_m;           ///< Tagging value for CHARGE_DENSITY
    
    size_t regridFreq_m;        ///< After how many time steps to regrid
    
};

#endif
