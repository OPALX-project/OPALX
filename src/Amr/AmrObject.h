#ifndef AMR_OBECT_H
#define AMR_OBECT_H

#include "Index/NDIndex.h"

#include "Algorithms/PBunchDefs.h"

// #include "Algorithms/AmrPartBunch.h"

// class AmrPartBunch;

#include <array>

class AmrObject {
    
public:
    // FIXME Why not using typedef of PartBunchBase
    typedef std::pair<Vector_t, Vector_t> VectorPair_t;
//     using VectorPair_t = typename AmrPartBunch::VectorPair_t;
    
//     typedef std::array<short, 3> RefineRatios_t;
//     typedef std::array<double, 3> DomainBoundary_t;
    
    /// Methods for tagging cells for refinement
    enum TaggingCriteria {
        CHARGE_DENSITY = 0, // default
        POTENTIAL,
        EFIELD
    };
    
public:
    
    AmrObject() : tagging_m(CHARGE_DENSITY),
                  scaling_m(0.75),
                  nCharge_m(1.0e-15)
    {}
    
    AmrObject(TaggingCriteria tagging,
              double scaling,
              double nCharge) : tagging_m(tagging),
                                scaling_m(scaling),
                                nCharge_m(nCharge)
    {}
    
//     /*!
//      * @param realbox is the physical box extent
//      * @param nGridPts in each dimension of base level.
//      * @param maxLevel of AMR
//      * @param refRatio in x, y, z between levels
//      */
//     AmrObject(const DomainBoundary_t& realbox,
//               const NDIndex<3>& nGridPts,
//               short maxLevel,
//               const RefineRatios_t& refRatio)
//         : tagging_m(CHARGE_DENSITY),
//           scaling_m(0.75),
//           nCharge_m(1.0e-15)
//     {}
    
    virtual ~AmrObject() {}
    
//     /*!
//      * This function should be called in case
//      * the AmrObject is constructed using either
//      * one of the following constructors:
//      * - AmrObject()
//      * - AmrObject(TaggingCriteria, scaling, nCharge)
//      * @param lower physical domain boundary
//      * @param upper physical domain boundary
//      * @param nGridPts in each dimension of base level
//      * @param maxLevel of AMR
//      * @param refRatio in x, y, z between levels
//      */
//     virtual void initialize(const DomainBoundary_t& lower,
//                             const DomainBoundary_t& upper,
//                             const NDIndex<3>& nGridPts,
//                             short maxLevel,
//                             const RefineRatios_t& refRatio)
//     {}
    
    virtual void regrid(int lbase, int lfine, double time) = 0;
    
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
    
    /* Methods that are needeb by the
     * bunch
     */
    virtual VectorPair_t getEExtrema() = 0;
    
    virtual double getRho(int x, int y, int z) = 0;
    
    virtual void computeSelfFields() = 0;
    
    virtual void computeSelfFields(int b) = 0;
    
    virtual void computeSelfFields_cycl(double gamma) = 0;
    
    virtual void computeSelfFields_cycl(int b) = 0;
    
protected:
    TaggingCriteria tagging_m;
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (POTENTIAL, EFIELD)
    double nCharge_m;           ///< Tagging value for CHARGE_DENSITY
    
};

#endif
