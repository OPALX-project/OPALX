#ifndef AMR_OBECT_H
#define AMR_OBECT_H

#include "Index/NDIndex.h"

#include <array>

class AmrObject {
    
public:
    
//     typedef std::array<short, 3> RefineRatios_t;
//     typedef std::array<double, 3> DomainBoundary_t;
    
    /// Methods for tagging cells for refinement
    enum TaggingCriteria {
        kChargeDensity = 0, // default
        kPotentialStrength,
        kEfieldGradient
    };
    
public:
    
    AmrObject() : tagging_m(kChargeDensity),
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
//         : tagging_m(kChargeDensity),
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
     * It is used with kPotentialStrength and kEfieldGradient
     * @param scaling factor in [0, 1]
     */
    void setScalingFactor(double scaling) {
        scaling_m = scaling;
    }
    
    /*!
     * Charge for tagging with kChargeDensity
     * @param charge >= 0.0 (e.g. 1e-14)
     */
    void setCharge(double charge) {
        nCharge_m = charge;
    }
    
protected:
    TaggingCriteria tagging_m;
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (kPotentialStrength, kEfieldGradient)
    double nCharge_m;           ///< Tagging value for kChargeDensity
    
};

#endif
