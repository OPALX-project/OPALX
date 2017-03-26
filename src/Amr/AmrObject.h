#ifndef AMR_OBECT_H
#define AMR_OBECT_H

class AmrObject {
    
public:
    /// Methods for tagging cells for refinement
    enum TaggingCriteria {
        kChargeDensity = 0, // default
        kPotentialStrength,
        kEfieldGradient
    };
    
    virtual void initialize() {}
    
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
