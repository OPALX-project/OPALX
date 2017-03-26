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
    
    virtual void regrid() = 0;
    
    void setTagging(TaggingCriteria tagging) {
        tagging_m = tagging;
    }
    
protected:
    TaggingCriteria tagging_m;
    
};

#endif
