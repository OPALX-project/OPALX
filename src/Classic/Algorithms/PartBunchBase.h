#ifndef PartBunchBase_H
#define PartBunchBase_H

class PartBunchBase {
    
public:
    
    virtual void computeSelfFields() = 0;

    /** /brief used for self fields with binned distribution */
    virtual void computeSelfFields(int b) = 0;

    virtual void computeSelfFields_cycl(double gamma) = 0;
    virtual void computeSelfFields_cycl(int b) = 0;
    
    
};

#endif