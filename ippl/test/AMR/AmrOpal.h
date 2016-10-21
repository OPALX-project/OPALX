#ifndef AMROPAL_H
#define AMROPAL_H

#include <AmrCore.H>

class AmrOpal : public AmrCore {
    
public:
    
    AmrOpal();
    virtual ~AmrOpal();
    
    void ErrorEst(int lev, MultiFab& mf, TagBoxArray& tags, Real time);

protected:
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;
    
};

#endif