#ifndef AMROPAL_H
#define AMROPAL_H

#include <AmrCore.H>
#include "AmrPartBunch.h"

class AmrOpal : public AmrCore {
    
public:
    
//     AmrOpal();
    AmrOpal(const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord, PartBunchBase* bunch);

    virtual ~AmrOpal();
    
    // defines BoxArray, DistributionMapping of level 0
    void initBaseLevel();
    
//     void initFineLevels();
    
    void regrid (int lbase, Real time);
    
    void RemakeLevel (int lev, Real time,
                      const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap);
    
    
    void writePlotFile();

protected:
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;
    
private:
    AmrPartBunch* bunch_m;
    PArray<MultiFab> nPartPerCell_m;
    
};

#endif