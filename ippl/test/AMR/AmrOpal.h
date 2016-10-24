#ifndef AMROPAL_H
#define AMROPAL_H

#include <AmrCore.H>
#include "AmrPartBunch.h"

#include <MultiFabUtil.H>

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
    
    
    void info() {
        for (int i = 0; i < finest_level; ++i)
            std::cout << "density level " << i << ": "
                      << nPartPerCell_m[i].min(0) << " "
                      << nPartPerCell_m[i].max(0) << std::endl;
    }
    
    void writePlotFile(std::string filename, int step);

protected:
    virtual void ErrorEst(int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;
    
private:
    AmrPartBunch* bunch_m;
    PArray<MultiFab> nPartPerCell_m;
    
};

#endif