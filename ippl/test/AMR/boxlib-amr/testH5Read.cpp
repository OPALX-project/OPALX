#include <iostream>

#include "Ippl.h"
#include <memory>
#include <MultiFab.H>

#include "H5Reader.h"

#include <cassert>

// #include "Distribution.h"

int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc, argv, false);
    
    assert(BL_SPACEDIM == 3);
    
    int nr = 8;
    
    IntVect low(0, 0, 0);
    IntVect high(nr - 1, nr - 1, nr - 1);    
    Box bx(low, high);
    
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    int bc[BL_SPACEDIM] = {0, 0, 0};
    
    Array<Geometry> geom;
    geom.resize(1);
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    Array<BoxArray> ba;
    
    ba.resize(1);
    
    ba[0].define(bx);
    ba[0].maxSize(nr);
    
    Array<DistributionMapping> dmap;
    dmap.resize(1);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    
    container_t rhs(1);
    rhs.set(0, new MultiFab(ba[0], 1, 0));
    rhs[0].setVal(1.0);
    
    H5Reader h5reader;
    
    h5reader.writeScalarField(rhs, geom);
    
    rhs.clear();
    
    return 0;
}