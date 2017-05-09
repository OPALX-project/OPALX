#ifndef AMR_DEFS_H
#define AMR_DEFS_H

#include <MultiFab.H>
#include <Geometry.H>
#include <memory>

namespace amr {
    typedef MultiFab                                AmrField_t;
    typedef Array< std::unique_ptr<AmrField_t>  >   AmrFieldContainer_t;
    typedef Array<Geometry>                         AmrGeomContainer_t;
    typedef Array<BoxArray>                         AmrGridContainer_t;
    typedef Array<DistributionMapping>              AmrProcMapContainer_t;
    typedef RealBox                                 AmrDomain_t;
    typedef Array<int>                              AmrIntArray_t;
};

#endif
