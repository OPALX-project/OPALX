
#include "Ippl.h"

#include "AmrParticleBase1.h"
#include "ParticleAmrLayout.h"
#include "PartBunchAmr.h"

#define Dim 3

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase1<amrplayout_t> amrbase_t;
typedef PartBunchAmr<amrplayout_t> amrbunch_t;


template <class PLayout>
class BoxLibLayout : public ParticleAmrLayout

class AbstractBunch {
    
    
    
};

template <class Bunch>
class PartBunchBase : public AbstractBunch
{
public:
    
    
private:
    Bunch* bunch_mp;
};

class PartBunch : public PartBunchBase<PartBunch>,
                  public IpplParticleBase<ParticleSpatialLayout<double, 3> >
{
public:
    
    
    
};

class AmrPartBunch : public PartBunchBase<AmrPartBunch>,
                     public AmrParticleBase1<BoxLibLayout<double, 3> >
{
public:
    
    
};


void doIppl(Array<Geometry> &geom, Array<BoxArray> &ba, 
            Array<DistributionMapping> &dmap, Array<int> &rr, 
            size_t nLevels)
{
    amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);
    
    //create a particle bunch
    PartBunchAmr<amrplayout_t>* pbase = new PartBunchAmr<amrplayout_t>();
    pbase->initialize(PL);
    pbase->initializeAmr();
    
    pbase->update();
    
}


int main(int argc, char** argv) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc,argv, false);
    
    
    size_t nLevels = 2;
    size_t maxBoxSize = 8;
    
    //set up the geometry
    int n_cell = 16;
    IntVect low(0, 0, 0);
    IntVect high(n_cell - 1, n_cell - 1, n_cell - 1);    
    Box bx(low, high);
    
    //physical domain boundaries
    RealBox domain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        domain.setLo(i, 0.0);
        domain.setHi(i, 1.0);
    }
    
    RealBox fine_domain;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        fine_domain.setLo(i, 0.0);
        fine_domain.setHi(i, 0.5);
    }
    
    //periodic boundary conditions in all directions
    int bc[AMREX_SPACEDIM] = {1, 1, 1};
    
    //Container for geometry at all levels
    Array<Geometry> geom;
    geom.resize(nLevels);
    
    // Container for boxes at all levels
    Array<BoxArray> ba;
    ba.resize(nLevels);    
    
    // level 0 describes physical domain
    geom[0].define(bx, &domain, 0, bc);
    
    //refinement for each level
    Array<int> rr(nLevels - 1);
    for (unsigned int lev = 0; lev < rr.size(); ++lev)
        rr[lev] = 2;
    
    // geometries of refined levels
    for (unsigned int lev = 1; lev < nLevels; ++lev)
        geom[lev].define(BoxLib::refine(geom[lev - 1].Domain(), rr[lev - 1]), &domain, 0, bc);
    
    // box at level 0
    ba[0].define(bx);
    ba[0].maxSize(maxBoxSize);

    //box at level 1
    
    //build boxes at finer levels
    if (nLevels > 1) {
        int n_fine = n_cell * rr[0];
        IntVect refined_lo(0, 0, 0);
        IntVect refined_hi(15, 15, 15);
        
        Box refined_box(refined_lo, refined_hi);
        ba[1].define(refined_box);
        ba[1].maxSize(maxBoxSize);
    }

    /*
     * distribution mapping
     */
    Array<DistributionMapping> dmap;
    dmap.resize(nLevels);
    dmap[0].define(ba[0], ParallelDescriptor::NProcs() /*nprocs*/);
    if (nLevels > 1)
    dmap[1].define(ba[1], ParallelDescriptor::NProcs() /*nprocs*/);
    
    doIppl(geom, ba, dmap, rr, nLevels);
    
    return 0;
}