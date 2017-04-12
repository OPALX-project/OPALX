#include "Ippl.h"

#include "PartBunch.h"
#include "AmrPartBunch.h"

// #include "BoxLibClasses.h"

#include <memory>


void createRandomParticles(PartBunchBase<double, 3>* pbase, int N, int myNode, int seed = 1) {
    srand(seed);
    for (int i = 0; i < N; ++i) {
        pbase->createWithID(myNode * N + i + 1);
        pbase->Q[i] = (double)rand() / RAND_MAX;
    
        pbase->R[i][0] = (double)rand() / RAND_MAX;
        pbase->R[i][1] = (double)rand() / RAND_MAX;
        pbase->R[i][2] = (double)rand() / RAND_MAX;  
    }
}

void initBunch(std::unique_ptr<PartBunchBase<double, 3> >& bunch) {
    
    size_t nLevels = 2;
    size_t maxBoxSize = 8;
    
    //set up the geometry
    int n_cell = 32;
    IntVect low(0, 0, 0);
    IntVect high(n_cell - 1, n_cell - 1, n_cell - 1);    
    Box bx(low, high);
    
    //physical domain boundaries
    RealBox domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        domain.setLo(i, -2.0);
        domain.setHi(i, 2.0);
    }
    
    RealBox fine_domain;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        fine_domain.setLo(i, 0.0);
        fine_domain.setHi(i, 0.5);
    }
    
    //periodic boundary conditions in all directions
    int bc[BL_SPACEDIM] = {1, 1, 1};
    
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
    
    bunch.reset(
        new AmrPartBunch(
            new BoxLibParticle<BoxLibLayout<double, 3> >(
                new BoxLibLayout<double, 3>(geom, dmap, ba, rr)
            )
        )
    );
}



int main(int argc, char** argv) {
    
    Ippl ippl(argc, argv);
    
    std::cout << "Init PartBunch" << std::endl
              << "--------------" << std::endl;
    
    std::unique_ptr<PartBunchBase<double, 3> > bunch(
        new PartBunch(
            new IpplParticleBase<ParticleSpatialLayout<double, 3> >(
                new ParticleSpatialLayout<double, 3>()
            )
        )
    );
    
//     std::unique_ptr<PartBunchBase<double, 3> > bunch(new PartBunch(new IpplParticleBase<ParticleSpatialLayout<double, 3> >()));
    
    
    
    
//     dynamic_cast<PartBunch::pbase_t*>(bunch->pbase)->initialize(new PartBunch::pbase_t::Layout_t());
    
//     bunch->R = *bunch->pbase->R_p;
//     bunch->ID = *bunch->pbase->ID_p;
    
    createRandomParticles(bunch.get(), 10, Ippl::myNode());
    
//     FieldLayout_t fl = bunch->getFieldLayout();
    
    
    
    e_dim_tag decomp[3] = {SERIAL, SERIAL, SERIAL};

    NDIndex<3> domain;
    domain[0] = Index(16 + 1);
    domain[1] = Index(16 + 1);
    domain[2] = Index(16 + 1);

    
        decomp[0] = PARALLEL;
    
//         decomp[1] = PARALLEL;
    
//         decomp[2] = PARALLEL;
    // create prototype mesh and layout objects for this problem domain
    Mesh_t* mesh_m   = new Mesh_t(domain);
    FieldLayout_t* FL_m     = new FieldLayout_t(*mesh_m, decomp);
//     PL_m     = new Layout_t(*FL_m, *mesh_m);
//     FieldLayout_t fl =
    bunch->setFieldLayout(*FL_m);
    
//     bunch.reset(
//         new PartBunch(
//             new IpplParticleBase<ParticleSpatialLayout<double, 3> >(
//                 new ParticleSpatialLayout<double, 3>(*FL_m, *mesh_m)
//             )
//         )
//     );
    
    std::cout << "Hallo" << std::endl;
    
    
    PartBunch::pbase_t::Layout_t* layout = static_cast<PartBunch::pbase_t::Layout_t*>(&bunch->getLayout());
    
    RegionLayout<double ,3,Mesh_t>& RL = layout->getLayout();
//     RegionLayout<double ,3,Mesh_t>& RL = dynamic_cast<PartBunch::pbase_t*>(bunch.get())->getLayout().getLayout();
    
    if ( ! RL.initialized()) {
        ERRORMSG("Cannot repartition particles: uninitialized layout." << endl);
        return false;
    }
    
    FieldLayout<3>& FL = RL.getFieldLayout();
    Mesh_t& mesh = RL.getMesh();
    
    const NDIndex<3>& TotalDomain = FL.getDomain();
//     NDIndex<Dim> indx;
    
    bool CenterOffset[3];
    int CenteringTotal = 0;
    unsigned int d;
    for (d=0; d<3; ++d) {
        CenterOffset[d] = (TotalDomain[d].length() < mesh.gridSizes[d]);
        CenteringTotal += CenterOffset[d];
    }
    
    
    std::cout << "Hallo" << std::endl;
    
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
//     bunch->update();
    
//     bunch->getFieldLayout();
    
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
//     
//     std::cout << "bunch->R[0] = " << bunch->R[0] << std::endl;
//     
//     double two = 2.0;
//     bunch->R *= Vector_t(two);
//     
//     std::cout << "2 * bunch->R[0] = " << bunch->R[0] << std::endl;
//     std::cout << "max(bunch->R) = " << max(bunch->R) << std::endl;
//     
//     std::cout << "Destroy R[0]" << std::endl;
//     bunch->destroy(1, 0);
//     bunch->update();
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
//     std::cout << std::endl << "Init AmrPartBunch" << std::endl
//               << "-----------------" << std::endl;
//     
//     BoxLib::Initialize(argc,argv, false);
//     
//     bunch.reset(nullptr);
//     initBunch(bunch);
//     
//     createRandomParticles(bunch.get(), 10, Ippl::myNode());
//     
//     
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
//     
//     bunch->update();
//     
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
//     
//     std::cout << "bunch->R[0] = " << bunch->R[0] << std::endl;
//     
//     bunch->R *= Vector_t(2.0);
//     
//     std::cout << "2 * bunch->R[0] = " << bunch->R[0] << std::endl;
//     std::cout << "max(bunch->R) = " << max(bunch->R) << std::endl;
//     
//     bunch->update();
//     
//     std::cout << "Destroy R[0]" << std::endl;
//     std::cout << "level[0] = " << dynamic_cast<BoxLibParticle<BoxLibLayout<double, 3> >*>(bunch->pbase)->level[0] << std::endl;
//     std::cout << "grid[0] = " << dynamic_cast<BoxLibParticle<BoxLibLayout<double, 3> >*>(bunch->pbase)->grid[0] << std::endl;
//     bunch->destroy(1, 0);
//     bunch->update();
//     std::cout << bunch->getTotalNum() << " " << bunch->getLocalNum() << std::endl;
    
    return 0;
}