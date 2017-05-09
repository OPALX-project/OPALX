/*!
 * @file testScatterAMReX.cpp
 * @author Matthias Frey, Andrew Myers
 * 
 * @details Compare the fields of AMReX with OPAL. We use the
 * test example of the AMReX repository
 * amrex/Tests/Particles/AssignMultiLevelDensity
 * 
 * https://Yavin_4@bitbucket.org/berkeleylab/amrex.git
 */
#include <iostream>

#include "Ippl.h"

#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>

#include "AmrParticleBase.h"
#include "ParticleAmrLayout.h"
#include "PartBunchAmr.h"

#define Dim 3

typedef ParticleAmrLayout<double,Dim> amrplayout_t;
typedef AmrParticleBase<amrplayout_t> amrbase_t;
typedef PartBunchAmr<amrplayout_t> amrbunch_t;

typedef std::deque<Particle<1,0> > PBox;
typedef typename std::map<int,PBox> PMap;

typedef Array<std::unique_ptr<MultiFab> > container_t;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
};


void createParticles(TestParams& parms,
                     PartBunchAmr<amrplayout_t>* pbase,
                     const Array<Geometry>& geom,
                     const Array<DistributionMapping>& dmap,
                     const Array<BoxArray>& ba,
                     const Array<int>& rr)
{
    typedef ParticleContainer<1> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);
    
    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    //    myPC.InitRandom(num_particles, iseed, mass, serialize, fine_box);
    myPC.InitRandom(num_particles, iseed, mass, serialize);
//     std::string distr = "amrex_particle_distribution.ascii";
//     myPC.InitFromAsciiFile(distr, 0, 0);
    
//     myPC.WriteAsciiFile("new_boxlib_particle_distribution.ascii");
    
    int nLocParticles = myPC.TotalNumberOfParticles(true, true);
    int nGlobParticles = myPC.TotalNumberOfParticles();
    std::cout << "#Local Particles: " << nLocParticles << std::endl;
    std::cout << "#Global Particles: " << nGlobParticles << std::endl;
    pbase->create(nLocParticles);
    
    int i = 0;
    for (unsigned int lev = 0; lev < dmap.size(); lev++) {
        PMap pmap = myPC.GetParticles(lev);
        for (unsigned int grid = 0; grid < pmap.size(); ++grid) {
            for (unsigned int p = 0; p < pmap[grid].size(); ++p) {
                pbase->qm[i] = 10.0;//pmap[grid][p].m_data[0];
                pbase->R[i](0) = pmap[grid][p].m_pos[0];
                pbase->R[i](1) = pmap[grid][p].m_pos[1];
                pbase->R[i++](2) = pmap[grid][p].m_pos[2];
            }
        }
    }
}

void doTestScatter(TestParams& parms) {
    
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);
    
    int nlevs = parms.nlevs;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Array<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(BoxLib::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, coord, is_per);
    }

    Array<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(n_fine/4,n_fine/4,n_fine/4); 
        IntVect refined_hi(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Array<DistributionMapping> dmap(nlevs);

    Array<std::unique_ptr<MultiFab> > partMF(nlevs);
//     Array<std::unique_ptr<MultiFab> > partMF_old(nlevs);
//     Array<std::unique_ptr<MultiFab> > density(nlevs);
//     Array<std::unique_ptr<MultiFab> > acceleration(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping(ba[lev], ParallelDescriptor::NProcs());
        partMF[lev].reset(new MultiFab(ba[lev], 1, 2, dmap[lev]));
        partMF[lev]->setVal(0.0, 2);
        
//         partMF_old[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 2));
//         partMF_old[lev]->setVal(0.0, 2);
//         density[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 0));
//         density[lev]->setVal(0.0);
//         acceleration[lev].reset(new MultiFab(ba[lev], dmap[lev], 3, 1));
//         acceleration[lev]->setVal(5.0, 1);
    }
    
    
    
    
    //create a new layout using ParticleAmrLayout class
    amrplayout_t* PL = new amrplayout_t(geom, dmap, ba, rr);
    
    //create a particle bunch
    PartBunchAmr<amrplayout_t>* pbase = new PartBunchAmr<amrplayout_t>();
    pbase->initialize(PL);
    pbase->initializeAmr();
    
    createParticles(parms, pbase, geom, dmap, ba, rr);
    
    //update redistributes particles among the cores
    
    
    //call assign density to scatter the paarticle attribute qm on the grid
//     pbase->setAllowParticlesNearBoundary(true);
    pbase->update();
    
    std::cout << "#Local Particles = " << pbase->getLocalNum() << std::endl;
    std::cout << "#Particles = " << pbase->getTotalNum() << std::endl;
    
//     pbase->AssignCellDensitySingleLevelFort(pbase->qm, *(partMF[0].get()), 0);
    pbase->AssignDensityFort(pbase->qm, partMF, 0, 1, nlevs-1);
    
    // single core
    unsigned int lev = 0, nLocParticles = 0;
    for (unsigned int ip = 0; ip < pbase->getLocalNum(); ++ip) {
        if ( lev != pbase->m_lev[ip] ) {
            std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
            nLocParticles = 0;
            lev++;
        }
        nLocParticles++;
    }
    
    if ( lev == 0 )
        std::cout << "#Local Particles at level " << lev << ": " << nLocParticles << std::endl;
    
    for (unsigned int lev = 0; lev < partMF.size(); ++lev) {
        //calculate the sum of all the components in multifab
        double sum = partMF[lev]->sum();
    
        std::cout << "sum = " << sum << std::endl;
    }
    
//     delete pbase;
    
    IpplTimings::stopTimer(mainTimer);
}


int main(int argc, char* argv[]) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc,argv/*, false*/);
    
    ParmParse pp;
  
    TestParams parms;
    
    pp.get("nx", parms.nx);
    pp.get("ny", parms.ny);
    pp.get("nz", parms.nz);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nlevs", parms.nlevs);
    pp.get("nppc", parms.nppc);
    if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("Must specify at least one particle per cell");
    
    parms.verbose = false;
    pp.query("verbose", parms.verbose);
    
    if (parms.verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl;
        std::cout << "Number of particles per cell : ";
        std::cout << parms.nppc  << std::endl;
        std::cout << "Size of domain               : ";
        std::cout << "Num levels: ";
        std::cout << parms.nlevs << std::endl;
        std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
    }
    
    doTestScatter(parms);
    
    IpplTimings::print();
    
    return 0;
}