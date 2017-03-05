/*!
 * @file testDeposition.cpp
 * @details Fully BoxLib charge deposition timing. It initializes
 *          Gaussian distribution and refines the centere eighth.
 * @authors Andrew Myers
 * @date February 2017
 * @brief Charge deposition timing
 */

#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <Particles.H>
#include <PlotFileUtil.H>

#include <chrono>

typedef std::chrono::time_point<std::chrono::high_resolution_clock> chrono_t;

#include <random>

class MyParticleContainer
    : public ParticleContainer<1+2*BL_SPACEDIM>
{
 public:

    MyParticleContainer(const Array<Geometry>            & geom, 
                        const Array<DistributionMapping> & dmap,
                        const Array<BoxArray>            & ba,
                        const Array<int>                 & rr)
	: ParticleContainer<1+2*BL_SPACEDIM> (geom, dmap, ba, rr)
        {
        }

    void InitParticles(unsigned long icount,
                       unsigned long iseed,
                       Real          mass)
    {
        BL_PROFILE("MyParticleContainer::InitParticles()");
        BL_ASSERT(m_gdb != 0);
        
        const int       MyProc   = ParallelDescriptor::MyProc();
        const int       NProcs   = ParallelDescriptor::NProcs();
        const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
        const Real      strttime = ParallelDescriptor::second();
        const Geometry& geom     = m_gdb->Geom(0);
        
        Real r, x, len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                               geom.ProbLength(1),
                                               geom.ProbLength(2)) };
        
        RealBox containing_bx = geom.ProbDomain();
        const Real* xlo = containing_bx.lo();
        const Real* xhi = containing_bx.hi();
        
        m_particles.resize(m_gdb->finestLevel()+1);
        
        for (int lev = 0; lev < m_particles.size(); lev++)
            BL_ASSERT(m_particles[lev].empty());
        
        std::mt19937_64 mt(0/*seed*/ /*0*/);
        std::normal_distribution<double> dist(0.0, 0.07);
        
        // We'll let IOProc generate the particles so we get the same
        // positions no matter how many CPUs we have.  This is here
        // mainly for debugging purposes.  It's not really useful for
        // very large numbers of particles.
        Array<ParticleBase::RealType> pos(icount*BL_SPACEDIM);
        
        if (ParallelDescriptor::IOProcessor()) {
            for (long j = 0; j < icount; j++) {
                for (int i = 0; i < BL_SPACEDIM; i++) {
                    do {
                        x = dist(mt) + 0.5;
                    }
                    while (x < xlo[i] || x > xhi[i]);
                    
                    pos[j*BL_SPACEDIM + i] = x;
                }
            }
        }
        
        // Send the particle positions to other procs (this is very slow)
        ParallelDescriptor::Bcast(pos.dataPtr(), icount*BL_SPACEDIM, IOProc);
        
        int cnt = 0;
        for (long j = 0; j < icount; j++) {
            ParticleType p;
            
            for (int i = 0; i < BL_SPACEDIM; i++)
                p.m_pos[i] = pos[j*BL_SPACEDIM + i];
        
            p.m_data[0] = mass;
            
            for (int i = 1; i < 1 + 2*BL_SPACEDIM; i++)
                p.m_data[i] = 0;

            if (!ParticleBase::Where(p,m_gdb))
                BoxLib::Abort("invalid particle");
        
            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
            
            const int who = m_gdb->ParticleDistributionMap(p.m_lev)[p.m_grid];
            
            if (who == MyProc) {
                p.m_id  = ParticleBase::NextID();
                p.m_cpu = MyProc;
                
                m_particles[p.m_lev][p.m_grid].push_back(p);
                
                cnt++;
            }
        }
    
        BL_ASSERT(OK());
    }
};

struct TestParams {
    int nx;
    int ny;
    int nz;
    int max_grid_size;
    int nppc;
    int nlevs;
    bool verbose;
};

void test_assign_density(TestParams& parms)
{
    
    int nlevs = parms.nlevs;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
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
        IntVect refined_lo(3*n_fine/8,3*n_fine/8,3*n_fine/8); 
        IntVect refined_hi(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Array<DistributionMapping> dmap(nlevs);
    PArray<MultiFab> density;
    PArray<MultiFab> partMF;  

    density.resize(nlevs,PArrayManage);
    partMF.resize(nlevs,PArrayManage);
    for (int lev = 0; lev < nlevs; lev++) {
	density.set     (lev,new MultiFab(ba[lev],1          ,0));
	partMF.set      (lev,new MultiFab(ba[lev],1          ,1));

	density[lev].setVal(0.0);
	partMF[lev].setVal(0.0, 1);

        dmap[lev] = density[lev].DistributionMap();
    }

    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * parms.nx * parms.ny * parms.nz;
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    myPC.InitParticles(num_particles, iseed, mass);
    
    chrono_t start, end;
    if ( parms.verbose && ParallelDescriptor::IOProcessor() )
        start = std::chrono::high_resolution_clock::now();

    myPC.AssignDensity(0, false, partMF, 0, 1, 1);
    
    if ( parms.verbose && ParallelDescriptor::IOProcessor() ) {
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << "Time for AssignDensity: " << diff.count()
                  << " s" << std::endl;
    }
}

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nppc", parms.nppc);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    BoxLib::Abort("Must specify at least one particle per cell");
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  pp.get("nlevs", parms.nlevs);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }
  
  test_assign_density(parms);
  
  BoxLib::Finalize();
}
