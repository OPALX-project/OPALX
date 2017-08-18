/*!
 * @file testGeometry.cpp
 * @author Matthias Frey
 * @date May 2017
 * @details Init particles inside the domain [0, 0.5]^3 and the update
 * the geometry to [0.15, 0.35]^3. Some particles will be lost due to
 * the geometry change.
 */
#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>

#include <random>

/* NStructReal = 3 (mass, position are already reserved)
 * 
 * p.m_rdata.arr[0] = pos x
 * p.m_rdata.arr[1] = pos y
 * p.m_rdata.arr[2] = pos z
 * p.m_rdata.arr[3] = mass
 * p.m_rdata.arr[4] = momenta x
 * p.m_rdata.arr[5] = momenta y
 * p.m_rdata.arr[6] = momenta z
 */

using namespace amrex;

class MyParticleContainer
    : public ParticleContainer<BL_SPACEDIM>
{
public:

    MyParticleContainer(const Array<Geometry>            & geom, 
                        const Array<DistributionMapping> & dmap,
                        const Array<BoxArray>            & ba,
                        const Array<int>                 & rr)
        : ParticleContainer<BL_SPACEDIM> (geom, dmap, ba, rr), particles_rm(GetParticles())
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
        
        particles_rm.resize(m_gdb->finestLevel()+1);
        
        for (int lev = 0; lev < particles_rm.size(); lev++)
            BL_ASSERT(particles_rm[lev].empty());
        
        std::mt19937_64 mt(0/*seed*/ /*0*/);
        std::normal_distribution<double> dist(0.0, 0.07);
        
        // We'll let IOProc generate the particles so we get the same
        // positions no matter how many CPUs we have.  This is here
        // mainly for debugging purposes.  It's not really useful for
        // very large numbers of particles.
        Array<typename ParticleType::RealType> pos(icount*BL_SPACEDIM);
        
        if (ParallelDescriptor::IOProcessor()) {
            for (long j = 0; j < icount; j++) {
                for (int i = 0; i < BL_SPACEDIM; i++) {
                    do {
                        x = dist(mt) + 0.25;
                    }
                    while (x < xlo[i] || x > xhi[i]);
                    
                    pos[j*BL_SPACEDIM + i] = x;
                }
            }
        }
        
        // Send the particle positions to other procs (this is very slow)
        ParallelDescriptor::Bcast(pos.dataPtr(), icount*BL_SPACEDIM, IOProc);
        
        ParticleLocData pld;
        
        int cnt = 0;
        for (long j = 0; j < icount; j++) {
            ParticleType p;
            
            for (int i = 0; i < BL_SPACEDIM; i++)
                p.m_rdata.arr[i] = pos[j*BL_SPACEDIM + i];
        
            p.m_rdata.arr[BL_SPACEDIM] = mass;
            
            // momenta
//             for (int i = 1; i < 1 + BL_SPACEDIM; i++)
                p.m_rdata.arr[4/*BL_SPACEDIM + i*/] = 0.0;
                p.m_rdata.arr[5/*BL_SPACEDIM + i*/] = 0.0;
                p.m_rdata.arr[6/*BL_SPACEDIM + i*/] = 0.0;

            if (!this->Where(p, pld))
                amrex::Abort("invalid particle");
        
            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= m_gdb->finestLevel());
            
            const int who = ParticleDistributionMap(pld.m_lev)[pld.m_grid];
            
            if (who == MyProc) {
                p.m_idata.id  = ParticleType::NextID();
                p.m_idata.cpu = MyProc;
                
                particles_rm[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
                
                cnt++;
            }
        }
    
        BL_ASSERT(OK());
    }
    
    void printParticles() {
        std::cout << "Geometry:" << std::endl;
        Geometry geom = m_gdb->Geom(0);
        
        RealBox rb = geom.ProbDomain();
        for (int n = 0; n < BL_SPACEDIM; n++) {
            std::cout << rb.lo(n) << " " << rb.hi(n) << std::endl;
        }
        
        for (int lev = 0; lev < particles_rm.size();  lev++) {
            std::cout << "level = " << lev << std::endl;
            auto& pmap = particles_rm[lev];
            for (const auto& kv : pmap) {
                const auto& aos = kv.second.GetArrayOfStructs();
                
                for (auto it = aos.cbegin(); it != aos.cend(); ++it) {
                    if (it->m_idata.id > 0) {

                        // write out the particle struct first... 
                        std::cout << it->m_rdata.pos[0] << ' '
                                  << it->m_rdata.pos[1] << ' '
                                  << it->m_rdata.pos[2] << ' '
                                  << std::endl;	     
                    }
                }
            }
        }
    }
    
    void updateGeometry(double lower, double upper) {
        RealBox real_box;
        for (int n = 0; n < BL_SPACEDIM; n++) {
            real_box.setLo(n, lower);
            real_box.setHi(n, upper);
        }
        
        Geometry::ProbDomain(real_box);
        
        // This says we are using Cartesian coordinates
        int coord = 0;
    
        // This sets the boundary conditions to be doubly or triply periodic
        int is_per[BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM; i++) 
            is_per[i] = 0; 
        
        int nlevs = GetParticles().size();
        int ngrids = 32;
        int max_grid_size = 8;
        
        // Define the refinement ratio
        Array<int> rr(nlevs-1);
        for (int lev = 1; lev < nlevs; lev++)
            rr[lev-1] = 2;
        
        IntVect domain_lo(0 , 0, 0); 
        IntVect domain_hi(ngrids - 1, ngrids - 1, ngrids - 1); 
        const Box domain(domain_lo, domain_hi);
        
        // This defines a Geometry object which is useful for writing the plotfiles  
        Array<Geometry> geom(nlevs);
        geom[0].define(domain, &real_box, coord, is_per);
        for (int lev = 1; lev < nlevs; lev++) {
            geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                            &real_box, coord, is_per);
        }
    
        Array<BoxArray> ba(nlevs);
        ba[0].define(domain);
        
        // Now we make the refined level be the center eighth of the domain
        if (nlevs > 1) {
            int n_fine = ngrids*rr[0];
            IntVect refined_lo(3*n_fine/8,3*n_fine/8,3*n_fine/8); 
            IntVect refined_hi(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1);
    
            // Build a box for the level 1 domain
            Box refined_patch(refined_lo, refined_hi);
            ba[1].define(refined_patch);
        }
        
        // break the BoxArrays at both levels into max_grid_size^3 boxes
        Array<DistributionMapping> dmap(nlevs);
        for (int lev = 0; lev < nlevs; lev++) {
            ba[lev].maxSize(max_grid_size);
            dmap[lev].define(ba[lev], ParallelDescriptor::NProcs() /*nprocs*/);
        }
        
        this->Define(geom, dmap, ba, rr);
    }
    
    
private:
    Array<ParticleLevel>& particles_rm;
};

void testGeometry() {
    
    int nlevs = 2;
    int ngrids = 32;
    int num_particles = 10;
    int max_grid_size = 8;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 0.5);
    }

    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(ngrids - 1, ngrids - 1, ngrids - 1); 
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
        is_per[i] = 0; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, coord, is_per);
    }

    Array<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = ngrids*rr[0];
        IntVect refined_lo(3*n_fine/8,3*n_fine/8,3*n_fine/8); 
        IntVect refined_hi(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1);

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    Array<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(max_grid_size);
        dmap[lev].define(ba[lev]);
    }
    
    
    
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);
    
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    myPC.InitParticles(num_particles, iseed, mass);
    
    myPC.printParticles();
    
    myPC.Redistribute();
    
    myPC.updateGeometry(0.15, 0.35);
    
    myPC.Redistribute();
    
    myPC.printParticles();
    
}

int main(int argc, char* argv[]) {
    
    amrex::Initialize(argc,argv);
  
    testGeometry();
  
    amrex::Finalize();
    
    return 0;
}