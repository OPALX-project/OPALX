/*!
 * @file testDeposition.cpp
 * @details Fully AMReX charge deposition timing. It initializes
 *          Gaussian distribution and refines the centere eighth.
 * @authors Andrew Myers
 * @date February 2017
 * @brief Charge deposition timing
 */

#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>

#include <AMReX_AmrMesh.H>
#include "../AmrOpal_F.h"

#include <random>

using namespace amrex;


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

class MyParticleContainer
    : public ParticleContainer<BL_SPACEDIM /* momenta */>
{
 public:

    MyParticleContainer(const Geometry            & geom, 
                        const DistributionMapping & dmap,
                        const BoxArray            & ba)
        : ParticleContainer<BL_SPACEDIM> (geom, dmap, ba), particles_rm(GetParticles())
        {}

    // init particles within sphere of radius r
    void initSphere(unsigned long nParticles,
                    double r)
    {
        BL_PROFILE("MyParticleContainer::initSphere()");
        BL_ASSERT(m_gdb != 0);
        
        const int       MyProc   = ParallelDescriptor::MyProc();
        const int       NProcs   = ParallelDescriptor::NProcs();
        const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
        const Real      strttime = ParallelDescriptor::second();
        
        particles_rm.resize(m_gdb->finestLevel()+1);
        
        for (int lev = 0; lev < particles_rm.size(); lev++)
            BL_ASSERT(particles_rm[lev].empty());
        
        int nLocParticles = nParticles / NProcs;
        
        std::mt19937_64 eng[3];
    
        for (int i = 0; i < 3; ++i) {
            eng[i].seed(42 + 3 * i);
            eng[i].discard( nLocParticles * MyProc);
        }
        
        std::uniform_real_distribution<> ph(-1.0, 1.0);
        std::uniform_real_distribution<> th(0.0, 2.0 * M_PI);
        std::uniform_real_distribution<> u(0.0, 1.0);
    
        long double qi = 4.0 * M_PI * r * r / double(nParticles);
        
        ParticleLocData pld;
        
        for (uint i = 0; i < nLocParticles; ++i) {
            double phi = std::acos( ph(eng[0]) );
            double theta = th(eng[1]);
            double radius = r * std::cbrt( u(eng[2]) );
        
            double x = radius * std::cos( theta ) * std::sin( phi );
            double y = radius * std::sin( theta ) * std::sin( phi );
            double z = radius * std::cos( phi );
        
            ParticleType p;
            
            p.m_rdata.arr[0] = x;
            p.m_rdata.arr[1] = y;
            p.m_rdata.arr[2] = z;
            
            // charge instead of mass
            p.m_rdata.arr[BL_SPACEDIM] = qi; // C
            
            // momenta
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
            }
        }
    
        BL_ASSERT(OK());
    }
    
private:
    Array<ParticleLevel>& particles_rm;
};


class MyAmr : public AmrMesh {
    
public:
    MyAmr(const RealBox* rb,
          int max_level_in,
          const Array<int>& n_cell_in,
          int coord,
          MyParticleContainer* bunch)
        : AmrMesh(rb, max_level_in, n_cell_in, coord),
          bunch_mp(bunch),
          nCharge_m(1.0e-9)
    {
        // setup base level
        this->finest_level = 0;
        
        const ParGDBBase* gdb = bunch->GetParGDB();
        
        const BoxArray& ba = gdb->ParticleBoxArray(0);
        const DistributionMapping& dm = gdb->ParticleDistributionMap(0);
    
        this->SetBoxArray(0, ba);
        this->SetDistributionMap(0, dm);
        
        nChargePerCell_m.resize(max_level_in + 1);
        nChargePerCell_m[0] = std::unique_ptr<MultiFab>(
            new MultiFab(this->boxArray(0),
                         this->DistributionMap(0),
                         1, 1)
                                                       );
        nChargePerCell_m[0]->setVal(0.0, 1);
    }
    
    
    void regrid (int lbase, Real time)
    {
        int new_finest = 0;
        Array<BoxArray> new_grids(finest_level+2);
        
        MakeNewGrids(lbase, time, new_finest, new_grids);
    
        BL_ASSERT(new_finest <= finest_level+1);
        
        for (int lev = lbase+1; lev <= new_finest; ++lev)
        {
            if (lev <= finest_level) // an old level
            {
                if (new_grids[lev] != grids[lev]) // otherwise nothing
                {
                    DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
                    RemakeLevel(lev, time, new_grids[lev], new_dmap);
                }
            }
            else  // a new level
            {
                DistributionMapping new_dmap(new_grids[lev], ParallelDescriptor::NProcs());
                MakeNewLevel(lev, time, new_grids[lev], new_dmap);
            }
        }
        
        for (int lev = new_finest+1; lev <= finest_level; ++lev) {
            ClearLevel(lev);
            ClearBoxArray(lev);
            ClearDistributionMap(lev);
        }
        
        finest_level = new_finest;
        
        // update to multilevel
        bunch_mp->Redistribute();
    }
    
    void ErrorEst(int lev, TagBoxArray& tags, Real time, int ngrow) {
        for (int i = lev; i <= finest_level; ++i) {
            nChargePerCell_m[i]->setVal(0.0, 1);
        }
        
        bunch_mp->AssignDensityFort(0, nChargePerCell_m, lev, 1, finest_level);
        
        const int clearval = TagBox::CLEAR;
        const int   tagval = TagBox::SET;
    
        const Real* dx      = geom[lev].CellSize();
        const Real* prob_lo = geom[lev].ProbLo();
        
        #ifdef _OPENMP
        #pragma omp parallel
        #endif
        {
            Array<int>  itags;
            for (MFIter mfi(*nChargePerCell_m[lev],false/*true*/); mfi.isValid(); ++mfi) {
                const Box&  tilebx  = mfi.validbox();//mfi.tilebox();
                
                TagBox&     tagfab  = tags[mfi];
                
                // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
                // So we are going to get a temporary integer array.
                tagfab.get_itags(itags, tilebx);
                
                // data pointer and index space
                int*        tptr    = itags.dataPtr();
                const int*  tlo     = tilebx.loVect();
                const int*  thi     = tilebx.hiVect();
    
                state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
                            BL_TO_FORTRAN_3D((*nChargePerCell_m[lev])[mfi]),
                            &tagval, &clearval, 
                            ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
                            ZFILL(dx), ZFILL(prob_lo), &time, &nCharge_m);
                //
                // Now update the tags in the TagBox.
                //
                tagfab.tags_and_untags(itags, tilebx);
            }
        }
    }

    void RemakeLevel (int lev, Real time,
                      const BoxArray& new_grids, const DistributionMapping& new_dmap)
    {
        SetBoxArray(lev, new_grids);
        SetDistributionMap(lev, new_dmap);
        
        bunch_mp->SetParticleBoxArray(lev, new_grids);
        bunch_mp->SetParticleDistributionMap(lev, new_dmap);
        
        nChargePerCell_m[lev].reset(new MultiFab(new_grids, new_dmap, 1, 1));
    }

    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap)
    {
        SetBoxArray(lev, new_grids);
        SetDistributionMap(lev, new_dmap);
        
        bunch_mp->SetParticleBoxArray(lev, new_grids);
        bunch_mp->SetParticleDistributionMap(lev, new_dmap);
        
        nChargePerCell_m[lev] = std::unique_ptr<MultiFab>(new MultiFab(new_grids, new_dmap, 1, 1));
    }
    
    void ClearLevel(int lev) {
        
        nChargePerCell_m[lev].reset(nullptr);
        
        ClearBoxArray(lev);
        ClearDistributionMap(lev);
    }
    
private:
    MyParticleContainer* bunch_mp;
    
    Array<std::unique_ptr<MultiFab> > nChargePerCell_m;
    
    double nCharge_m;
};


struct TestParams {
    int nx;
    int ny;
    int nz;
    int max_grid_size;
    int nppc;
    int nparticles;
    double sphere_radius;
    int nlevs;
    bool verbose;
};

void doTest(TestParams& parms)
{
    int nlevs = parms.nlevs;
    double halflength = parms.sphere_radius * 1.05;
    
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, -halflength);
        real_box.setHi(n, halflength);
    }
    
    IntVect domain_lo(0 , 0, 0); 
    IntVect domain_hi(parms.nx - 1, parms.ny - 1, parms.nz-1); 
    const Box domain(domain_lo, domain_hi);
    
    // use Cartesian coordinates
    int coord = 0;

    // Dirichlet boundary conditions
    Array<int> is_per = { 0, 0, 0};
    
    Geometry geom;
    geom.define(domain, &real_box, coord, &is_per[0]);
    
    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);
    DistributionMapping dm(ba);
    
    
    
    // init bunch
    MyParticleContainer myPC(geom, dm, ba);
    myPC.SetVerbose(false);
    
    myPC.initSphere(parms.nparticles, parms.sphere_radius);
    
    // init Amr object
    
    ParmParse pp("amr");
    pp.add("max_grid_size", parms.max_grid_size);
    
    Array<int> error_buf(nlevs, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    ParmParse pgeom("geometry");
    pgeom.addarr("is_periodic", is_per);
    
    Array<int> n_cell_in = { parms.nx, parms.ny, parms.nz  };
    
    
    MyAmr myAmr(&real_box, nlevs, n_cell_in, 0, &myPC);
    
    const Array<Geometry>& gv = myAmr.Geom();
    const Array<BoxArray>& bmv = myAmr.boxArray();
    const Array<DistributionMapping>& dmv = myAmr.DistributionMap();
    Array<int> rv;
    
    for (int i = 0; i < myAmr.maxLevel(); ++i)
        rv.push_back( myAmr.MaxRefRatio(i) );
    
    myPC.Define(gv, dmv, bmv, rv);
    
    myPC.Redistribute();
    
    for (int i = 0; i <= myAmr.finestLevel() && i < myAmr.maxLevel(); ++i)
        myAmr.regrid(i, 0.0);
    
    
    Array<std::unique_ptr<MultiFab> > partMF(nlevs);
//     Array<std::unique_ptr<MultiFab> > density(nlevs);
    
    for (int lev = 0; lev < nlevs; lev++) {
        partMF[lev].reset(new MultiFab(myAmr.boxArray(lev),
                                       myAmr.DistributionMap(lev),
                                       1, 2));
//         density[lev].reset(new MultiFab(bmv[lev], dmv[lev], 1, 0));

        partMF[lev]->setVal(0.0, 2);
//         density[lev]->setVal(0.0);

    }

    myPC.AssignDensityFort(0, partMF, 0, 1, nlevs - 1);
    
    for (int i = 0; i < nlevs; ++i) {
        
        std::cout << i << " " << partMF[i]->min(0) << " " << partMF[i]->max(0) << std::endl;
        
        if (partMF[i]->contains_nan()) {
            if ( ParallelDescriptor::IOProcessor() )
                std::cout << "Grid contains NAN(s)" << std::endl;
        } else if (partMF[i]->contains_inf()) {
            if ( ParallelDescriptor::IOProcessor() )
                std::cout << "Grid contains INF(s)" << std::endl;
        }
    }
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("nppc", parms.nppc);
  pp.get("max_grid_size", parms.max_grid_size);
  pp.get("nparticles", parms.nparticles);
  pp.get("sphere_radius", parms.sphere_radius);
  
  parms.verbose = false;
  pp.query("verbose", parms.verbose);
  
  pp.get("nlevs", parms.nlevs);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles : ";
    std::cout << parms.nparticles  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }
  
  doTest(parms);
  
  amrex::Finalize();
}
