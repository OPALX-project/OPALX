/*!
 * @file testSolverFail.cpp
 * @details Solves Poisson equation
 * @authors Matthias Frey
 * @date July 2017
 * @brief
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

#include "../Solver.h"

using namespace amrex;


/* NStructReal = 4 (mass, position are already reserved)
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
    : public ParticleContainer<BL_SPACEDIM + 1 /* momenta */>
{
 public:

    MyParticleContainer(const Geometry            & geom, 
                        const DistributionMapping & dmap,
                        const BoxArray            & ba)
        : ParticleContainer<BL_SPACEDIM + 1> (geom, dmap, ba), particles_rm(GetParticles())
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
        
        std::cout << nParticles << " " << nLocParticles << std::endl;
    
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
            
            p.m_idata.id  = ParticleType::NextID();
            p.m_idata.cpu = MyProc;
                
            particles_rm[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
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

	this->MakeNewLevel(0, 0.0, ba, dm);
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

	Array<std::unique_ptr<MultiFab> > nChargePerCell_m(max_level + 1);

        for (int i = 0; i <= finest_level; ++i) {
	    nChargePerCell_m[i].reset(new MultiFab(grids[i], dmap[i], 1, 2));
            nChargePerCell_m[i]->setVal(0.0, 2);
        }

        bunch_mp->AssignDensityFort(0, nChargePerCell_m, 0, 1, finest_level);
        
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
    }

    void MakeNewLevel (int lev, Real time,
                       const BoxArray& new_grids, const DistributionMapping& new_dmap)
    {
        SetBoxArray(lev, new_grids);
        SetDistributionMap(lev, new_dmap);
        
        bunch_mp->SetParticleBoxArray(lev, new_grids);
        bunch_mp->SetParticleDistributionMap(lev, new_dmap);
    }
    
    void ClearLevel(int lev) {
        ClearBoxArray(lev);
        ClearDistributionMap(lev);
    }
    
private:
    MyParticleContainer* bunch_mp;
    
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

void solve(MyAmr& myAmr, MyParticleContainer& myPC) {
    int maxiter = 100;
    int maxiter_b = 100;
    int verbose = 0;
    bool usecg = true;
    double bottom_solver_eps = 1.0e-4;
    int max_nlevel = 1024;
    
    /* MG_SMOOTHER_GS_RB  = 1
     * MG_SMOOTHER_JACOBI = 2
     * MG_SMOOTHER_MINION_CROSS = 5
     * MG_SMOOTHER_MINION_FULL = 6
     * MG_SMOOTHER_EFF_RB = 7
     */
    int smoother = 1;
    
    // #smoothings at each level on the way DOWN the V-cycle
    int nu_1 = 2;
    
    // #smoothings at each level on the way UP the V-cycle
    int nu_2 = 2;
    
    // #smoothings before and after the bottom solver
    int nu_b = 0;
    
    // #smoothings
    int nu_f = 8;

    /* MG_FCycle = 1  (full multigrid)
     * MG_WCycle = 2
     * MG_VCycle = 3
     * MG_FVCycle = 4
     */
    int cycle = 1;
    
    bool cg_solver = true;
    
    /* if cg_solver == true:
     * - BiCG --> 1
     * - CG --> 2
     * - CABiCG --> 3
     * 
     * else if cg_solver == false
     * - CABiCG is taken
     */
    int bottom_solver = 1;
    
    ParmParse pp("mg");

    pp.add("maxiter", maxiter);
    pp.add("maxiter_b", maxiter_b);
    pp.add("nu_1", nu_1);
    pp.add("nu_2", nu_2);
    pp.add("nu_b", nu_b);
    pp.add("nu_f", nu_f);
    pp.add("v"   , verbose);
    pp.add("usecg", usecg);
    pp.add("cg_solver", cg_solver);

    pp.add("rtol_b", bottom_solver_eps);
    pp.add("numLevelsMAX", max_nlevel);
    pp.add("smoother", smoother);
    pp.add("cycle_type", cycle); // 1 -> F, 2 -> W, 3 -> V, 4 -> F+V
    //
    // The C++ code usually sets CG solver type using cg.cg_solver.
    // We'll allow people to also use mg.cg_solver but pick up the former as well.
    //
    if (!pp.query("cg_solver", cg_solver))
    {
        ParmParse pp("cg");

        pp.add("cg_solver", cg_solver);
    }

    pp.add("bottom_solver", bottom_solver);
    
    
    /*
     * setup containers
     */
    int base_level = 0;
    int finest_level = myAmr.finestLevel();
    int nlevs = finest_level + 1;
    
    Array<std::unique_ptr<MultiFab> > partMF(nlevs);
    Array<std::unique_ptr<MultiFab> > rho(nlevs), phi(nlevs), efield(nlevs);
    
    for (int i = 0; i < nlevs; ++i) {
        const DistributionMapping& dm = myAmr.DistributionMap(i);
        const BoxArray& ba = myAmr.boxArray(i);
        rho[i] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm,    1,           0));
        phi[i] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm,    1,           1));
        efield[i] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm, BL_SPACEDIM, 1));
    
        rho[i]->setVal(0.0);
        phi[i]->setVal(0.0, 1);
        efield[i]->setVal(0.0, 1);
    }
    

    for (int lev = 0; lev < nlevs; lev++) {
        partMF[lev].reset(new MultiFab(myAmr.boxArray(lev),
                                       myAmr.DistributionMap(lev),
                                       1, 2));
        partMF[lev]->setVal(0.0, 2);
    }

    myPC.AssignDensityFort(0, partMF, 0, 1, nlevs - 1);
    
    for (int i = 0; i < nlevs; ++i)
        MultiFab::Copy(*rho[i], *partMF[i], 0, 0, 1, 0);
    
    
    double l0norm = rho[finest_level]->norm0(0);
    for (int i = 0; i <= finest_level; ++i) {
        rho[i]->mult(1.0 / l0norm, 0, 1);
    }
    
    const Array<Geometry>& geom = myAmr.Geom();
    
    Solver sol;
    
    sol.solve_for_accel(rho,            // [V m]
                        phi,            // [V m^3]
                        efield,       // [V m^2]
                        geom,
                        base_level,
                        finest_level,
                        0.0);
    
}

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
    
    std::cout << "Total num: " << myPC.TotalNumberOfParticles() << std::endl;
    
    for (int i = 0; i <= myAmr.finestLevel() && i < myAmr.maxLevel(); ++i)
        myAmr.regrid(i, 0.0);
    
    solve(myAmr, myPC);
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
