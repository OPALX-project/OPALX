#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>

#include <boost/filesystem.hpp>


#include <ParmParse.H>

#include "PartBunch.h"
#include "AmrPartBunch.h"

#include "Distribution.h"
#include "Solver.h"
#include "AmrOpal.h"

#include "helper_functions.h"

#include "writePlotFile.H"

#include <PlotFileUtil.H>

#include <cmath>

#include "Physics/Physics.h"

#include "Particle/BoxParticleCachingPolicy.h"
#define Dim 3
typedef UniformCartesian<Dim, double>                                 Mesh_t;
typedef BoxParticleCachingPolicy<double, Dim, Mesh_t>                 CachingPolicy_t;
typedef ParticleSpatialLayout<double, Dim, Mesh_t, CachingPolicy_t>   playout_t;
typedef Cell                                                          Center_t;
typedef CenteredFieldLayout<Dim, Mesh_t, Center_t>                    FieldLayout_t;
typedef Field<double, Dim, Mesh_t, Center_t>                          Field_t;

typedef IntCIC                                                        IntrplCIC_t;

typedef UniformCartesian<2, double>                                   Mesh2d_t;
typedef CenteredFieldLayout<2, Mesh2d_t, Center_t>                    FieldLayout2d_t;
typedef Field<double, 2, Mesh2d_t, Center_t>                          Field2d_t;

template<class PL>
class TempParticle : public IpplParticleBase<PL> {
    
public:
    ParticleAttrib<double>     	Q;
    ParticleAttrib<Vektor<double,2> > Rphase; //velocity of the particles
    ParticleAttrib<Vector_t>	v; //velocity of the particles
    
    TempParticle(PL* pl,
                 Vektor<double,3> nr,
                 e_dim_tag decomp[Dim],
                 Vektor<double,3> extend_l_,
                 Vektor<double,3> extend_r_,
                 Vektor<int,3> Nx_,
                 Vektor<int,3> Nv_,
                 Vektor<double,3> Vmax_) :
        IpplParticleBase<PL>(pl),
        nr_m(nr),
        extend_l(extend_l_),
        extend_r(extend_r_),
        Nx(Nx_), Nv(Nv_), Vmax(Vmax_)
    {
        this->addAttribute(Q);
        this->addAttribute(Rphase);
        this->addAttribute(v);
    }
    
    void interpolate_distribution(Vektor<double,3> dx, Vektor<double,3> dv, int step, std::string dir){
        
        double spacings[2] = {(extend_r[2]-extend_l[2])/(Nx[2]),2.*Vmax[2]/(Nv[2])};
        Vektor<double,2> origin;
        origin(0) = extend_l[2]; origin(1) = -Vmax[2]; 
        Index I(Nx[2]+1); Index J(Nv[2]+1);
        
        NDIndex<2> domain2d_m;
        domain2d_m[0]=I; domain2d_m[1]=J;

        Mesh2d_t mesh2d_m = Mesh2d_t(domain2d_m, spacings, origin);
        std::unique_ptr<FieldLayout2d_t>  layout2d_m = std::unique_ptr<FieldLayout2d_t>( new FieldLayout2d_t(mesh2d_m) );
        
        //set origin and spacing is needed for correct results, even if mesh was created with these paremeters ?!
        mesh2d_m.set_meshSpacing(&(spacings[0]));
        mesh2d_m.set_origin(origin);

        domain2d_m = layout2d_m->getDomain();
        //f_m is used for twostream instability as 2D phase space mesh
        
        BConds<double, 2, Mesh2d_t, Center_t> BC;
        if (Ippl::getNodes()>1) {
            BC[0] = new ParallelInterpolationFace<double,2,Mesh2d_t,Cell>(0);
            BC[1] = new ParallelInterpolationFace<double,2,Mesh2d_t,Cell>(1);
            BC[2] = new ParallelInterpolationFace<double,2,Mesh2d_t,Cell>(2);
            BC[3] = new ParallelInterpolationFace<double,2,Mesh2d_t,Cell>(3);
        }
        else {
            BC[0] = new InterpolationFace<double,2,Mesh2d_t,Cell>(0);
            BC[1] = new InterpolationFace<double,2,Mesh2d_t,Cell>(1);
            BC[2] = new InterpolationFace<double,2,Mesh2d_t,Cell>(2);
            BC[3] = new InterpolationFace<double,2,Mesh2d_t,Cell>(3);
        }
        
        Field2d_t f_m;
        f_m.initialize(mesh2d_m, *(layout2d_m.get()), GuardCellSizes<2>(1),BC);
        for (unsigned i=0; i<this->getLocalNum(); ++i)
            Rphase[i]=Vektor<double,2>(this->R[i][2],v[i][2]);
        this->Q.scatter(f_m, this->Rphase, IntrplCIC_t());
        
        std::ofstream csvout;
        csvout.precision(10);
        csvout.setf(std::ios::scientific, std::ios::floatfield);
        
        std::cout << "Hi 6b" << std::endl;
        
        std::stringstream fname;
        fname << dir << "/f_mesh_";
        fname << std::setw(4) << std::setfill('0') << step;
        fname << ".csv";
        
        std::cout << "Hi 6c" << std::endl;
        
        // open a new data file for this iteration
        // and start with header
        csvout.open(fname.str().c_str(), std::ios::out);
        csvout << "z, vz, f" << std::endl;
        NDIndex<2> lDom = domain2d_m;
        
        std::cout << "Hi 6d" << std::endl;
        
        std::cout << lDom[0].first() << " " << lDom[0].last() << std::endl
                  << lDom[1].first() << " " << lDom[1].last() << std::endl;
        
        for (int i=lDom[0].first(); i<=lDom[0].last(); i++) {
        
            for (int j=lDom[1].first(); j<=lDom[1].last(); j++) {
            
                csvout << (i+0.5)*dx[2] << ","
                       << (j+0.5)*dv[2]-Vmax[2]
                       << "," << f_m[i][j].get() << std::endl;
            }
        }
        // close the output file for this iteration:
        csvout.close();
        
//         fname.str("");
        
//         delete BC[0];
    }
    
public:
    Vektor<int,Dim> nr_m;
    Vektor<double,Dim> extend_l;
    Vektor<double,Dim> extend_r;
    Vektor<int,Dim> Nx;
    Vektor<int,Dim> Nv; 
    Vektor<double,Dim> Vmax;
};

void writeEnergy(PartBunchBase* bunch,
                 container_t& rho,
                 container_t& phi,
                 container_t& efield,
                 const Array<int>& rr,
                 const Array<Geometry>& geom,
                 int step,
                 std::string dir = "./")
{
    for (int lev = efield.size() - 2; lev >= 0; lev--)
        BoxLib::average_down(efield[lev+1], efield[lev], 0, 3, rr[lev]);
    
    double cell_volume_0 = geom[0].CellSize(0) *
                           geom[0].CellSize(1) *
                           geom[0].CellSize(2);
    
    // field energy (Ulmer version, i.e. cell_volume instead #points)
    double field_energy = 0.5 * cell_volume_0 * MultiFab::Dot(efield[0], 0, efield[0], 0, 3, 0);
    
    // kinetic energy
    double ekin_local = 0.0;
    double ekin_global = 0.0;
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        Vector_t PP = bunch->getP(i);
        ekin_local += 0.5 * (
            PP(0) * PP(0) +
            PP(1) * PP(1) +
            PP(2) * PP(2)
        );
    }
    
    reduce(ekin_local, ekin_global, OpAddAssign());
    
    // potential energy
    double integral_phi = 0.0;
    for (int lev = 0; lev < rho.size(); ++lev) {
        double cell_volume = geom[lev].CellSize(0) *
                             geom[lev].CellSize(1) *
                             geom[lev].CellSize(2);
        rho[lev].mult(cell_volume, 0, 1);
        MultiFab::Multiply(phi[lev], rho[lev], 0, 0, 1, 0);
        integral_phi += 0.5 * phi[lev].sum(0);
    }
    
    if(Ippl::myNode()==0) {
        std::ofstream csvout;
        csvout.precision(10);
        csvout.setf(std::ios::scientific, std::ios::floatfield);

        std::stringstream fname;
        fname << dir << "/energy";
        fname << ".csv";

        // open a new data file for this iteration
        // and start with header
        csvout.open(fname.str().c_str(), std::ios::out | std::ofstream::app);
        
        if (step == 0) {
            csvout << "it,Efield,Ekin,Etot,Epot" << std::endl;
        }
        
        csvout << step << ", "
               << field_energy << ","
               << ekin_global << ","
               << field_energy + ekin_global << "," 
               << integral_phi << std::endl;
        
        csvout.close();
    }
}

void ipplProjection(std::unique_ptr<TempParticle<playout_t> >& P,
                    PartBunchBase* bunch, int step, std::string dir = "./") {
    
    //**************
    // copy the particles to Ippl container
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        P->Q[i] = bunch->getQM(i);
        P->v[i] = bunch->getP(i);
        P->R[i] = bunch->getR(i);
    }
    //**************
    std::cout << "Hi 5" << std::endl;
    
    P->interpolate_distribution((P->extend_r-P->extend_l)/(P->Nx),2.*P->Vmax/(P->Nv), step, dir);
}
    
    

void printLongitudinalPhaseSpace(PartBunchBase* bunch, const AmrOpal& myAmrOpal, int step) {
    
//     // charge projection onto the (z, vz) phase space
//     /*
//      * We assign the charge, z, vz onto the grid and then sum up the x and y directions of the grid
//      */
//     container_t data, z, vz;
//     data.resize(myAmrOpal.finestLevel() + 1);
//     
//     for (int lev = 0; lev < myAmrOpal.finestLevel() + 1; ++lev) {
// #ifdef UNIQUE_PTR
//         //                                                  # component # ghost cells                                                                                                                                          
//         data[lev] = std::unique_ptr<MultiFab>(new MultiFab(myAmrOpal.boxArray()[lev], 8, 0));
//         data[lev]->setVal(0.0);
// #else
//         //                       # component # ghost cells                                                                                                                                          
//         data.set(lev, new MultiFab(myAmrOpal.boxArray()[lev], 8, 0));
//         data[lev].setVal(0.0);
// #endif
//     }
//     
//     dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, data, 0, 8, myAmrOpal.finestLevel() + 1);
    
    
    
    
    
    // 03. February 2017,
    // http://stackoverflow.com/questions/225362/convert-a-number-to-a-string-with-specified-length-in-c
    std::stringstream num;
    num << "density"  << std::setw(4) << std::setfill('0') << step << ".dat";
    std::string plotfilename = num.str();
    
    std::ofstream out(plotfilename);
    
    
//     struct triple_t {
//         int i, j;
//         double f;
//     };
    
    out << "z vz f" << std::endl;
    for (std::size_t j = 0; j < bunch->getLocalNum(); ++j)
        out << bunch->getR(j)(2) << ", " << bunch->getP(j)(2) << ", " << bunch->getQM(j) << std::endl;
    
    
    
//     for (int lev = 0; lev < myAmrOpal.finestLevel() + 1; ++lev) {
//         
//         for (MFIter mfi(data[lev], false); mfi.isValid(); ++mfi) {
//         
//             const Box&  bx  = mfi.validbox();
//             FArrayBox& fab = data[lev][mfi];
//             
//             std::cout << bx.loVect()[0] << " " << bx.hiVect()[0] << std::endl;
//             std::cout << bx.loVect()[1] << " " << bx.hiVect()[1] << std::endl;
//             std::cout << bx.loVect()[2] << " " << bx.hiVect()[2] << std::endl;
//             
//             int iext = bx.hiVect()[0] - bx.loVect()[0] + 1;
//             int jext = bx.hiVect()[1] - bx.loVect()[1] + 1;
//             std::vector<triple_t> bla(iext * jext);
//             
//             for (int i = bx.loVect()[0]; i <= bx.hiVect()[0]; ++i) {
//                 for (int j = bx.loVect()[1]; j <= bx.hiVect()[1]; ++j) {
//                     double z = 0.0, vz = 0.0, f = 0.0;
//                     bla[j + i * jext].i = i;
//                     bla[j + i * jext].j = j;
//                     
//                     for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; ++k ) {
//                     
//                         IntVect iv(i, j, k);
//                         
// //                         z  += fab(iv, 7);
// //                         vz += fab(iv, 4);
// //                         f  += fab(iv, 0);
//                         std::cout << i << " " << j << " " << k << " " << fab(iv, 0) << std::endl; std::cin.get(); 
//                         bla[j + i * jext].f += fab(iv, 0);
//                         
// //                         for (int ii = 0; ii < 8; ++ii)
// //                             std::cout << fab(iv, ii) << " ";
// //                         std::cout << std::endl; //std::cin.get();
//                     }
//                 }
// //                 out << i << ", " << j << ", " << f << std::endl;
//             }
//             for (std::vector<triple_t>::size_type i = 0; i < bla.size(); ++i) {
//                 out << bla[i].i << ", " << bla[i].j << ", " << bla[i].f << std::endl;
//             }
//         }
//     }
    
    out.close();
}

void doSolve(AmrOpal& myAmrOpal, PartBunchBase* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels,
             Inform& msg,
             std::string dir = "./")
{
    // =======================================================================                                                                                                                                   
    // 4. prepare for multi-level solve                                                                                                                                                                          
    // =======================================================================
    
    rhs.resize(nLevels);
    phi.resize(nLevels);
    grad_phi.resize(nLevels);
    
    for (int lev = 0; lev < nLevels; ++lev) {
        initGridData(rhs, phi, grad_phi, myAmrOpal.boxArray()[lev], lev);
    }

    // Define the density on level 0 from all particles at all levels                                                                                                                                            
    int base_level   = 0;
    int finest_level = myAmrOpal.finestLevel();

    dynamic_cast<AmrPartBunch*>(bunch)->AssignDensity(0, false, rhs, base_level, 1, finest_level);
    
    // eps in C / (V * m)
    double constant = -1.0;
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(constant, 0, 1);
#endif
    }
    
    // **************************************************************************                                                                                                                                
    // Compute the total charge of all particles in order to compute the offset                                                                                                                                  
    //     to make the Poisson equations solvable                                                                                                                                                                
    // **************************************************************************                                                                                                                                

    Real offset = 0.;
    
    if ( geom[0].isAllPeriodic() ) {
        double sum = 0, max = 0;
        for (int lev = 0; lev < rhs.size(); ++lev) {
            sum = rhs[lev].sum(0);
            max = rhs[lev].max(0);
        }
        msg << "total charge in density field before ion subtraction is " << sum << endl;
        msg << "max total charge in density field before ion subtraction is " << max << endl;
        for (std::size_t i = 0; i < bunch->getLocalNum(); ++i)
                offset += bunch->getQM(i);
        offset /= geom[0].ProbSize();
//         offset = -1.0;
    }


    // solve                                                                                                                                                                                                     
    Solver sol;
    sol.solve_for_accel(rhs,
                        phi,
                        grad_phi,
                        geom,
                        base_level,
                        finest_level,
                        offset,
                        false);
    
    // for plotting unnormalize
    for (int i = 0; i <=finest_level; ++i) {
#ifdef UNIQUE_PTR
        rhs[i]->mult(1.0 / constant, 0, 1);       // in [V m]
#else
        rhs[i].mult(1.0 / constant, 0, 1);
#endif
    }
}

void doPlasma(Vektor<std::size_t, 3> nr,
              std::size_t nLevels,
              std::size_t maxBoxSize,
              double dt,
              std::size_t nIter,
              Distribution::Type type,
              Inform& msg)
{
    std::array<double, BL_SPACEDIM> lower = {{0.0, 0.0, 0.0}}; // m
    std::array<double, BL_SPACEDIM> upper = {{4.0 * Physics::pi,
                                              4.0 * Physics::pi,
                                              4.0 * Physics::pi}
    }; // m
    
    RealBox domain;
    
    // in helper_functions.h
    init(domain, nr, lower, upper);
    
    
    /*
     * create an Amr object
     */
    ParmParse pp("amr");
    pp.add("max_grid_size", int(maxBoxSize));
    
    Array<int> error_buf(nLevels, 0);
    
    pp.addarr("n_error_buf", error_buf);
    pp.add("grid_eff", 0.95);
    
    ParmParse pgeom("geometry");
    Array<int> is_per = { 1, 1, 1};
    pgeom.addarr("is_periodic", is_per);
    
    Array<int> nCells(3);
    for (int i = 0; i < 3; ++i)
        nCells[i] = nr[i];
    
    AmrOpal myAmrOpal(&domain, nLevels - 1, nCells, 0 /* cartesian */);
    
    
    // ========================================================================
    // 2. initialize all particles (just single-level)
    // ========================================================================
    
    const Array<BoxArray>& ba = myAmrOpal.boxArray();
    const Array<DistributionMapping>& dmap = myAmrOpal.DistributionMap();
    const Array<Geometry>& geom = myAmrOpal.Geom();
    
    Array<int> rr(nLevels);
    for (std::size_t i = 0; i < nLevels; ++i)
        rr[i] = 2;
    
    PartBunchBase* bunch = new AmrPartBunch(geom, dmap, ba, rr);
    
    
    // initialize a particle distribution
    Distribution dist;
    
    Vektor<double, 3> extend_l = Vektor<double, 3>(lower[0], lower[1], lower[2]);
    Vektor<double, 3> extend_r = Vektor<double, 3>(upper[0], upper[1], upper[2]);
    
    Vektor<std::size_t, 3> Nx, Nv;
    Vektor<double, 3> Vmax;
    
    std::string dirname = "";
    
    if ( type == Distribution::Type::kTwoStream ) {
        dirname = "twostream";
        Nx = Vektor<std::size_t, 3>(4, 4, 32);
        Nv = Vektor<std::size_t, 3>(8, 8, 128);
        Vmax = Vektor<double, 3>(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     type,
                     0.05);
    } else if ( type == Distribution::Type::kRecurrence ) {
        dirname = "recurrence";
        Nx = Vektor<std::size_t, 3>(8, 8, 8);
        Nv = Vektor<std::size_t, 3>(32, 32, 32);
        Vmax = Vektor<double, 3>(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     type,
                     0.01);
    } else if ( type == Distribution::Type::kLandauDamping ) {
        dirname = "landau";
        Nx = Vektor<std::size_t, 3>(8, 8, 8);
        Nv = Vektor<std::size_t, 3>(32, 32, 32);
        Vmax = Vektor<double, 3>(6.0, 6.0, 6.0);
        
        dist.special(extend_l,
                     extend_r,
                     Nx,
                     Nv,
                     Vmax,
                     type,
                     0.05);
    }
    
    dirname += (
        "-data-grid-" +
        std::to_string(nr[0]) + "-" + 
        std::to_string(nr[1]) + "-" + 
        std::to_string(nr[2])
    );
    
    boost::filesystem::path dir(dirname);
    if ( Ippl::myNode() == 0 )
        boost::filesystem::create_directory(dir);
    
    // copy particles to the PartBunchBase object.
    dist.injectBeam(*bunch);
    
    // redistribute on single-level
    bunch->myUpdate();
    
//     dynamic_cast<AmrPartBunch*>(bunch)->SetAllowParticlesNearBoundary(true);
    
//     writePlotFile(plotsolve, rhs, phi, grad_phi, rr, geoms, 0);
    
    myAmrOpal.setBunch(dynamic_cast<AmrPartBunch*>(bunch));
    
    const Array<Geometry>& geoms = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    // -----------------------------------------------------
    
    e_dim_tag decomp[Dim];
    Mesh_t *mesh;
    FieldLayout_t *FL;

    NDIndex<Dim> domain_ippl;
    for (unsigned i=0; i<Dim; i++)
        domain_ippl[i] = domain_ippl[i] = Index(nr[i]+1);

    for (unsigned d=0; d < Dim; ++d)
        decomp[d] = PARALLEL;
    
    std::cout << "Hi 2" << std::endl;
    
    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain_ippl);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);
    
    std::unique_ptr<TempParticle<playout_t> > P( new TempParticle<playout_t>(PL, nr,
                                                                             decomp,
                                                                             Vektor<double, 3>(lower[0], lower[1], lower[2]),
                                                                             Vektor<double, 3>(upper[0], upper[1], upper[2]),
                                                                             Vektor<std::size_t, 3>(4, 4, 32),
                                                                             Vektor<std::size_t, 3>(8, 8, 128),
                                                                             Vektor<double, 3>(6.0, 6.0, 6.0)) );
    P->create(bunch->getLocalNum());
    // ----------------------------------------------------
    
    
    container_t rhs(PArrayManage);
    container_t phi(PArrayManage);
    container_t grad_phi(PArrayManage);
    doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geoms, rr, nLevels, msg, dir.string());
    
    writeEnergy(bunch, rhs, phi, grad_phi, rr, geom, 0, dir.string());
    
    for (std::size_t i = 0; i < nIter; ++i) {
        msg << "Processing step: " << i << endl;
        
        if ( type == Distribution::Type::kTwoStream ) {
            ipplProjection(P, bunch, i, dir.string());
            msg << "Done writing projection." << endl;
        }
        
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            int level = dynamic_cast<AmrPartBunch*>(bunch)->getLevel(j);
            bunch->setR( bunch->getR(j) + dt * bunch->getP(j), j );
        }
        
        msg << "Done pushing particles." << endl;
            
        if ( myAmrOpal.maxLevel() > 0 )
            for (int j = 0; j <= myAmrOpal.finestLevel() && j < myAmrOpal.maxLevel(); ++j)
                myAmrOpal.regrid(j /*lbase*/, 0.0 /*time*/);
        else
            bunch->myUpdate();
        
        msg << "Done updating." << endl;
        
        doSolve(myAmrOpal, bunch, rhs, phi, grad_phi, geoms, rr, nLevels, msg, dir.string());
        
        msg << "Done solving Poisson's equation." << endl;
        
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            int level = dynamic_cast<AmrPartBunch*>(bunch)->getLevel(j);
            Vector_t Ef = dynamic_cast<AmrPartBunch*>(bunch)->interpolate(j, grad_phi[level]);
            
            bunch->setP( bunch->getP(j) + dt * bunch->getQM(j) / bunch->getMass(j) * Ef, j);
        }
        
        writeEnergy(bunch, rhs, phi, grad_phi, rr, geom, i + 1, dir.string());
        
        msg << "Done with step: " << i << endl;
    }
    
    // ----------------------------------------------
    P->destroy(bunch->getLocalNum(), 0);
    P.reset(nullptr);
    delete FL;
    delete mesh;
    
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc,argv, false);
    
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels]"
         << " [max. box size] [timstep] [#iterations]" "[test: twostream, recurrence, landau]"
         << " [out: timing file name (optiona)]";
    
    if ( argc < 9 ) {
        std::cerr << call.str() << std::endl;
        return -1;
    }
    
    // number of grid points in each direction
    Vektor<std::size_t, 3> nr(std::atoi(argv[1]),
                              std::atoi(argv[2]),
                              std::atoi(argv[3]));
    
    
    std::size_t nLevels    = std::atoi( argv[4] ) + 1; // i.e. nLevels = 0 --> only single level
    std::size_t maxBoxSize = std::atoi( argv[5] );
    double dt              = std::atof( argv[6] );
    std::size_t nIter      = std::atof( argv[7] );
    std::string test       = argv[8];
    
    Inform msg(argv[8]);
    
    Distribution::Type type = Distribution::Type::kTwoStream;
    if ( test == "twostream" )
        type = Distribution::Type::kTwoStream;
    else if ( test == "recurrence" )
        type = Distribution::Type::kRecurrence;
    else if ( test == "landau" )
        type = Distribution::Type::kLandauDamping;
    else {
        msg << "No such plasma test." << endl;
        return -1;
    }
    
    msg << "Particle test running with" << endl
        << "- #level     = " << nLevels << endl
        << "- grid       = " << nr << endl
        << "- max. size  = " << maxBoxSize << endl
        << "- time step  = " << dt << endl
        << "- #steps     = " << nIter << endl
        << "- test       = " << test << endl;
    
    doPlasma(nr, nLevels, maxBoxSize,
             dt, nIter, type, msg);
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    if ( argc == 10 ) {
        timefile.str(std::string());
        timefile << std::string(argv[8]);
    }
    
    IpplTimings::print(timefile.str());
    
    return 0;
}