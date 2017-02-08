#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <set>
#include <sstream>
#include <iomanip>


#include <ParmParse.H>

#include "../Distribution.h"
#include "../Solver.h"
#include "../AmrOpal.h"

#include "../helper_functions.h"

// #include "writePlotFile.H"

#include <cmath>

#include "Physics/Physics.h"

typedef AmrOpal::amrplayout_t amrplayout_t;
typedef AmrOpal::amrbase_t amrbase_t;
typedef AmrOpal::amrbunch_t amrbunch_t;

typedef Vektor<double, BL_SPACEDIM> Vector_t;

// ----------------------------------------------------------------------------

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
        
        double spacings[2] = {(extend_r[2]-extend_l[2])/(Nx[2]),2.*Vmax[2]/(Nv[2])};
        Vektor<double,2> origin;
        origin(0) = extend_l[2]; origin(1) = -Vmax[2]; 
        Index I(Nx[2]+1); Index J(Nv[2]+1);
        domain2d_m[0]=I; domain2d_m[1]=J;

        mesh2d_m=Mesh2d_t(domain2d_m, spacings, origin);
        layout2d_m = new FieldLayout2d_t(mesh2d_m);
        
        //set origin and spacing is needed for correct results, even if mesh was created with these paremeters ?!
        mesh2d_m.set_meshSpacing(&(spacings[0]));
        mesh2d_m.set_origin(origin);

        domain2d_m = layout2d_m->getDomain();
        //f_m is used for twostream instability as 2D phase space mesh
        
        BConds<double,2,UniformCartesian<2,double>,Cell> BC;
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
        
        f_m.initialize(mesh2d_m, *layout2d_m, GuardCellSizes<2>(1),BC);
    }
    
    void interpolate_distribution(Vektor<double,3> dx, Vektor<double,3> dv){
        f_m=0;
        for (unsigned i=0; i<this->getLocalNum(); ++i) {
//             SpaceQ[i]=Q[i]/(dx[0]*dx[1]*dv[0]*dv[1]);
            Rphase[i]=Vektor<double,2>(this->R[i][2],v[i][2]);
        }
        //this->SpaceQ.scatter(this->f_m, this->Rphase, IntrplCIC_t());
        this->Q.scatter(this->f_m, this->Rphase, IntrplCIC_t());
    }
    
public:
    Vektor<int,Dim> nr_m;
    Vektor<double,Dim> extend_l;
    Vektor<double,Dim> extend_r;
    Vektor<int,Dim> Nx;
    Vektor<int,Dim> Nv; 
    Vektor<double,Dim> Vmax;
    //Fields for tracking distribution function
    Field2d_t f_m;
    Mesh2d_t mesh2d_m;
    NDIndex<2> domain2d_m;	
    FieldLayout2d_t *layout2d_m;
    
};

void ipplProjection(Vektor<int,Dim> nr,
                    Vektor<double,Dim> extend_l,
                    Vektor<double,Dim> extend_r,
                    Vektor<int,Dim> Nx,
                    Vektor<int,Dim> Nv,
                    Vektor<double,Dim> Vmax,
                    amrbunch_t* bunch, int step) {
    e_dim_tag decomp[Dim];
    Mesh_t *mesh;
    FieldLayout_t *FL;

    NDIndex<Dim> domain;
    for (unsigned i=0; i<Dim; i++)
        domain[i] = domain[i] = Index(nr[i]+1);

    for (unsigned d=0; d < Dim; ++d)
        decomp[d] = PARALLEL;
    
    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);
    
    TempParticle<playout_t>  *P = new TempParticle<playout_t>(PL, nr, decomp, extend_l, extend_r,Nx,Nv,Vmax);
    
    //**************
    // copy the particles to Ippl container
    for (std::size_t i = 0; i < bunch->getLocalNum(); ++i) {
        P->create(1);
        P->Q[i] = bunch->qm[i];
        P->v[i] = bunch->P[i];
        P->R[i] = bunch->R[i];
    }
    //**************
    
    
    P->interpolate_distribution((extend_r-extend_l)/(Nx),2.*Vmax/(Nv));
    
    // ---------------
    Vektor<double,3> dx = (P->extend_r-P->extend_l)/(P->Nx);
    
    Vektor<double,3> dv = 2.*P->Vmax/(P->Nv);
    std::ofstream csvout;
    csvout.precision(10);
    csvout.setf(std::ios::scientific, std::ios::floatfield);

    std::stringstream fname;
    fname << "data/f_mesh_";
    fname << std::setw(4) << std::setfill('0') << step;
    fname << ".csv";
    
    // open a new data file for this iteration
    // and start with header
    csvout.open(fname.str().c_str(), std::ios::out);
    csvout << "z, vz, f" << std::endl;
    NDIndex<2> lDom = P->domain2d_m;
    
    for (int i=lDom[0].first(); i<=lDom[0].last(); i++) {
    
        for (int j=lDom[1].first(); j<=lDom[1].last(); j++) {
        
            csvout << (i+0.5)*dx[2] << ","
                   << (j+0.5)*dv[2]-P->Vmax[2]
                   << "," << P->f_m[i][j].get() << std::endl;
        }
    }
    // close the output file for this iteration:
    csvout.close();
    // -------------------
    
    delete P;
//     delete PL;
    delete FL;
    delete mesh;
}

// ------------------------------------------------------------------------------------------

void doSolve(AmrOpal& myAmrOpal, amrbunch_t* bunch,
             container_t& rhs,
             container_t& phi,
             container_t& grad_phi,
             const Array<Geometry>& geom,
             const Array<int>& rr,
             int nLevels)
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
    
    bunch->AssignDensity(bunch->qm, false, rhs, base_level, finest_level);
    
    // eps in C / (V * m)
    double constant = -1.0 / Physics::epsilon_0;  // in [V m / C]
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

void doTwoStream(Vektor<std::size_t, 3> nr,
                 std::size_t nLevels,
                 std::size_t maxBoxSize,
                 double dt,
                 std::size_t nIter,
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
    
    
    amrplayout_t* playout = new amrplayout_t(geom, dmap, ba, rr);
    
    std::unique_ptr<amrbunch_t> bunch( new amrbunch_t() );
    bunch->initialize(playout);
    bunch->initializeAmr(); // add attributes: level, grid
    
    
    // initialize a particle distribution
    Distribution dist;
    dist.twostream(Vektor<double, 3>(lower[0], lower[1], lower[2]),
                   Vektor<double, 3>(upper[0], upper[1], upper[2]),
                   Vektor<std::size_t, 3>(4, 4, 32),
                   Vektor<std::size_t, 3>(8, 8, 128),
                   Vektor<double, 3>(6.0, 6.0, 6.0), 0.05);
    
    // copy particles to the PartBunchAmr object.
    dist.injectBeam( *(bunch.get()) );
    
    // redistribute on single-level
    bunch->update();
    
    myAmrOpal.setBunch( bunch.get() );
    
    const Array<Geometry>& geoms = myAmrOpal.Geom();
    
    for (int i = 0; i <= myAmrOpal.finestLevel() && i < myAmrOpal.maxLevel(); ++i)
        myAmrOpal.regrid(i /*lbase*/, 0.0 /*time*/);
    
    msg << "Multi-level statistics" << endl;
    bunch->gatherStatistics();
    
    for (std::size_t i = 0; i < nIter; ++i) {
        
        ipplProjection(nr,
                       Vektor<double, 3>(lower[0], lower[1], lower[2]),
                       Vektor<double, 3>(upper[0], upper[1], upper[2]),
                       Vektor<std::size_t, 3>(4, 4, 32),
                       Vektor<std::size_t, 3>(8, 8, 128),
                       Vektor<double, 3>(6.0, 6.0, 6.0),
                       bunch.get(), i);
        
        
        assign(bunch->R, bunch->R + dt * bunch->P);
        
        // periodic shift
        double up = 4.0 * Physics::pi;
        for (std::size_t j = 0; j < bunch->getLocalNum(); ++j) {
            for (int d = 0; d < BL_SPACEDIM; ++d) {
                if ( bunch->R[j](d) > up )
                    bunch->R[j](d) = bunch->R[j](d) - up;
                else if ( bunch->R[j](d) < 0.0 )
                    bunch->R[j](d) = up + bunch->R[j](d);
            }
        }
        
        if ( myAmrOpal.maxLevel() > 0 )
            for (int k = 0; k <= myAmrOpal.finestLevel() && k < myAmrOpal.maxLevel(); ++k)
                myAmrOpal.regrid(k /*lbase*/, i /*time*/);
        else
            bunch->update();
        
        container_t rhs;
        container_t phi;
        container_t grad_phi;
        doSolve(myAmrOpal, bunch.get(), rhs, phi, grad_phi, geoms, rr, nLevels);
        
        bunch->GetGravity(bunch->E, grad_phi);
        
        /* epsilon_0 not used by Ulmer --> multiply it away */
        assign(bunch->P, bunch->P + dt * bunch->qm / bunch->mass * bunch->E * Physics::epsilon_0);
    }
}


int main(int argc, char *argv[]) {
    
    Ippl ippl(argc, argv);
    BoxLib::Initialize(argc,argv, false);
    
    Inform msg("TwoStream");
    

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("main");
    IpplTimings::startTimer(mainTimer);

    std::stringstream call;
    call << "Call: mpirun -np [#procs] " << argv[0]
         << " [#gridpoints x] [#gridpoints y] [#gridpoints z] [#levels]"
         << " [max. box size] [timstep] [#iterations]"
         << " [out: timing file name (optiona)]";
    
    if ( argc < 7 ) {
        msg << call.str() << endl;
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
    
    msg << "Particle test running with" << endl
        << "- #level     = " << nLevels << endl
        << "- grid       = " << nr << endl
        << "- max. size  = " << maxBoxSize << endl
        << "- time step  = " << dt << endl
        << "- #steps     = " << nIter << endl;
    
    doTwoStream(nr, nLevels, maxBoxSize,
                dt, nIter, msg);
    
    IpplTimings::stopTimer(mainTimer);

    IpplTimings::print();
    
    std::stringstream timefile;
    timefile << std::string(argv[0]) << "-timing-cores-"
             << std::setfill('0') << std::setw(6) << Ippl::getNodes()
             << "-threads-1.dat";
    
    if ( argc == 8 ) {
        timefile.str(std::string());
        timefile << std::string(argv[7]);
    }
    
    IpplTimings::print(timefile.str());
    
    return 0;
}