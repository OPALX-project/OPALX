
// -*- C++ -*-
/**************************************************************************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *

Arguments: nx ny nz qi
Example: ./FFTFieldSolver 4 4 4 1e-10 --commlib mpi --info 9 | tee output

 *************************************************************************************************************************************/

#include "Ippl.h"
#include <string>
#include <vector>
#include <limits>
#include "Particle/NoParticleCachingPolicy.h"

#include "mpi.h"

// some typedefs
typedef UniformCartesian<3, double> 								Mesh_t;
typedef ParticleSpatialLayout<double,3,Mesh_t,
        NoParticleCachingPolicy<double, 3, Mesh_t> >                playout_t;
typedef playout_t::SingleParticlePos_t    						    Vector_t;
typedef Cell                                                        Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t>                    FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>      					Field_t;
typedef Field<int, 3, Mesh_t, Center_t>                             IField_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>    					VField_t;
typedef Field<dcomplex, 3, Mesh_t, Center_t>    					CxField_t;
typedef FFT<RCTransform, 3, double>                  				FFT_t;
typedef IntCIC  													IntrplCIC_t;

template<unsigned int Dim>
struct SpecializedGreensFunction { };

template<>
struct SpecializedGreensFunction<3> {
    template<class T, class FT, class FT2>
        static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI) {
            grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
            grn = 1.0 / sqrt(grn);
            grn[0][0][0] = grn[0][0][1];
        }
};


template<class PL>
class ChargedParticles : public IpplParticleBase<PL> {

public:

    ParticleAttrib<double>     	Q;
    ParticleAttrib<Vector_t> 	EF;

    ChargedParticles(PL* pl, Vector_t nr, e_dim_tag decomp[3]) :
        IpplParticleBase<PL>(pl),
        nr_m(nr)
    {
        this->addAttribute(Q);

        for(int i = 0; i < 2 * 3; ++i) {
            bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
            getBConds()[i] = ParticleNoBCond;
        }

        for(int i=0; i<3; i++)
            decomp_m[i]=decomp[i];

        for(int d = 0;d<3;++d)
        {
            rmax_m[d] = 1.0;
            rmin_m[d] = 0.0;
            hr_m[d] = (rmax_m[d] - rmin_m[d]) / (nr_m[d] - 1.0);
        }
        getMesh().set_meshSpacing(&(hr_m[0]));
        getMesh().set_origin(rmin_m);
    }

    inline const Mesh_t& getMesh() const { return this->getLayout().getLayout().getMesh(); }
    inline Mesh_t& getMesh() { return this->getLayout().getLayout().getMesh(); }
    inline const FieldLayout_t& getFieldLayout() const {
        return dynamic_cast<FieldLayout_t&>( this->getLayout().getLayout().getFieldLayout());
    }
    inline FieldLayout_t& getFieldLayout() {
        return dynamic_cast<FieldLayout_t&>(this->getLayout().getLayout().getFieldLayout());
    }

    void setRMin(Vector_t x) { rmin_m = x; }
    void setHr(Vector_t x) { hr_m = x; }

    void update() {
        Inform msg("update");
        bounds(this->R, rmin_m, rmax_m);

        double stretch = 8*std::numeric_limits<double>::epsilon();
        Vector_t diff = rmax_m - rmin_m;
        rmax_m += stretch*diff;
        rmin_m -= stretch*diff;

        for(int d = 0;d<3;++d)
        {			
            hr_m[d] = (rmax_m[d] - rmin_m[d]) / (nr_m[d] - 1.0);
        }
        getMesh().set_meshSpacing(&(hr_m[0]));
        getMesh().set_origin(rmin_m);

        IpplParticleBase<PL>::update();
    }

    //FIXME: explain in more detail why we double fields and why we sometimes
    //       use "+1".
    void initialize_fields() {
        Inform msg("initialize_fields");

        density.initialize(getMesh(), getFieldLayout(), GuardCellSizes<3>(1), bc_m);
        gradient.initialize(getMesh(), getFieldLayout(), GuardCellSizes<3>(1), vbc_m);

        domain_m = getFieldLayout().getDomain();
        e_dim_tag decomp2[3];
        for(int d = 0; d < 3; ++d) {
            decomp2[d] = getFieldLayout().getRequestedDistribution(d);
        }

        // The FFT's require double-sized field sizes in order to (more closely
        // do not understand this ...)
        // simulate an isolated system.  The FFT of the charge density field, rho,
        // would otherwise mimic periodic boundary conditions, i.e. as if there were
        // several beams set a periodic distance apart.  The double-sized fields
        // alleviate this problem.
        for(int i = 0; i < 3; i++)
            domain2_m[i] = Index(2 * nr_m[i] + 1);

        // create double sized mesh and layout objects for the use in the FFT's
        mesh2_m = new Mesh_t(domain2_m);
        layout2_m = new FieldLayout_t(*mesh2_m, decomp2);
        density2.initialize(*mesh2_m, *layout2_m);

        NDIndex<3> tmpdomain;
        // Create the domain for the transformed (complex) fields.  Do this by
        // taking the domain from the doubled mesh, permuting it to the right, and
        // setting the 2nd dimension to have n/2 + 1 elements.
        domain3_m[0] = Index(2 * nr_m[3-1] + 1);
        domain3_m[1] = Index(nr_m[0] + 2);

        for(int i = 2; i < 3; ++i)
            domain3_m[i] = Index(2 * nr_m[i-1] + 1);

        // create mesh and layout for the new real-to-complex FFT's, for the
        // complex transformed fields
        mesh3_m = new Mesh_t(domain3_m);
        layout3_m = new FieldLayout_t(*mesh3_m, decomp2);
        densitytt.initialize(*mesh3_m, *layout3_m);
        greentt.initialize(*mesh3_m, *layout3_m);

        // create a domain used to indicate to the FFT's how to construct it's
        // temporary fields.  This is the same as the complex field's domain,
        // but permuted back to the left.
        tmpdomain = layout3_m->getDomain();
        for(int i = 0; i < 3; ++i)
            domainFFTConstruct_m[i] = tmpdomain[(i+1) % 3];

        // create the FFT object
        fft_m = new FFT_t(layout2_m->getDomain(), domainFFTConstruct_m);
        // these are fields that are used for calculating the Green's function.
        // they eliminate some calculation at each time-step.
        for(int i = 0; i < 3; ++i) {
            grnIField_m[i].initialize(*mesh2_m, *layout2_m);
            grnIField_m[i][domain2_m] = where(lt(domain2_m[i], nr_m[i]),
                    domain2_m[i] * domain2_m[i],
                    (2 * nr_m[i] - domain2_m[i]) *
                    (2 * nr_m[i] - domain2_m[i]));
        }
    }

    void computeElectricField() {
        density = 0.0;
        this->Q.scatter(this->density, this->R, IntrplCIC_t());

        // use grid of complex doubled in both dimensions
        // and store rho in lower left quadrant of doubled grid
        density2 = 0.0;

        density2[domain_m] = density[domain_m];

        // FFT double-sized charge density
        // we do a backward transformation so that we dont have to account for the normalization factor
        // that is used in the forward transformation of the IPPL FFT
        fft_m->transform(-1, density2, densitytt);

        // must be called if the mesh size has changed
        // have to check if we can do G with h = (1,1,1)
        // and rescale later

        Vector_t hrsq(hr_m * hr_m);
        SpecializedGreensFunction<3>::calculate(hrsq, density2, grnIField_m);
        fft_m->transform(-1, density2, greentt);

        densitytt *= greentt;

        // inverse FFT, rho2_m equals to the electrostatic gradient
        fft_m->transform(+1, densitytt, density2);
        // end convolution

        // back to physical grid
        // reuse the charge density field to store the electrostatic potential
        density[domain_m] = density2[domain_m];

        gradient = -Grad(density, gradient);
        EF.gather(gradient, this->R,  IntrplCIC_t());
    }

private:

    BConds<double, 3, Mesh_t, Center_t> bc_m;
    BConds<Vector_t, 3, Mesh_t, Center_t> vbc_m;

    Field_t density, density2;
    VField_t gradient;
    CxField_t densitytt, greentt;
    IField_t grnIField_m[3];

    Vektor<int,3> nr_m;
    Vector_t hr_m;
    Vector_t rmax_m;
    Vector_t rmin_m;

    NDIndex<3> domain_m, domain2_m, domain3_m, domainFFTConstruct_m;

    Mesh_t *mesh2_m, *mesh3_m;
    FieldLayout_t *layout2_m, *layout3_m;
    FFT_t *fft_m;

    e_dim_tag decomp_m[3];
};

int main(int argc, char *argv[]) {

    Ippl ippl(argc, argv);
    Inform msg(argv[0]);
    Inform msg2all(argv[0],INFORM_ALL_NODES);

    IpplTimings::TimerRef allTimer = IpplTimings::getTimer("AllTimer");
    IpplTimings::startTimer(allTimer);

    Vektor<int,3> nr;

    unsigned param = 1;

    nr = Vektor<int,3>(atoi(argv[param++]),atoi(argv[param++]),atoi(argv[param++]));

    double qi = atof(argv[param++]);

    e_dim_tag decomp[3];
    Mesh_t *mesh;
    FieldLayout_t *FL;
    ChargedParticles<playout_t>  *P;

    NDIndex<3> domain;
    for(int i=0; i<3; i++)
        domain[i] = domain[i] = Index(nr[i]);

    for (int d=0; d < 3; ++d)
        decomp[d] = PARALLEL;

    // create mesh and layout objects for this problem domain
    mesh          = new Mesh_t(domain);
    FL            = new FieldLayout_t(*mesh, decomp);
    playout_t* PL = new playout_t(*FL, *mesh);

    Vector_t hr(1.0);
    Vector_t rmin(0.0);
    Vector_t rmax(nr - 1.0); //XXX: rmax never used...


    P = new ChargedParticles<playout_t>(PL, nr, decomp);
    INFOMSG(P->getMesh() << endl);
    INFOMSG(P->getFieldLayout() << endl);
    msg << endl << endl; 
    Ippl::Comm->barrier();

    if(Ippl::myNode()==0)
    {
        P->create((nr[0]-1)*(nr[1]-1)*(nr[2]-1));
        for(int i = 0;i<(nr[0]-1)*(nr[1]-1)*(nr[2]-1);++i)
        {
            P->R[i](0) = double(0.5+i%(nr[0]-1));
            P->R[i](1) = double(0.5+i/(nr[0]-1));
            P->R[i](2) = double(0.5+i/(nr[0]-1)*(nr[1]-1));
            P->Q = qi;
        }	
    }

    P->update();

    msg << "initializing fields..";
    P->initialize_fields();
    msg << " done." << endl;

    msg << "computing electric field..";
    P->computeElectricField();
    msg << " done." << endl;

    IpplTimings::stopTimer(allTimer);
    IpplTimings::print();

    delete P;
    delete FL;
    delete mesh;

    return 0;
}

