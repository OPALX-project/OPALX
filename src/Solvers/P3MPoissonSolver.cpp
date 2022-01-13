//
// Class P3MPoissonSolver
//   This class contains methods for solving Poisson's equation for the
//   space charge portion of the calculation.
//
// Copyright (c) 2016, Benjamin Ulmer, ETH Zürich
// All rights reserved
//
// Implemented as part of the Master thesis
// "The P3M Model on Emerging Computer Architectures With Application to Microbunching"
// (http://amas.web.psi.ch/people/aadelmann/ETH-Accel-Lecture-1/projectscompleted/cse/thesisBUlmer.pdf)
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Solvers/P3MPoissonSolver.h"
#include "Algorithms/PartBunch.h"
#include "Particle/BoxParticleCachingPolicy.h"
#include "Particle/PairBuilder/HashPairBuilderPeriodic.h"
#include "Particle/PairBuilder/HashPairBuilderPeriodicParallel.h"
#include "Particle/PairBuilder/HashPairBuilder.h"
#include "Particle/PairBuilder/HashPairBuilderParallel.h"
//#include "Particle/PairBuilder/HashPairBuilderPeriodicParallel_globCHaining.h"
#include "Particle/PairBuilder/PairConditions.h"
#include "Structure/DataSink.h"
#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include <fstream>
#include <cmath>
//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.

//const double ke_m=1./(4.*M_PI*8.8e-14);
//const double ke_m=2.532638e8;

template<unsigned int Dim>
struct P3MGreensFunctionPeriodic { };

template<>
struct P3MGreensFunctionPeriodic<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI, double alpha,double eps) {
        double r;
        NDIndex<3> elem0=NDIndex<3>(Index(0,0), Index(0,0),Index(0,0));
        grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
        NDIndex<3> lDomain_m = grn.getLayout().getLocalNDIndex();
        NDIndex<3> elem;
        for (int i=lDomain_m[0].min(); i<=lDomain_m[0].max(); ++i) {
            elem[0]=Index(i,i);
            for (int j=lDomain_m[1].min(); j<=lDomain_m[1].max(); ++j) {
                elem[1]=Index(j,j);
                for (int k=lDomain_m[2].min(); k<=lDomain_m[2].max(); ++k) {
                    elem[2]=Index(k,k);
                    r = std::real(std::sqrt(grn.localElement(elem)));
                    if(elem==elem0) {
                        grn.localElement(elem) = 0 ;
                    }
                    else
                        grn.localElement(elem) = std::complex<double>(std::erf(alpha*r)/(r+eps));
                }
            }
        }
    }
};

template<unsigned int Dim>
struct P3MGreensFunction { };

template<>
struct P3MGreensFunction<3> {
    template<class T, class FT, class FT2>
    static void calculate(Vektor<T, 3> &hrsq, FT &grn, FT2 *grnI, double alpha, double eps) {
        grn = grnI[0] * hrsq[0] + grnI[1] * hrsq[1] + grnI[2] * hrsq[2];
        grn = erf(alpha*sqrt(grn))/(sqrt(grn)+eps);
        grn[0][0][0] = grn[0][0][1];
    }
};

template<class T>
struct ApplyField {
    ApplyField(T c, double r, double epsilon, double alpha, double ke_) : C(c), R(r), eps(epsilon), a(alpha), ke(ke_) {}
    void operator()(std::size_t i, std::size_t j, PartBunch &P,Vektor<double,3> &shift) const
    {
        Vector_t diff = P.R[i] - (P.R[j]+shift);
        double sqr = 0;

        for (unsigned d = 0; d<Dim; ++d)
            sqr += diff[d]*diff[d];

        //compute r with softening parameter, unsoftened r is obtained by sqrt(sqr)
        if(sqr!=0) {
            double r = std::sqrt(sqr+eps*eps);
            //for order two transition
            if (P.Q[i]!=0 && P.Q[j]!=0) {
                //compute potential energy
                //double phi =ke*(1.-std::erf(a*std::sqrt(sqr)))/r;

                //compute force
                Vector_t Fij = ke*C*(diff/std::sqrt(sqr))*((2.*a*std::exp(-a*a*sqr))/(std::sqrt(M_PI)*r)+(1.-std::erf(a*std::sqrt(sqr)))/(r*r));

                //Actual Force is F_ij multiplied by Qi*Qj
                //The electrical field on particle i is E=F/q_i and hence:
                P.Ef[i] -= P.Q[j]*Fij;
                P.Ef[j] += P.Q[i]*Fij;
                //update potential per particle
                //P.Phi[i] += P.Q[j]*phi;
                //P.Phi[j] += P.Q[i]*phi;
            }
        }
    }
    T C;
    double R;
    double eps;
    double a;
    double ke;
};


////////////////////////////////////////////////////////////////////////////

// constructor


P3MPoissonSolver::P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, double interaction_radius, double alpha, double eps, bool isTest):
    mesh_m(mesh),
    layout_m(fl),
    interaction_radius_m(interaction_radius),
    alpha_m(alpha),
    eps_m(eps),
    isTest_m(isTest)
{
    Inform msg("P3MPoissonSolver::P3MPoissonSolver ");
    if(isTest_m) {
        initFieldsTest();
        ke_m=2.532638e8;
    }
    else {
        initializeFields();
        ke_m = 1.0 / (4 * Physics::pi * Physics::epsilon_0);
    }

    GreensFunctionTimer_m = IpplTimings::getTimer("GreensFTotalP3M");
    ComputePotential_m = IpplTimings::getTimer("ComputePotentialP3M");
}

void P3MPoissonSolver::initFieldsTest() {
    
    Vector_t ll(-0.005);
    Vector_t ur(0.005);

    domain_m = layout_m->getDomain();
    for(unsigned int i = 0; i < 3; i++) {
        nr_m[i] = domain_m[i].length();
    }
    
    for (int i = 0; i < 3; i++)
        hr_m[i] = (ur[i] - ll[i]) / nr_m[i];

    mesh_m->set_meshSpacing(&(hr_m[0]));
    mesh_m->set_origin(ll);

    for (unsigned int i = 0; i < 2 * Dim; ++i) {
        if (Ippl::getNodes()>1) {
            bc_m[i] = new ParallelInterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new ParallelPeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new ParallelPeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
        else {
            bc_m[i] = new InterpolationFace<double, Dim, Mesh_t, Center_t>(i);
            //std periodic boundary conditions for gradient computations etc.
            vbc_m[i] = new PeriodicFace<Vector_t, Dim, Mesh_t, Center_t>(i);
            bcp_m[i] = new PeriodicFace<double, Dim, Mesh_t, Center_t>(i);
        }
    }

    rhocmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));
    grncmpl_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1));

    rho_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bc_m);
    phi_m.initialize(*mesh_m, *layout_m, GuardCellSizes<Dim>(1), bcp_m);
    eg_m.initialize (*mesh_m, *layout_m, GuardCellSizes<Dim>(1), vbc_m);

    bool compressTemps = true;
    if (fft_m)
        fft_m.reset();

    // create the FFT object
    fft_m = std::unique_ptr<FFTC_t>(new FFTC_t(layout_m->getDomain(), compressTemps));
    fft_m->setDirectionName(+1, "forward");
    fft_m->setDirectionName(-1, "inverse");

}

void P3MPoissonSolver::initializeFields() {

    domain_m = layout_m->getDomain();

    // For efficiency in the FFT's, we can use a parallel decomposition
    // which can be serial in the first dimension.
    e_dim_tag decomp[3];
    e_dim_tag decomp2[3];
    for(int d = 0; d < 3; ++ d) {
        decomp[d] = layout_m->getRequestedDistribution(d);
        decomp2[d] = layout_m->getRequestedDistribution(d);
    }

    // The FFT's require double-sized field sizes in order to
    // simulate an isolated system.  The FFT of the charge density field, rho,
    // would otherwise mimic periodic boundary conditions, i.e. as if there were
    // several beams set a periodic distance apart.  The double-sized fields
    // alleviate this problem.
    for (int i = 0; i < 3; ++ i) {
        hr_m[i] = mesh_m->get_meshSpacing(i);
        nr_m[i] = domain_m[i].length();
        domain2_m[i] = Index(2 * nr_m[i] + 1);
    }

    for (int i = 0; i < 2 * 3; ++ i) {
        bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
        vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
    }

    // create double sized mesh and layout objects for the use in the FFT's
    mesh2_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain2_m));
    layout2_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh2_m, decomp));

    rho2_m.initialize(*mesh2_m, *layout2_m);

    NDIndex<3> tmpdomain;
    // Create the domain for the transformed (complex) fields.  Do this by
    // taking the domain from the doubled mesh, permuting it to the right, and
    // setting the 2nd dimension to have n/2 + 1 elements.
    domain3_m[0] = Index(2 * nr_m[2] + 1);
    domain3_m[1] = Index(nr_m[0] + 2);
    domain3_m[2] = Index(2 * nr_m[1] + 1);

    // create mesh and layout for the new real-to-complex FFT's, for the
    // complex transformed fields
    mesh3_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain3_m));
    layout3_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh3_m, decomp2));

    rho2tr_m.initialize(*mesh3_m, *layout3_m);
    imgrho2tr_m.initialize(*mesh3_m, *layout3_m);
    grntr_m.initialize(*mesh3_m, *layout3_m);

    // helper field for sin
    greentr_m.initialize(*mesh3_m, *layout3_m);

    for (int i = 0; i < 3; ++ i) {
        domain4_m[i] = Index(nr_m[i] + 2);
    }
    mesh4_m = std::unique_ptr<Mesh_t>(new Mesh_t(domain4_m));
    layout4_m = std::unique_ptr<FieldLayout_t>(new FieldLayout_t(*mesh4_m, decomp));

    tmpgreen_m.initialize(*mesh4_m, *layout4_m);

    // create a domain used to indicate to the FFT's how to construct it's
    // temporary fields.  This is the same as the complex field's domain,
    // but permuted back to the left.
    tmpdomain = layout3_m->getDomain();
    for (int i = 0; i < 3; ++ i)
        domainFFTConstruct_m[i] = tmpdomain[(i+1) % 3];

    // create the FFT object
    fftrc_m = std::unique_ptr<FFTRC_t>(new FFTRC_t(layout2_m->getDomain(), domainFFTConstruct_m));

    // these are fields that are used for calculating the Green's function.
    // they eliminate some calculation at each time-step.
    for (int i = 0; i < 3; ++ i) {
        grnIField_m[i].initialize(*mesh2_m, *layout2_m);
        grnIField_m[i][domain2_m] = where(lt(domain2_m[i], nr_m[i]),
                                          domain2_m[i] * domain2_m[i],
                                          (2 * nr_m[i] - domain2_m[i]) *
                                          (2 * nr_m[i] - domain2_m[i]));
    }
}


////////////////////////////////////////////////////////////////////////////
// destructor
P3MPoissonSolver::~P3MPoissonSolver() {
}

void P3MPoissonSolver::calculatePairForcesPeriodic(PartBunchBase<double, 3> *bunch) {
    if (interaction_radius_m>0){
        if (Ippl::getNodes() > 1) {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilderPeriodicParallel<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius_m), ApplyField<double>(-1,interaction_radius_m,eps_m,alpha_m,ke_m),extend_l, extend_r);
        }
        else {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilderPeriodic<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius_m), ApplyField<double>(-1,interaction_radius_m,eps_m,alpha_m,ke_m),extend_l, extend_r);
        }
    }

}

void P3MPoissonSolver::calculatePairForces(PartBunchBase<double, 3> *bunch) {
    if (interaction_radius_m>0){
        if (Ippl::getNodes() > 1) {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilderParallel<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius_m), ApplyField<double>(-1,interaction_radius_m,eps_m,alpha_m,ke_m));
        }
        else {
            PartBunch &tmpBunch = *(dynamic_cast<PartBunch*>(bunch));
            HashPairBuilder<PartBunch> HPB(tmpBunch);
            HPB.for_each(RadiusCondition<double, Dim>(interaction_radius_m), ApplyField<double>(-1,interaction_radius_m,eps_m,alpha_m,ke_m));
        }
    }
}

void P3MPoissonSolver::calculateGridForces(PartBunchBase<double, 3> *bunch){

    Inform msg ("calculateGridForces ");
    Vector_t l,h;

    // (1) scatter charge to charge density grid and transform to fourier space
    rho_m[domain_m]=0;
    bunch->Q.scatter(rho_m, bunch->R, IntrplCIC_t());

    rhocmpl_m[domain_m] = rho_m[domain_m]/(hr_m[0]*hr_m[1]*hr_m[2]);

    fft_m->transform("inverse", rhocmpl_m);

    // (2) compute Greens function in real space and transform to fourier space
    IField_t grnIField[3];

    // This loop stores in grnIField_m[i] the index of the ith dimension mirrored at the central axis. e.g.
    // grnIField_m[0]=[(0 1 2 3 ... 3 2 1) ; (0 1 2 3 ... 3 2 1; ...)]
    for (int i = 0; i < 3; ++i) {
        grnIField[i].initialize(*mesh_m, *layout_m);
        grnIField[i][domain_m] = where(lt(domain_m[i], nr_m[i]/2),
                                       domain_m[i] * domain_m[i],
                                       (nr_m[i]-domain_m[i]) *
                                       (nr_m[i]-domain_m[i]));
    }
    Vector_t hrsq(hr_m * hr_m);
    P3MGreensFunctionPeriodic<3>::calculate(hrsq, grncmpl_m, grnIField, alpha_m, eps_m);

    //transform G -> Ghat and store in grncmpl_m
    fft_m->transform("inverse", grncmpl_m);
    //multiply in fourier space and obtain PhiHat in rhocmpl_m
    rhocmpl_m *= grncmpl_m;

    // (3) Backtransformation: compute potential field in real space and E=-Grad Phi
    //compute electrostatic potential Phi in real space by FFT PhiHat -> Phi and store it in rhocmpl_m
    fft_m->transform("forward", rhocmpl_m);

    //take only the real part and store in phi_m (has periodic bc instead of interpolation bc)
    phi_m = real(rhocmpl_m)*hr_m[0]*hr_m[1]*hr_m[2];
    phi_m *= ke_m;
    //dumpVTKScalar( phi_m, this,it, "Phi_m") ;


    //compute Electric field on the grid by -Grad(Phi) store in eg_m
    eg_m = -Grad1Ord(phi_m, eg_m);

    //interpolate the electric field to the particle positions
    bunch->Ef.gather(eg_m, bunch->R,  IntrplCIC_t());
    //interpolate electrostatic potenital to the particle positions
    bunch->Phi.gather(phi_m, bunch->R, IntrplCIC_t());
}



////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void P3MPoissonSolver::computePotential(Field_t &/*rho*/, Vector_t /*hr*/, double /*zshift*/) {


}

void P3MPoissonSolver::computeAvgSpaceChargeForces(PartBunchBase<double, 3> *bunch) {

    Inform m("computeAvgSpaceChargeForces ");

    const double N =  static_cast<double>(bunch->getTotalNum());
    double locAvgEf[Dim]={};
    for (unsigned i=0; i<bunch->getLocalNum(); ++i) {
        locAvgEf[0]+=std::abs(bunch->Ef[i](0));
        locAvgEf[1]+=std::abs(bunch->Ef[i](1));
        locAvgEf[2]+=std::abs(bunch->Ef[i](2));
    }

    reduce(&(locAvgEf[0]), &(locAvgEf[0]) + Dim,
           &(globSumEf_m[0]), OpAddAssign());

    //    m << "globSumEF = " << globSumEf_m[0] << "\t" << globSumEf_m[1] << "\t" << globSumEf_m[2] << endl;

    avgEF_m[0]=globSumEf_m[0]/N;
    avgEF_m[1]=globSumEf_m[1]/N;
    avgEF_m[2]=globSumEf_m[2]/N;

}


void P3MPoissonSolver::applyConstantFocusing(PartBunchBase<double, 3> *bunch, double f, double r){
    Vektor<double,Dim> beam_center(0,0,0);
    Vector_t Rrel;
    double scFoc = std::sqrt(dot(avgEF_m,avgEF_m));
    for (unsigned i=0; i<bunch->getLocalNum(); ++i) {
        Rrel=bunch->R[i] - beam_center;
        bunch->Ef[i] += Rrel/r*f*scFoc;
    }

}

// given a charge-density field rho and a set of mesh spacings hr,
// compute the scalar potential in open space
void P3MPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {
    IpplTimings::startTimer(ComputePotential_m);

    // use grid of complex doubled in both dimensions
    // and store rho in lower left quadrant of doubled grid
    rho2_m = 0.0;

    rho2_m[domain_m] = rho[domain_m];

    // needed in greens function
    hr_m = hr;

    // FFT double-sized charge density
    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    fftrc_m->transform(-1, rho2_m, rho2tr_m);

    // must be called if the mesh size has changed
    // have to check if we can do G with h = (1,1,1)
    // and rescale later
    IpplTimings::startTimer(GreensFunctionTimer_m);
    greensFunction();
    IpplTimings::stopTimer(GreensFunctionTimer_m);
    // multiply transformed charge density
    // and transformed Green function
    // Don't divide by (2*nx_m)*(2*ny_m), as Ryne does;
    // this normalization is done in POOMA's fft routine.
    rho2tr_m *= grntr_m;

    // inverse FFT, rho2_m equals to the electrostatic potential
    fftrc_m->transform(+1, rho2tr_m, rho2_m);
    // end convolution

    // back to physical grid
    // reuse the charge density field to store the electrostatic potential
    rho[domain_m] = rho2_m[domain_m];


    rho *= hr[0] * hr[1] * hr[2];
    IpplTimings::stopTimer(ComputePotential_m);
}

void P3MPoissonSolver::greensFunction() {

    Vector_t hrsq(hr_m * hr_m);
    P3MGreensFunction<3>::calculate(hrsq, rho2_m, grnIField_m, alpha_m, eps_m);
    // Green's function calculation complete at this point.
    // The next step is to FFT it.
    // FFT of Green's function

    // we do a backward transformation so that we dont have to account for the normalization factor
    // that is used in the forward transformation of the IPPL FFT
    fftrc_m->transform(-1, rho2_m, grntr_m);
}


void P3MPoissonSolver::compute_temperature(PartBunchBase<double, 3> *bunch) {

    Inform m("compute_temperature ");
    // double loc_temp[Dim]={0.0,0.0,0.0};
    double loc_avg_vel[Dim]={0.0,0.0,0.0};
    double avg_vel[Dim]={0.0,0.0,0.0};

    for(unsigned long k = 0; k < bunch->getLocalNum(); ++k) {
        for(unsigned i = 0; i < Dim; i++) {
            loc_avg_vel[i]   += bunch->P[k](i);
        }
    }
    reduce(&(loc_avg_vel[0]), &(loc_avg_vel[0]) + Dim,
           &(avg_vel[0]), OpAddAssign());

    const double N =  static_cast<double>(bunch->getTotalNum());
    avg_vel[0]=avg_vel[0]/N;
    avg_vel[1]=avg_vel[1]/N;
    avg_vel[2]=avg_vel[2]/N;

    m << "avg_vel[0]= " << avg_vel[0] << " avg_vel[1]= " << avg_vel[1]  << " avg_vel[2]= " << avg_vel[2] << endl;

    /*
      for(unsigned long k = 0; k < bunch.getLocalNum(); ++k) {
      for(unsigned i = 0; i < Dim; i++) {
      loc_temp[i]   += (bunch.P[k](i)-avg_vel[i])*(bunch.P[k](i)-avg_vel[i]);
      }
      }
      reduce(&(loc_temp[0]), &(loc_temp[0]) + Dim,
      &(temperature[0]), OpAddAssign());
      temperature[0]=temperature[0]/N;
      temperature[1]=temperature[1]/N;
      temperature[2]=temperature[2]/N;
    */
}

void P3MPoissonSolver::test(PartBunchBase<double, 3> *bunch) {
    Inform msg("P3MPoissonSolver::test ");

    // set special conditions for this test
    const double mi = 1.0;
    const double qi = -1.0;
    const double qom = qi/mi;
    const double beam_radius = 0.001774;
    const double f = 1.5;
    const double dt = bunch->getdT();

   
    OpalData *opal = OpalData::getInstance();
    DataSink *ds = opal->getDataSink();

    Vector_t FDext[6];

    bunch->Q = qi;
    bunch->M = mi;

    ds->dumpSDDS(bunch, FDext, 0);

    //initFieldsTest();

    for (int i=0; i<3; i++) {
        extend_r[i] =  hr_m[i]*nr_m[i]/2;
        extend_l[i] = -hr_m[i]*nr_m[i]/2;
    }

    msg << *this << endl;

    // calculate initial space charge forces
    bunch->update();
    calculateGridForces(bunch);
    calculatePairForcesPeriodic(bunch);


    //avg space charge forces for constant focusing
    computeAvgSpaceChargeForces(bunch);

    for (int it=0; it<1000; it++) {

        // advance the particle positions
        // basic leapfrogging timestep scheme.  velocities are offset
        // by half a timestep from the positions.

        assign(bunch->R, bunch->R + dt * bunch->P);

        bunch->update();

        calculateGridForces(bunch);
        calculatePairForcesPeriodic(bunch);
        applyConstantFocusing(bunch,f,beam_radius);

        assign(bunch->P, bunch->P + dt * qom * bunch->Ef);

        ds->dumpSDDS(bunch, FDext, it+1);
        bunch->incrementT();
        msg << "Finished iteration " << it+1 << endl;
    }
}


Inform &P3MPoissonSolver::print(Inform &os) const {
    os << "* ************* P 3 M - P o i s s o n S o l v e r *************** " << endl;
    os << "* h        " << hr_m << '\n';
    os << "* RC       " << interaction_radius_m << '\n';
    os << "* ALPHA    " << alpha_m << '\n';
    os << "* EPSILON  " << eps_m << '\n';
    os << "* Extend L " << extend_l << '\n';
    os << "* Extend R " << extend_r << '\n';
    os << "* nr       " << nr_m << '\n';
    os << "* *************************************************************** " << endl;
    return os;
}
