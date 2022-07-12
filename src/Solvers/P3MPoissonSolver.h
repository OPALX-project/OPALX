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
#ifndef P3M_POISSON_SOLVER_H_
#define P3M_POISSON_SOLVER_H_
const unsigned Dim = 3;

#ifdef dontOPTIMIZE_FIELD_ASSIGNMENT
#define FIELDASSIGNOPTIMIZATION __attribute__((optimize(0)))
#else
#define FIELDASSIGNOPTIMIZATION
#endif

#include <memory>
//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"

#include "FFT/FFT.h"

//#include "Algorithms/PartBunchBase.h"

template <class T, unsigned Dim>
class PartBunchBase;

//////////////////////////////////////////////////////////////

class P3MPoissonSolver : public PoissonSolver {
public:

    typedef FFT<CCTransform, 3, double>              FFTC_t;
    typedef FFT<RCTransform, 3, double>              FFTRC_t;

    // constructor and destructor
    P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, 
                     double interaction_radius, 
                     double alpha, double eps, bool isTest);

    ~P3MPoissonSolver();

    void initFieldsTest();
    
    void initializeFields();

    void calculateGridForces(PartBunchBase<double, 3> *bunch);

    void calculatePairForcesPeriodic(PartBunchBase<double, 3> *bunch);
    
    void calculatePairForces(PartBunchBase<double, 3> *bunch, double gammaz);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential with image charges at  -z
    void computePotential(Field_t &rho, Vector_t hr, double zshift);

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    void computePotential(Field_t &rho, Vector_t hr);

    void greensFunction();

    void integratedGreensFunction();

    void mirrorRhoField();

    void applyConstantFocusing(PartBunchBase<double, 3> *bunch, double f, double r);
    void test(PartBunchBase<double, 3> *bunch);

    double getXRangeMin(unsigned short /*level*/) {return 1.0;}
    double getXRangeMax(unsigned short /*level*/) {return 1.0;}
    double getYRangeMin(unsigned short /*level*/) {return 1.0;}
    double getYRangeMax(unsigned short /*level*/) {return 1.0;}
    double getZRangeMin(unsigned short /*level*/) {return 1.0;}
    double getZRangeMax(unsigned short /*level*/) {return 1.0;}
    double getinteractionRadius() override {return interaction_radius_m;}
    bool isTest() override {return isTest_m;}
    void setinteractionRadius(double r) {interaction_radius_m = r;}
    void setAlpha(double alpha) {alpha_m = alpha;}

    void computeAvgSpaceChargeForces(PartBunchBase<double, 3> *bunch);
    void compute_temperature(PartBunchBase<double, 3> *bunch);
    Inform &print(Inform &os) const;
    
private:

    BConds<double, Dim, Mesh_t, Center_t> bc_m;
    BConds<double, Dim, Mesh_t, Center_t> bcp_m;
    BConds<Vector_t, Dim, Mesh_t, Center_t> vbc_m;

    // rho_m is the charge-density field with mesh doubled in each dimension
    Field_t rho_m;
    Field_t phi_m;
    Field_t rho2_m;

    VField_t eg_m;

    // real field with layout of complex field: domain3_m
    Field_t greentr_m;
    
    // rho2tr_m is the Fourier transformed charge-density field
    // domain3_m and mesh3_ are used
    CxField_t rho2tr_m;
    CxField_t imgrho2tr_m;
    
    // Fields used to eliminate excess calculation in greensFunction()
    // mesh2_m and layout2_m are used
    IField_t grnIField_m[3];

    CxField_t rhocmpl_m;
    CxField_t grncmpl_m;

    // grntr_m is the Fourier transformed Green's function
    // domain3_m and mesh3_ are used
    CxField_t grntr_m;

    // the FFT object
    std::unique_ptr<FFTC_t> fft_m;
    std::unique_ptr<FFTRC_t> fftrc_m;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_m;
    FieldLayout_t *layout_m;

    // mesh and layout objects for rho2_m
    std::unique_ptr<Mesh_t> mesh2_m;
    std::unique_ptr<FieldLayout_t> layout2_m;

    //
    std::unique_ptr<Mesh_t> mesh3_m;
    std::unique_ptr<FieldLayout_t> layout3_m;

    // mesh and layout for integrated greens function
    std::unique_ptr<Mesh_t> mesh4_m;
    std::unique_ptr<FieldLayout_t> layout4_m;

    // tmp
    Field_t tmpgreen_m;

    // domains for the various fields
    NDIndex<3> domain_m;             // original domain, gridsize
    // mesh and gridsize defined outside of P3M class, given as
    NDIndex<3> domain2_m;            // doubled gridsize (2*Nx,2*Ny,2*Nz)
    NDIndex<3> domain3_m;            // field for the complex values of the RC transformation
    NDIndex<3> domain4_m;
    // (2*Nx,Ny,2*Nz)
    NDIndex<3> domainFFTConstruct_m;


    double interaction_radius_m;
    double alpha_m;
    double eps_m;
    bool isTest_m;

    Vector_t hr_m;
    Vektor<int, 3> nr_m;
    double ke_m;

    // for tests
    Vektor<double,Dim> avgEF_m;
    double globSumEf_m[Dim];


    IpplTimings::TimerRef GreensFunctionTimer_m;
    IpplTimings::TimerRef ComputePotential_m;


public:
    Vektor<double,3> extend_l;
    Vektor<double,3> extend_r;


};

inline Inform &operator<<(Inform &os, const P3MPoissonSolver &fs) {
    return fs.print(os);
}



#endif
