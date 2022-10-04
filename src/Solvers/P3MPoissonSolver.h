//
// Class P3MPoissonSolver
//   This class contains methods for solving Poisson's equation for the
//   space charge portion of the calculation including collisions.
//
// Copyright (c) 2016, Benjamin Ulmer, ETH Zürich
//               2022, Sriramkrishnan Muralikrishnan, PSI
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

//////////////////////////////////////////////////////////////
#include "PoissonSolver.h"

#include "FFT/FFT.h"

#include <memory>

//#include "Algorithms/PartBunchBase.h"

template <class T, unsigned Dim>
class PartBunchBase;

//////////////////////////////////////////////////////////////

class P3MPoissonSolver : public PoissonSolver {
public:

    typedef FFT<RCTransform, 3, double> FFT_t;

    // constructor and destructor
    P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, 
                     double interaction_radius, 
                     double alpha,
                     std::string greensFunction);

    ~P3MPoissonSolver();

    void initializeFields();
    
    void calculatePairForces(PartBunchBase<double, 3> *bunch, double gammaz) override;

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    void computePotential(Field_t &rho, Vector_t hr) override;

    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential with image charges at  -z
    void computePotential(Field_t &rho, Vector_t hr, double zshift) override;

    void greensFunction();

    void integratedGreensFunction();

    void mirrorRhoField();

    void test(PartBunchBase<double, 3> */*bunch*/) override {};

    double getXRangeMin(unsigned short /*level*/) override {return 1.0;}
    double getXRangeMax(unsigned short /*level*/) override {return 1.0;}
    double getYRangeMin(unsigned short /*level*/) override {return 1.0;}
    double getYRangeMax(unsigned short /*level*/) override {return 1.0;}
    double getZRangeMin(unsigned short /*level*/) override {return 1.0;}
    double getZRangeMax(unsigned short /*level*/) override {return 1.0;}
    double getinteractionRadius() const override {return interaction_radius_m;}

    Inform &print(Inform &os) const;
    
private:
    // original charge density
    Field_t rho_m;
    // rho2_m is the charge-density field with mesh doubled in each dimension
    Field_t rho2_m;

    // rho2tr_m is the Fourier transformed charge-density field
    // domain3_m and mesh3_mp are used
    CxField_t rho2tr_m;
    
    // Fields used to eliminate excess calculation in greensFunction()
    // mesh2_mp and layout2_mp are used
    IField_t grnIField_m[3];

    // grntr_m is the Fourier transformed Green's function
    // domain3_m and mesh3_mp are used
    CxField_t grntr_m;

    // the FFT object
    std::unique_ptr<FFT_t> fft_mp;

    // mesh and layout objects for rho_m
    Mesh_t *mesh_mp;
    FieldLayout_t *layout_mp;

    // mesh and layout objects for rho2_m
    std::unique_ptr<Mesh_t> mesh2_mp;
    std::unique_ptr<FieldLayout_t> layout2_mp;

    std::unique_ptr<Mesh_t> mesh3_mp;
    std::unique_ptr<FieldLayout_t> layout3_mp;

    // mesh and layout for integrated greens function
    std::unique_ptr<Mesh_t> mesh4_mp;
    std::unique_ptr<FieldLayout_t> layout4_mp;

    // tmp
    Field_t tmpgreen_m;

    // domains for the various fields
    NDIndex<3> domain_m;              // original domain, gridsize
    NDIndex<3> domain2_m;             // doubled gridsize (2*Nx,2*Ny,2*Nz)
    NDIndex<3> domain3_m;             // field for the complex values of the RC transformation
    NDIndex<3> domain4_m;             // domain for tmp in integrated Greens function
    NDIndex<3> domainFFTConstruct_m;  // domain for output of FFT 

    double interaction_radius_m;
    double alpha_m;

    Vector_t hr_m;
    Vektor<int, 3> nr_m;
    double ke_m;

    bool integratedGreens_m;

    IpplTimings::TimerRef GreensFunctionTimer_m;
    IpplTimings::TimerRef ComputePotential_m;
    IpplTimings::TimerRef CalculatePairForces_m;

};

inline Inform &operator<<(Inform &os, const P3MPoissonSolver &fs) {
    return fs.print(os);
}



#endif
