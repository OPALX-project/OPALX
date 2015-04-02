/*
 * yfld-q1
 * Hal Finkel
 * Copyright (C) 2008 Yale Univeristy. All rights reserved.
 */

/*
 * This code was inspired by defrost.f90 version 0.17
 * which was written by Andrei Frolov <frolov@sfu.ca>
 * http://www.sfu.ca/physics/cosmology/defrost
 */

#define _XOPEN_SOURCE 500

#ifndef __cplusplus
#error standard C++ required
#endif

#include <cmath>
#include <complex>
#include <cfloat>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
using namespace std;

typedef complex<double> double_complex;

// This code requires the IPPL library from PSI: http://amas.web.psi.ch/
#include <Ippl.h>
#define Vector Vektor

/* **** Compile-Time Simulation Parameters **** */
const long n = 32;				// The size of the grid
/* **** END Compile-Time Simulation Parameters **** */

const long field_dim = 3;			// Dimension of the grid
const long nn = n/2+1;

const long fields = 2;				// Number of fields

// Create the geometry objects
const long cells = n-1;	
typedef UniformCartesian<field_dim> mesh_t;
typedef BConds<double, field_dim, mesh_t, Cell> bc_t;
typedef ParallelPeriodicFace<double, field_dim, mesh_t, Cell> face_t;
typedef CenteredFieldLayout<field_dim, mesh_t, Cell> fieldlayout_t;
typedef Field<double, field_dim, mesh_t, Cell> field_t;
typedef Field<bool, field_dim, mesh_t, Cell> fieldmask_t;

int main(int argc, char* argv[]) {
	long d;

	Ippl ippl(argc, argv);

	Ippl::Comm->barrier();

	NDIndex<field_dim> vertex_domain;
	for (d = 0; d < field_dim; ++d) {
		vertex_domain[d] = Index(n);
	}

	Index I = vertex_domain[0];
	Index J = vertex_domain[1];
	Index K = vertex_domain[2];

	mesh_t mesh(vertex_domain);

	fieldlayout_t layout(mesh, PARALLEL, PARALLEL, PARALLEL);

	bc_t bc;
	for (d = 0; d < field_dim; ++d) {
		bc[d] = new face_t(d);
	}

	const long os = 16;
	const long nos = n * os*os;

	// A dummy factor
	const double fl = 0.01;

	// Kernel domain and layout
	NDIndex<1> kernel_domain;
	Index kernelI(nos);
	kernel_domain[0] = kernelI;
	FieldLayout<1> kernel_layout(kernel_domain);
	Field<double, 1> kernel(kernel_layout);

	// Radial profile
	kernel[kernelI] = pow(kernelI*fl*(kernelI*fl*kernelI*fl + fl), fl) * exp(-(kernelI*fl)*(kernelI*fl));
	FFT<SineTransform,1,double> kernel_transform(kernel_domain);
	kernel_transform.transform(+1, kernel);

	kernel[kernelI] = fl * kernel[kernelI]/(kernelI + 1);

	// Initialize the convolution kernel (using linear interpolation of a radial profile)
	field_t conv_kernel(layout);

#define kk	(sqrt((I - (double) nn)*(I - (double) nn) + (J - (double) nn)*(J - (double) nn) + (K - (double) nn)*(K - (double) nn)) * (double) os)
#define kk_idx	(floor(kk) + 1)

	conv_kernel[I][J][K] = where(
			gt(kk_idx, 0),
			kernel[kk] + (kk-kk_idx)*(kernel[kk_idx+1]-kernel[kk]),
			(4.0*kernel[0] - kernel[1])/3.0
		);

#undef kk_idx
#undef kk

	return 0;
}
