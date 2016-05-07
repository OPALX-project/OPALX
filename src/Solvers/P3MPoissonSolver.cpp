// -*- C++ -*-
/***************************************************************************
 *
 *
 * P3MPoissonSolver.cc
 *
 *
 *
 *
 *
 *
 *
 ***************************************************************************/

// include files
#include "Solvers/P3MPoissonSolver.h"
#include "Algorithms/PartBunch.h"
#include "Physics/Physics.h"
#include <fstream>
//////////////////////////////////////////////////////////////////////////////
// a little helper class to specialize the action of the Green's function
// calculation.  This should be specialized for each dimension
// to the proper action for computing the Green's function.  The first
// template parameter is the full type of the Field to compute, and the second
// is the dimension of the data, which should be specialized.

#ifdef OPAL_NOCPLUSPLUS11_NULLPTR
#define nullptr NULL
#endif

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

////////////////////////////////////////////////////////////////////////////

// constructor


P3MPoissonSolver::P3MPoissonSolver(Mesh_t *mesh, FieldLayout_t *fl, std::string greensFunction, std::string bcz):
    mesh_m(mesh),
    layout_m(fl),
    mesh2_m(nullptr),
    layout2_m(nullptr),
    mesh3_m(nullptr),
    layout3_m(nullptr),
    mesh4_m(nullptr),
    layout4_m(nullptr),
    greensFunction_m(greensFunction)
{
    

}


P3MPoissonSolver::P3MPoissonSolver(PartBunch &beam, std::string greensFunction):
    mesh_m(&beam.getMesh()),
    layout_m(&beam.getFieldLayout()),
    mesh2_m(nullptr),
    layout2_m(nullptr),
    mesh3_m(nullptr),
    layout3_m(nullptr),
    mesh4_m(nullptr),
    layout4_m(nullptr),
    greensFunction_m(greensFunction)
{
   

}

////////////////////////////////////////////////////////////////////////////
// destructor
P3MPoissonSolver::~P3MPoissonSolver() {
}

////////////////////////////////////////////////////////////////////////////
// given a charge-density field rho and a set of mesh spacings hr,
// compute the electric potential from the image charge by solving
// the Poisson's equation

void P3MPoissonSolver::computePotential(Field_t &rho, Vector_t hr, double zshift) {
  
    
}


// given a charge-density field rho and a set of mesh spacings hr,
// compute the scalar potential in open space
void P3MPoissonSolver::computePotential(Field_t &rho, Vector_t hr) {
    
    
}



///////////////////////////////////////////////////////////////////////////
// calculate the FFT of the Green's function for the given field
void P3MPoissonSolver::greensFunction() {

    
}

void P3MPoissonSolver::mirrorRhoField(Field_t & ggrn2) {

}

Inform &P3MPoissonSolver::print(Inform &os) const {
    os << "* ************* F F T P o i s s o n S o l v e r ************************************ " << endl;
    os << "* h " << hr_m << '\n';
    os << "* ********************************************************************************** " << endl;
    return os;
}

/***************************************************************************
 * $RCSfile: P3MPoissonSolver.cc,v $   $Author: adelmann $
 * $Revision: 1.6 $   $Date: 2001/08/16 09:36:08 $
 ***************************************************************************/
