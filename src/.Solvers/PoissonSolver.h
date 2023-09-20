#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Algorithms/PBunchDefs.h"

#include "Field/Field.h"

//////////////////////////////////////////////////////////////
template <class T, unsigned Dim>
class PartBunch;

class PoissonSolver {
protected:
    typedef Field<int, 3, Mesh_t, Center_t> IField_t;
    typedef Field<std::complex<double>, 3, Mesh_t, Center_t> CxField_t;

public:
    // given a charge-density field rho and a set of mesh spacings hr,
    // compute the scalar potential in open space
    virtual void computePotential(Field_t& rho, Vector_t<double, 3> hr) = 0;

    virtual void computePotential(Field_t& rho, Vector_t<double, 3> hr, double zshift) = 0;

    virtual double getXRangeMin(unsigned short level = 0) = 0;
    virtual double getXRangeMax(unsigned short level = 0) = 0;
    virtual double getYRangeMin(unsigned short level = 0) = 0;
    virtual double getYRangeMax(unsigned short level = 0) = 0;
    virtual double getZRangeMin(unsigned short level = 0) = 0;
    virtual double getZRangeMax(unsigned short level = 0) = 0;
    virtual void test(PartBunch_t* bunch)    = 0;
    virtual ~PoissonSolver(){};

    virtual void resizeMesh(
        Vector_t<double, 3>& /*origin*/, Vector_t<double, 3>& /*hr*/, const Vector_t<double, 3>& /*rmin*/, const Vector_t<double, 3>& /*rmax*/,
        double /*dh*/){};
};

inline Inform& operator<<(Inform& os, const PoissonSolver& /*fs*/) {
    return os << "";
}
#endif
