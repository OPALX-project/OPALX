
#ifndef POISSON_SOLVER_H_
#define POISSON_SOLVER_H_

//////////////////////////////////////////////////////////////
#include "Ippl.h"
#include "Algorithms/PartBunch.h"
//////////////////////////////////////////////////////////////

typedef UniformCartesian<3,double> Mesh_t;
typedef ParticleSpatialLayout<double,3>::SingleParticlePos_t Vector_t;
typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;
typedef Cell Center_t;
typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t> Field_t;

//use Barton and Nackman Trick to avoid virtual functions
//template <class T_Leaftype> 
class PoissonSolver
{
public:

  //T_Leaftype& asLeaf() { return static_cast<T_Leaftype&>(*this); }

  // given a charge-density field rho and a set of mesh spacings hr,
  // compute the scalar potential in open space
  //void computePotential(Field_t &rho, Vector_t hr) {asLeaf().computePotential(rho, hr);}
  //void computePotential(Field_t &rho, Vector_t hr, double zshift) {asLeaf().computePotential(&rho, hr, zshift);}

  // given a charge-density field rho and a set of mesh spacings hr,
  // compute the scalar potential in open space
  virtual void computePotential(Field_t &rho, Vector_t hr) = 0;
  virtual void computePotential(Field_t &rho, Vector_t hr, double zshift) = 0;

  virtual double getXRangeMin() = 0;
  virtual double getXRangeMax() = 0;
  virtual double getYRangeMin() = 0;
  virtual double getYRangeMax() = 0;
};

inline Inform &operator<<(Inform &os, const PoissonSolver &fs)
{
  return os << "";
}

#endif /* #ifndef MG_POISSON_SOLVER_H_ */
