#ifndef FIELD_SOLVER_H_
#define FIELD_SOLVER_H_
#include <ChargedParticles.hh>
#include <iostream>
#include <fstream>
#include <strstream>
//////////////////////////////////////////////////////////////
//forward declaration
template<class T, unsigned int Dim> class ChargedParticles;

template <class T, unsigned int Dim>
class FieldSolver 
{
  
public:
  // some useful typedefs
  typedef typename ChargedParticles<T,Dim>::VField_t VField_t;
  typedef typename ChargedParticles<T,Dim>::Field_t  Field_t;
  typedef typename ChargedParticles<T,Dim>::Vector_t Vector_t;

  // constructor and destructor
  FieldSolver();
  ~FieldSolver();

  // given a charge-density field rho and a set of mesh spacings hr,
  // compute the electric field and put in eg by solving the Poisson's equation
  virtual void computeForceField(ChargedParticles<T,Dim> *beam) = 0;

  virtual void testConvolution(ChargedParticles<T,Dim> *beam) = 0;
};

// needed if we're not using implicit inclusion
#include "FieldSolver.cc"

#endif

/***************************************************************************
 * $RCSfile: FieldSolver.hh,v $   $Author: candel $
 * $Revision: 1.2 $   $Date: 2001/08/08 21:12:10 $
 ***************************************************************************/

