// -*- C++ -*-

/***************************************************************************
 *
 * 
 * FFTPoissonSolverPeriodic.hh 
 * 
 * Periodic BC in x,y and z. 
 * 
 *
 ***************************************************************************/


////////////////////////////////////////////////////////////////////////////
// This class contains methods for solving Poisson's equation for the 
// space charge portion of the calculation.
////////////////////////////////////////////////////////////////////////////

#ifndef FFT_POISSON_SOLVER_PERIODIC_H_
#define FFT_POISSON_SOLVER_PERIODIC_H_

//////////////////////////////////////////////////////////////

#include "FieldSolver.hh"
#include "AppTypes/dcomplex.h"
#include "FFT/FFT.h"
//////////////////////////////////////////////////////////////

template <class T, unsigned int Dim>
class FFTPoissonSolverPeriodic : public FieldSolver<T,Dim> 
{
public:
    // some useful typedefs
  typedef typename ChargedParticles<T,Dim>::Center_t      Center_t;
  typedef typename ChargedParticles<T,Dim>::Mesh_t        Mesh_t;
  typedef typename ChargedParticles<T,Dim>::FieldLayout_t FieldLayout_t;
  typedef typename ChargedParticles<T,Dim>::Vector_t      Vector_t;
  typedef typename ChargedParticles<T,Dim>::IntrplCIC_t   IntrplCIC_t;

  typedef Field<dcomplex, Dim, Mesh_t, Center_t>    CxField_t;
  typedef FFT<CCTransform, Dim, T>                  FFT_t;
public:
    // constructor and destructor
  FFTPoissonSolverPeriodic(ChargedParticles<T,Dim> *beam);
    
  virtual ~FFTPoissonSolverPeriodic();
    
  // given a density field rho 
  // compute the force field 
  virtual void computeForceField(ChargedParticles<T,Dim> *beam); 
    
  // precompute the green's function for a Poisson problem and put it in green_m
  void greensFunction();

  virtual void testConvolution(ChargedParticles<T,Dim> *beam);
private:    

  void saveGreensHatFunction();
    
  // green_m  is the Fourier transformed greens function
  CxField_t green_m;
  
  // rho_m is the Fourier transformed density field
  CxField_t rho_m;
    
   // the FFT object
  FFT_t *fft_m;

  // mesh and layout objects for rho_m
  Mesh_t *mesh_m;
  FieldLayout_t *layout_m;
      
  /// global domain for the various fields
  NDIndex<Dim> gDomain_m;  
  /// local domain for the various fields
  NDIndex<Dim> lDomain_m;             
};

// needed if we're not using implicit inclusion
#include "FFTPoissonSolverPeriodic.cc"

#endif

/***************************************************************************
 * $RCSfile: FFTPoissonSolverPeriodic.hh,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2001/08/08 11:21:48 $
 ***************************************************************************/
