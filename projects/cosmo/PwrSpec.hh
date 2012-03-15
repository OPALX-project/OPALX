#ifndef PWR_SPEC_H
#define PWR_SPEC_H

#include <ChargedParticles.hh>
#include <iostream>
#include <fstream>
#include <strstream>
#include "AppTypes/dcomplex.h"
#include "FFT/FFT.h"

template <class T, unsigned int Dim>
class PwrSpec {

  typedef typename ChargedParticles<T,Dim>::Center_t      Center_t;
  typedef typename ChargedParticles<T,Dim>::Mesh_t        Mesh_t;
  typedef typename ChargedParticles<T,Dim>::FieldLayout_t FieldLayout_t;
  typedef typename ChargedParticles<T,Dim>::Vector_t      Vector_t;
  typedef typename ChargedParticles<T,Dim>::IntrplCIC_t   IntrplCIC_t;

  typedef Field<dcomplex, Dim, Mesh_t, Center_t>    CxField_t;
  typedef FFT<CCTransform, Dim, T>                  FFT_t;

  typedef Field<T, 1, Mesh_t, Center_t>     Field1D_t;

  // Fourier transformed field 
  CxField_t fC_m;
  
  // the pwr spectrum 
  CxField_t pwrSpec_m;

  // the FFT object
  FFT_t *fft_m;

  // mesh and layout objects for rho_m
  Mesh_t *mesh_m;
  FieldLayout_t *layout_m;
      
  /// global domain for the various fields
  NDIndex<Dim> gDomain_m;  
  /// local domain for the various fields
  NDIndex<Dim> lDomain_m;             

public:

  PwrSpec(ChargedParticles<T,Dim> *beam);
  ~PwrSpec(); 
  
  void cacl1D(ChargedParticles<T,Dim> *b);
  void save1DSpectra(string fn);

};
#include "PwrSpec.cc"
#endif


