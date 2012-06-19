#ifndef PWR_SPEC_H
#define PWR_SPEC_H

#include <ChargedParticles/ChargedParticles.hh>
#include <iostream>
#include <fstream>
#include <strstream>

// #include "AppTypes/dcomplex.h"
// #include "FFT/FFT.h"

template <class T, unsigned int Dim>
class PwrSpec {

    typedef typename ChargedParticles<T,Dim>::Center_t      Center_t;
    typedef typename ChargedParticles<T,Dim>::Mesh_t        Mesh_t;
    typedef typename ChargedParticles<T,Dim>::FieldLayout_t FieldLayout_t;
    typedef typename ChargedParticles<T,Dim>::Vector_t      Vector_t;
    typedef typename ChargedParticles<T,Dim>::IntrplCIC_t   IntrplCIC_t;

    T*spectra1D_m;
    int *Nk_m;

    FieldLayout_t *layout_m;

    /// global domain for the various fields
    NDIndex<Dim> gDomain_m;
    /// local domain for the various fields
    NDIndex<Dim> lDomain_m;

    unsigned int ng_m;

    unsigned int kmax_m;

    SimData<T,Dim> simData_m;

    /// fortrans nint function
    inline T nint(T x)
    {
        return ceil(x + 0.5) - (fmod(x*0.5 + 0.25, 1.0) != 0);
    }

    void doInit();

public:

    PwrSpec(ChargedParticles<T,Dim> *beam, SimData<T,Dim> simData);

    ~PwrSpec();

    void calc1D(ChargedParticles<T,Dim> *b);

    void save1DSpectra(string fn);

};
#include "PwrSpec.cc"
#endif

