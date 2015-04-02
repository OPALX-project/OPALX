// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef IPPL_FFT_FFTPACK_FFT_H
#define IPPL_FFT_FFTPACK_FFT_H

// include files
#include "Utility/PAssert.h"
#include "Utility/IpplInfo.h"


/**************************************************************************
 * fftpack_FFT.h:  Prototypes for accessing Fortran 1D FFT routines from
 * Netlib, and definitions for templated class FFTPACK, which acts as an
 * FFT engine for the FFT class, providing storage for trigonometric
 * information and performing the 1D FFTs as needed.
 **************************************************************************/

// For platforms that do Fortran symbols in all caps.
#if defined(IPPL_T3E)

#define cffti_ CFFTI
#define cfftf_ CFFTF
#define cfftb_ CFFTB
#define rffti_ RFFTI
#define rfftf_ RFFTF
#define rfftb_ RFFTB
#define sinti_ SINTI
#define sint_ SINT
#define fcffti_ FCFFTI
#define fcfftf_ FCFFTF
#define fcfftb_ FCFFTB
#define frffti_ FRFFTI
#define frfftf_ FRFFTF
#define frfftb_ FRFFTB
#define fsinti_ FSINTI
#define fsint_ FSINT

#endif

// For platforms that do Fortran symbols just like C symbols.
#if defined(IPPL_SP2)

#define cffti_ cffti
#define cfftf_ cfftf
#define cfftb_ cfftb
#define rffti_ rffti
#define rfftf_ rfftf
#define rfftb_ rfftb
#define sinti_ sinti
#define sint_ sint
#define fcffti_ fcffti
#define fcfftf_ fcfftf
#define fcfftb_ fcfftb
#define frffti_ frffti
#define frfftf_ frfftf
#define frfftb_ frfftb
#define fsinti_ fsinti
#define fsint_ fsint

#endif

// FFTPACK function prototypes for Fortran routines
extern "C" {
  // double-precision CC FFT
  void cffti_(int& n, double& wsave);
  void cfftf_(int& n, double& r, double& wsave);
  void cfftb_(int& n, double& r, double& wsave);
  // double-precision RC FFT
  void rffti_(int& n, double& wsave);
  void rfftf_(int& n, double& r, double& wsave);
  void rfftb_(int& n, double& r, double& wsave);
  // double-precision sine transform
  void sinti_(int& n, double& wsave);
  void sint_(int& n, double& r, double& wsave);
  // single-precision CC FFT
  void fcffti_(int& n, float& wsave);
  void fcfftf_(int& n, float& r, float& wsave);
  void fcfftb_(int& n, float& r, float& wsave);
  // single-precision RC FFT
  void frffti_(int& n, float& wsave);
  void frfftf_(int& n, float& r, float& wsave);
  void frfftb_(int& n, float& r, float& wsave);
  // single-precision sine transform
  void fsinti_(int& n, float& wsave);
  void fsint_(int& n, float& r, float& wsave);
}


// FFTPACK_wrap provides static functions that wrap the Fortran functions
// in a common interface.  We specialize this class on precision type.
template <class T>
class FFTPACK_wrap {};

// Specialization for float
template <>
class FFTPACK_wrap<float> {

public:
  // interface functions used by class FFTPACK

  // initialization functions for CC FFT, RC FFT, and sine transform
  static void ccffti(int n, float* wsave) { fcffti_(n, *wsave); }
  static void rcffti(int n, float* wsave) { frffti_(n, *wsave); }
  static void rrffti(int n, float* wsave) { fsinti_(n, *wsave); }
  // forward and backward CC FFT
  static void ccfftf(int n, float* r, float* wsave) { fcfftf_(n, *r, *wsave); }
  static void ccfftb(int n, float* r, float* wsave) { fcfftb_(n, *r, *wsave); }
  // forward and backward RC FFT
  static void rcfftf(int n, float* r, float* wsave) { frfftf_(n, *r, *wsave); }
  static void rcfftb(int n, float* r, float* wsave) { frfftb_(n, *r, *wsave); }
  // sine transform
  static void rrfft(int n, float* r, float* wsave) { fsint_(n, *r, *wsave); }

};

// Specialization for double
template <>
class FFTPACK_wrap<double> {

public:
  // interface functions used by class FFTPACK

  // initialization functions for CC FFT, RC FFT, and sine transform
  static void ccffti(int n, double* wsave) { cffti_(n, *wsave); }
  static void rcffti(int n, double* wsave) { rffti_(n, *wsave); }
  static void rrffti(int n, double* wsave) { sinti_(n, *wsave); }
  // forward and backward CC FFT
  static void ccfftf(int n, double* r, double* wsave) {cfftf_(n, *r, *wsave);}
  static void ccfftb(int n, double* r, double* wsave) {cfftb_(n, *r, *wsave);}
  // forward and backward RC FFT
  static void rcfftf(int n, double* r, double* wsave) {rfftf_(n, *r, *wsave);}
  static void rcfftb(int n, double* r, double* wsave) {rfftb_(n, *r, *wsave);}
  // sine transform
  static void rrfft(int n, double* r, double* wsave) { sint_(n, *r, *wsave); }

};


// Definition of FFT engine class FFTPACK
template <class T>
class FFTPACK {

public:

  // definition of complex type
#ifdef IPPL_HAS_TEMPLATED_COMPLEX
  typedef std::complex<T> Complex_t;
#else
  typedef complex Complex_t;
#endif

  // Trivial constructor.  Do the real work in setup function.
  FFTPACK(void) {}

  // destructor
  ~FFTPACK(void);

  // setup internal storage and prepare to perform FFTs
  // inputs are number of dimensions to transform, the transform types,
  // and the lengths of these dimensions.
  void setup(unsigned numTransformDims, const int* transformTypes,
             const int* axisLengths);

  // invoke FFT on complex data for given dimension and direction
  void callFFT(unsigned transformDim, int direction, Complex_t* data);

  // invoke FFT on real data for given dimension and direction
  void callFFT(unsigned transformDim, int direction, T* data);

private:

  unsigned numTransformDims_m;  // number of dimensions to transform
  int* transformType_m;         // transform type for each dimension
  int* axisLength_m;            // length of each transform dimension
  T** trig_m;                   // trigonometric tables

};


// Inline member function definitions

// destructor
template <class T>
inline
FFTPACK<T>::~FFTPACK(void) {
  // delete storage
  for (unsigned d=0; d<numTransformDims_m; ++d)
    delete [] trig_m[d];
  delete [] trig_m;
  delete [] transformType_m;
  delete [] axisLength_m;
}

// setup internal storage to prepare for FFTs
template <class T>
inline void
FFTPACK<T>::setup(unsigned numTransformDims, const int* transformTypes,
                  const int* axisLengths) {

  // store transform types and lengths for each transform dim
  numTransformDims_m = numTransformDims;
  transformType_m = new int[numTransformDims_m];
  axisLength_m = new int[numTransformDims_m];
  unsigned d;
  for (d=0; d<numTransformDims_m; ++d) {
    transformType_m[d] = transformTypes[d];
    axisLength_m[d] = axisLengths[d];
  }

  // allocate and initialize trig table
  trig_m = new T*[numTransformDims_m];
  for (d=0; d<numTransformDims_m; ++d) {
    switch (transformType_m[d]) {
    case 0:  // CC FFT
      trig_m[d] = new T[4 * axisLength_m[d] + 15];
      FFTPACK_wrap<T>::ccffti(axisLength_m[d], trig_m[d]);
      break;
    case 1:  // RC FFT
      trig_m[d] = new T[2 * axisLength_m[d] + 15];
      FFTPACK_wrap<T>::rcffti(axisLength_m[d], trig_m[d]);
      break;
    case 2:  // Sine transform
      trig_m[d] = new T[static_cast<int>(2.5 * axisLength_m[d] + 0.5) + 15];
      FFTPACK_wrap<T>::rrffti(axisLength_m[d], trig_m[d]);
      break;
    default:
      ERRORMSG("Unknown transform type requested!!" << endl);
      break;
    }
  }

  return;
}

// invoke FFT on complex data for given dimension and direction
template <class T>
inline void
FFTPACK<T>::callFFT(unsigned transformDim, int direction,
                    FFTPACK<T>::Complex_t* data) {

  // check transform dimension and direction arguments
  PAssert(transformDim<numTransformDims_m);
  PAssert(direction==+1 || direction==-1);

  // cast complex number pointer to T* for calling Fortran routines
  T* rdata = reinterpret_cast<T*>(data);

  // branch on transform type for this dimension
  switch (transformType_m[transformDim]) {
  case 0:  // CC FFT
    if (direction == +1) {
      // call forward complex-to-complex FFT
      FFTPACK_wrap<T>::ccfftf(axisLength_m[transformDim], rdata,
                              trig_m[transformDim]);
    }
    else {
      // call backward complex-to-complex FFT
      FFTPACK_wrap<T>::ccfftb(axisLength_m[transformDim], rdata,
                              trig_m[transformDim]);
    }
    break;
  case 1:  // RC FFT
    if (direction == +1) {
      // call forward real-to-complex FFT
      FFTPACK_wrap<T>::rcfftf(axisLength_m[transformDim], rdata,
                              trig_m[transformDim]);
      // rearrange output to conform with SGI format for complex result
      int clen = axisLength_m[transformDim]/2 + 1;
      data[clen-1] = Complex_t(imag(data[clen-2]),0.0);
      for (int i = clen-2; i > 0; --i)
        data[i] = Complex_t(imag(data[i-1]),real(data[i]));
      data[0] = Complex_t(real(data[0]),0.0);
    }
    else {                
      // rearrange input to conform with Netlib format for complex modes
      int clen = axisLength_m[transformDim]/2 + 1;
      data[0] = Complex_t(real(data[0]),real(data[1]));
      for (int i = 1; i < clen-1; ++i)
        data[i] = Complex_t(imag(data[i]),real(data[i+1]));
      // call backward complex-to-real FFT
      FFTPACK_wrap<T>::rcfftb(axisLength_m[transformDim], rdata,
                              trig_m[transformDim]);
    }
    break;
  case 2:  // Sine transform
    ERRORMSG("Input for real-to-real FFT should be real!!" << endl);
    break;
  default:
    ERRORMSG("Unknown transform type requested!!" << endl);
    break;
  }

  return;
}

// invoke FFT on real data for given dimension and direction
template <class T>
inline void
FFTPACK<T>::callFFT(unsigned transformDim, int direction, T* data) {

  // check transform dimension and direction arguments
  PAssert(transformDim<numTransformDims_m);
  PAssert(direction==+1 || direction==-1);

  // branch on transform type for this dimension
  switch (transformType_m[transformDim]) {
  case 0:  // CC FFT
    ERRORMSG("Input for complex-to-complex FFT should be complex!!" << endl);
    break;
  case 1:  // RC FFT
    ERRORMSG("real-to-complex FFT uses complex input!!" << endl);
    break;
  case 2:  // Sine transform
    // invoke the real-to-real transform on the data
    FFTPACK_wrap<T>::rrfft(axisLength_m[transformDim], data,
                           trig_m[transformDim]);
    break;
  default:
    ERRORMSG("Unknown transform type requested!!" << endl);
    break;
  }

  return;
}


#endif // IPPL_FFT_FFTPACK_FFT_H

/***************************************************************************
 * $RCSfile: fftpack_FFT.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
 * IPPL_VERSION_ID: $Id: fftpack_FFT.h,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
 ***************************************************************************/

