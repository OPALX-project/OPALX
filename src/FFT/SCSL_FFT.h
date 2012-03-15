// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef IPPL_FFT_SCSL_FFT_H
#define IPPL_FFT_SCSL_FFT_H

// include files
#include "Utility/PAssert.h"
#include "Utility/IpplInfo.h"


/**************************************************************************
 * SCSL_FFT.h:  Prototypes for accessing Fortran 1D FFT routines from the
 * SGI/Cray Scientific Library, and definitions for templated class SCSL,
 * which acts as an FFT engine for the FFT class, providing storage for
 * trigonometric information and performing the 1D FFTs as needed.
 **************************************************************************/


// For platforms that do Fortran symbols in all caps.
#if defined(IPPL_T3E)

#define sgiccfft_ SGICCFFT
#define sgiscfft_ SGISCFFT
#define sgicsfft_ SGICSFFT

#endif

// SCSL Fortran wrapper function prototypes
// These functions are defined in file SCSL_FFT.F, and they just turn
// around and call the native library routines.

#if defined(IPPL_T3E)

// Cray T3E has only "real" routines, which are automatically double precision
extern "C" {
  // double-precision CC FFT
  void sgiccfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
  // double-precision RC FFT
  void sgiscfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
  // double-precision CR FFT
  void sgicsfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
}

#else

// SGI offers separate single- and double-precision routines
extern "C" {
  // double-precision CC FFT
  void sgizzfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
  // double-precision RC FFT
  void sgidzfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
  // double-precision CR FFT
  void sgizdfft_(int& isign, int& n, double& scale, double& x, double& y,
                 double& table, double& work, int& isys);
  // single-precision CC FFT
  void sgiccfft_(int& isign, int& n, float& scale, float& x, float& y,
                 float& table, float& work, int& isys);
  // single-precision RC FFT
  void sgiscfft_(int& isign, int& n, float& scale, float& x, float& y,
                 float& table, float& work, int& isys);
  // single-precision CR FFT
  void sgicsfft_(int& isign, int& n, float& scale, float& x, float& y,
                 float& table, float& work, int& isys);
}

#endif // defined(IPPL_T3E)


// SCSL_wrap provides static functions that wrap the Fortran functions
// in a common interface.  We specialize this class on precision type.
template <class T>
class SCSL_wrap {};

// Specialization for float
template <>
class SCSL_wrap<float> {

public:
  // interface functions used by class SCSL

  // complex-to-complex FFT
  static void ccfft(int isign, int n, float scale, float* x, float* y,
                    float* table, float* work, int isys) {
    sgiccfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // real-to-complex FFT
  static void rcfft(int isign, int n, float scale, float* x, float* y,
                    float* table, float* work, int isys) {
    sgiscfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // complex-to-real FFT
  static void crfft(int isign, int n, float scale, float* x, float* y,
                    float* table, float* work, int isys) {
    sgicsfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }

};


// Specialization for double

#if defined(IPPL_T3E)

// Cray T3E uses single-precision function names for double-precision FFTs
template <>
class SCSL_wrap<double> {

public:
  // interface functions used by class SCSL

  // complex-to-complex FFT
  static void ccfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgiccfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // real-to-complex FFT
  static void rcfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgiscfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // complex-to-real FFT
  static void crfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgicsfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }

};

#else

// SGI has special FFT routines for double-precision
template <>
class SCSL_wrap<double> {

public:
  // interface functions used by class SCSL

  // complex-to-complex FFT
  static void ccfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgizzfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // real-to-complex FFT
  static void rcfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgidzfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }
  // complex-to-real FFT
  static void crfft(int isign, int n, double scale, double* x, double* y,
                    double* table, double* work, int isys) {
    sgizdfft_(isign, n, scale, *x, *y, *table, *work, isys);
  }

};

#endif  // defined(IPPL_T3E)


// Definition of FFT engine class SCSL
template <class T>
class SCSL {

public:

  // definition of complex type
#ifdef IPPL_HAS_TEMPLATED_COMPLEX
  typedef complex<T> Complex_t;
#else
  typedef complex Complex_t;
#endif

  // Trivial constructor.  Do the real work in setup function.
  SCSL(void) {}

  // destructor
  ~SCSL(void);

  // setup internal storage and prepare to do FFTs
  void setup(unsigned numTransformDims, const int* transformTypes,
             const int* axisLengths);

  // invoke FFT on complex data for given dimension and direction
  void callFFT(unsigned transformDim, int direction, Complex_t* data);

  // invoke FFT on real data for given dimension and direction
  void callFFT(unsigned transformDim, int direction, T* data);

private:

  unsigned numTransformDims_m;  // number of dimensions to transform
  int* transformType_m;         // type of transform for each dimension
  int* axisLength_m;            // length along each dimension
  T** table_m;                  // tables of trigonometric factors
  T** work_m;                   // workspace arrays

};


// Inline member function definitions

// destructor
template <class T>
inline
SCSL<T>::~SCSL(void) {
  // delete storage
  for (unsigned d=0; d<numTransformDims_m; ++d) {
    delete [] table_m[d];
    delete [] work_m[d];
  }
  delete [] table_m;
  delete [] work_m;
  delete [] transformType_m;
  delete [] axisLength_m;
}

// setup internal storage and prepare for doing FFTs
template <class T>
inline void
SCSL<T>::setup(unsigned numTransformDims, const int* transformTypes,
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

  // allocate trig table and workspace arrays and initialize table
  // table and work array sizes vary for the Cray T3E and the SGI Origin 2000
  table_m = new T*[numTransformDims_m];
  work_m = new T*[numTransformDims_m];
  T dummy = T();
  for (d=0; d<numTransformDims_m; ++d) {
    switch (transformType_m[d]) {
    case 0:  // CC FFT
#if defined(IPPL_T3E)
      table_m[d] = new T[2 * axisLength_m[d]];
      work_m[d] = new T[4 * axisLength_m[d]];
#else
      table_m[d] = new T[2 * axisLength_m[d] + 30];
      work_m[d] = new T[2 * axisLength_m[d]];
#endif
      // this call initializes table
      SCSL_wrap<T>::ccfft(0, axisLength_m[d], dummy, &dummy, &dummy,
                          table_m[d], &dummy, 0);
      break;
    case 1:  // RC FFT
#if defined(IPPL_T3E)
      table_m[d] = new T[2 * axisLength_m[d]];
      work_m[d] = new T[2 * (axisLength_m[d] + 2)]; // extend length by two
#else
      table_m[d] = new T[axisLength_m[d] + 15];
      work_m[d] = new T[axisLength_m[d] + 2];  // extend length by two reals
#endif
      // this call initializes table
      SCSL_wrap<T>::rcfft(0, axisLength_m[d], dummy, &dummy, &dummy,
                          table_m[d], &dummy, 0);
      break;
    case 2:  // Sine transform
      ERRORMSG("No sine transforms available in SCSL!!" << endl);
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
SCSL<T>::callFFT(unsigned transformDim, int direction,
                 SCSL<T>::Complex_t* data) {

  // check transform dimension and direction arguments
  PAssert(transformDim<numTransformDims_m);
  PAssert(direction==+1 || direction==-1);

  // cast complex number pointer to T* for calling Fortran routines
  T* rdata = reinterpret_cast<T*>(data);

  // branch on type of transform for this dimension
  switch (transformType_m[transformDim]) {
  case 0:  // CC FFT
    SCSL_wrap<T>::ccfft(direction, axisLength_m[transformDim], 1.0, rdata,
                        rdata, table_m[transformDim], work_m[transformDim], 0);
    break;
  case 1:  // RC FFT
    if (direction == +1) {  // real-to-complex
      SCSL_wrap<T>::rcfft(+1, axisLength_m[transformDim], 1.0, rdata, rdata,
                          table_m[transformDim], work_m[transformDim], 0);
    }
    else {                  // complex-to-real
      SCSL_wrap<T>::crfft(-1, axisLength_m[transformDim], 1.0, rdata, rdata,
                          table_m[transformDim], work_m[transformDim], 0);
    }
    break;
  case 2:  // Sine transform
    ERRORMSG("No sine transforms available in SCSL!!" << endl);
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
SCSL<T>::callFFT(unsigned transformDim, int direction, T* data) {

  // check transform dimension argument
  PAssert(transformDim<numTransformDims_m);
  // branch on transform type for this dimension
  switch (transformType_m[transformDim]) {
  case 0:  // CC FFT
    ERRORMSG("complex-to-complex FFT uses complex data!!" << endl);
    break;
  case 1:  // RC FFT
    ERRORMSG("real-to-complex FFT uses complex data!!" << endl);
    break;
  case 2:  // Sine transform
    ERRORMSG("No sine transforms available in SCSL!!" << endl);
    break;
  default:
    ERRORMSG("Unknown transform type requested!!" << endl);
    break;
  }

  return;
}


#endif // IPPL_FFT_SCSL_FFT_H

/***************************************************************************
 * $RCSfile: SCSL_FFT.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:25 $
 * IPPL_VERSION_ID: $Id: SCSL_FFT.h,v 1.1.1.1 2003/01/23 07:40:25 adelmann Exp $ 
 ***************************************************************************/

