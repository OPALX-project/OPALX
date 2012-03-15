#ifndef _RANLIB_H_
#define _RANLIB_H_
//=================================================================
// class library for random number generation
// ------------------------------------------
// Random number generation is based on a C++ port of the RANLUX
// implementation by F. James with an improved padding scheme,
// yielding full 24-bit or 48-bit precision in the mantissa of
// random floats or doubles, respectively.
//
// Author: Michael Schmelling / May-1999
//=================================================================

#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <string.h>

class RANLIB_class {
private:

    // variables defining the internal state of the generator

    int    nwork;    // size of work space for ilux-kernel
    int    *work;    // work space for the ilux-kernel
    int    seed;     // seed used to initialize the generator
    int    luxl;     // lux-level of the ilux-kernel
    int    rpad;     // rndm-integer with padding bits
    int    npad;     // number of available bitsin rpad
    int    major;    // call-count to the kernel major*10^9+minor
    int    minor;    // call-count to the kernel major*10^9+minor

    // auxiliary data for bit padding

    int     maxpad;  // maximum number of padding bits
    float  *scal_f;  // conversion factors int -> float
    double *scal_d;  // conversion factors int -> double

    // private functions

    void   init();
    void   iluxinit(int seed, int luxlevel, int *work);
    int    ilux(int *work);

public:
    // destructor

    ~RANLIB_class();

    // constructors

    RANLIB_class();                    // constructor with defaults
    RANLIB_class(int seed);            // specify random seed
    RANLIB_class(int seed, int luxl);  // specify seed and luxlevel

    // access to private data

    int    val_seed()     {return seed ;};   // initial seed
    int    val_luxl()     {return luxl ;};   // luxlevel

    // return the number of calls to the ilux()-kernel

    double kernel_calls() {return (major * 1000000000. + minor);};

    // other member functions

    int    save(const char *filename);     // save generator status to a file
    int    restore(const char *filename);  // restore the status from a file
    int    selftest();                     // perform a selftest of the generator

    int    ilux();   // RANLUX integer generator [0,16777215]
    float  flux();   // basic float generator  for the interval ]0,1[
    double dlux();   // basic double generator for the interval ]0,1[


    //=============================================================================
    // extra member functions for RANLIB_class
    // ---------------------------------------
    // This header file defines the prototypes of random number generators
    // not defined as part of the basic functionality of the RANLIB_class,
    // which only provides random deviates in the range ]0,1[. It is activated
    // by #define RANLIBX before including ranlib.h.
    //
    // Author: Michael Schmelling / May-1999
    //=============================================================================

private:

    // auxiliary function to nbd()

    double gammln(double x);

public:

    // uniform generators for the interval ]low,high[

    float  uniform(float  low, float  high);
    double uniform(double low, double high);

    // gaussian random numbers with mean av and standard deviation sigma

    float  gauss(float  av, float  sigma);
    double gauss(double av, double sigma);

    // negative binomial generator

    int nbd(double av, double k);

    // power law generator x/(1+x)^a

    float  rnlat(float a);
    double rnlat(double a);

    // power law with a kink

    double lgPowerK(double lgxl, double lgxh, double lgxk, double gl, double gh);

};
#endif









