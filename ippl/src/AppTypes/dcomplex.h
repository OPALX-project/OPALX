// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef DCOMPLEX_H
#define DCOMPLEX_H

/***********************************************************************
 * 
 * Work around the lack of draft standard complex<T> in all compilers.
 * Correctly declare a dcomplex typedef based on the compiler capabilities
 * and available C++ standard library.  dcomplex is a complex number class
 * storing values as doubles.
 *
 ***********************************************************************/

// include standard complex header file
#ifdef IPPL_USE_STANDARD_HEADERS
#include <complex>
#else
#include <complex.h>
#endif


#ifdef IPPL_HAS_TEMPLATED_COMPLEX

// KAI and others have a templated complex class
#ifdef IPPL_USE_SINGLE_PRECISION
typedef complex<float> dcomplex;
#else // USE_DOUBLE_PRECISION
typedef std::complex<double> dcomplex;
#endif

typedef std::complex<float> fComplex;

#else

// This assumes that all other compilers have the old non-templated
// complex type which is like complex<double> in the draft standard.
typedef complex dcomplex;

#endif // IPPL_HAS_TEMPLATED_COMPLEX

#endif // DCOMPLEX_H

/***************************************************************************
 * $RCSfile: dcomplex.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:24 $
 * IPPL_VERSION_ID: $Id: dcomplex.h,v 1.1.1.1 2003/01/23 07:40:24 adelmann Exp $ 
 ***************************************************************************/
