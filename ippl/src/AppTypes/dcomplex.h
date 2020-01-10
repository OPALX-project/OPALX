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
#include <complex>

typedef std::complex<double> dcomplex;

#endif // DCOMPLEX_H
