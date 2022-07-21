
#ifndef PYOPAL_PyPolynomialCoefficient_H
#define PYOPAL_PyPolynomialCoefficient_H

#include <Python.h>


namespace interpolation {

// note following are in interpolation namespace
class PolynomialCoefficient;
}

namespace PyPolynomialCoefficient {

/** PyPolynomialMap is the python implementation of the C++ PolynomialMap class
 *
 *  Provides a multivariate polynomial object
 */
typedef struct {
    PyObject_HEAD;
    interpolation::PolynomialCoefficient* coeff;
} PyCoefficient;

}

#endif