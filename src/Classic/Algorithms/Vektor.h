#ifndef OPAL_VEKTOR_HH
#define OPAL_VEKTOR_HH

#include <limits>

#include "AppTypes/Vektor.h"

typedef Vektor<double, 3> Vector_t;

inline
double euclidian_norm(Vector_t a) {
    return sqrt(dot(a,a));
}

#endif