#ifndef OPAL_VEKTOR_HH
#define OPAL_VEKTOR_HH

#include <limits>

#include "AppTypes/Vektor.h"

template <class, unsigned>
class Tenzor;

typedef Vektor<double, 3> Vector_t;

void normalize(Vector_t & vec);

#include "Algorithms/Quaternion.h"

#endif