#include "Algorithms/Vektor.h"
#include "AppTypes/Tenzor.h"

void normalize(Vector_t & vec)
{
    double length = sqrt(dot(vec, vec));

    PAssert(length > 1e-12);

    vec /= length;
}