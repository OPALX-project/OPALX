#include "Algorithms/CoordinateSystemTrafo.h"

CoordinateSystemTrafo::CoordinateSystemTrafo()
    : origin_m(0.0), orientation_m(1.0, 0.0, 0.0, 0.0), rotationMatrix_m() {
    rotationMatrix_m(0, 0) = 1.0;
    rotationMatrix_m(0, 1) = 0.0;
    rotationMatrix_m(0, 2) = 0.0;
    rotationMatrix_m(1, 0) = 0.0;
    rotationMatrix_m(1, 1) = 1.0;
    rotationMatrix_m(1, 2) = 0.0;
    rotationMatrix_m(2, 0) = 0.0;
    rotationMatrix_m(2, 1) = 0.0;
    rotationMatrix_m(2, 2) = 1.0;
}

CoordinateSystemTrafo::CoordinateSystemTrafo(const CoordinateSystemTrafo& right)
    : origin_m(right.origin_m),
      orientation_m(right.orientation_m),
      rotationMatrix_m(right.rotationMatrix_m) {
}

CoordinateSystemTrafo::CoordinateSystemTrafo(
    const ippl::Vector<double, 3>& origin, const Quaternion& orientation)
    : origin_m(origin),
      orientation_m(orientation),
      rotationMatrix_m(orientation_m.getRotationMatrix()) {
}

void CoordinateSystemTrafo::invert() {
    origin_m         = -orientation_m.rotate(origin_m);
    orientation_m    = orientation_m.conjugate();
    rotationMatrix_m = get_transpose(rotationMatrix_m);
}

CoordinateSystemTrafo CoordinateSystemTrafo::operator*(const CoordinateSystemTrafo& right) const {
    CoordinateSystemTrafo result(*this);

    result *= right;
    return result;
}

void CoordinateSystemTrafo::operator*=(const CoordinateSystemTrafo& right) {
    origin_m = right.orientation_m.conjugate().rotate(origin_m) + right.origin_m;
    orientation_m *= right.orientation_m;
    orientation_m.normalize();
    rotationMatrix_m = orientation_m.getRotationMatrix();
}

/* ====================== Transformation Functions ========================== */


void CoordinateSystemTrafo::transformBunchTo(auto Rview)
{
    // local copies for the kernel 
    auto rotationMatrix = rotationMatrix_m;
    auto origin = origin_m;
    Kokkos::parallel_for("transformBunchTo()", ippl::getRangePolicy(Rview),
    KOKKOS_LAMBDA(const int i)
    {
        const ippl::Vector<double, 3> delta = Rview(i) - origin;
        Rview(i) = prod_vector(rotationMatrix,delta);
    });
}

void CoordinateSystemTrafo::transformBunchFrom(auto Rview)
{
     // local copies for the kernel 
    auto rotationMatrix = rotationMatrix_m;
    auto origin = origin_m;
    Kokkos::parallel_for("transformBunchFrom()", ippl::getRangePolicy(Rview),
    KOKKOS_LAMBDA(const int i)
    {
        Rview(i) = prod_vector_transpose(rotationMatrix, Rview(i)) + origin;
    });
}

void CoordinateSystemTrafo::rotateBunchTo(auto Pview)
{
    // local copy for the kernel
    auto rotationMatrix = rotationMatrix_m;
    Kokkos::parallel_for("rotateBunchTo()", ippl::getRangePolicy(Pview),
    KOKKOS_LAMBDA(const int i)
    {
        Pview(i) = prod_vector(rotationMatrix, Pview(i));
    });
} 

void CoordinateSystemTrafo::rotateBunchFrom(auto Pview)
{
    // local copy for the kernel
    auto rotationMatrix = rotationMatrix_m;
    Kokkos::parallel_for("rotateBunchFrom()", ippl::getRangePolicy(Pview),
    KOKKOS_LAMBDA(const int i)
    {
        Pview(i) = prod_vector_transpose(rotationMatrix, Pview(i));
    });
}

/* ========================================================================== */
/* ============================ Getters ===================================== */

/* ========================================================================== */