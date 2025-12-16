// This Class implements transformations between coordinate systems

#ifndef COORDINATESYSTEMTRAFO
#define COORDINATESYSTEMTRAFO

#include "Algorithms/Quaternion.hpp"

class CoordinateSystemTrafo {
public:
/* ============================== Constructors ============================== */
    CoordinateSystemTrafo();

    CoordinateSystemTrafo(
        const CoordinateSystemTrafo& right);

    CoordinateSystemTrafo(
        const ippl::Vector<double, 3>& origin, 
        const Quaternion& orientation);

    void invert();
    CoordinateSystemTrafo inverted() const;
/* ========================================================================== */    
/* ======================== Transformation Functions ======================== */
    ippl::Vector<double, 3> transformTo(
        const ippl::Vector<double, 3>& r) const;

    ippl::Vector<double, 3> transformFrom(
        const ippl::Vector<double, 3>& r) const;

    ippl::Vector<double, 3> rotateTo(
        const ippl::Vector<double, 3>& r) const;

    ippl::Vector<double, 3> rotateFrom(
        const ippl::Vector<double, 3>& r) const;

    void transformBunchTo(auto Rview) const;

    void transformBunchFrom(auto Rview) const;

    void rotateBunchTo(auto Pview) const;

    void rotateBunchFrom(auto Pview) const;
/* ========================================================================== */
/* =============================== Operators ================================ */
    CoordinateSystemTrafo& operator=(
        const CoordinateSystemTrafo& right) = default;

    CoordinateSystemTrafo operator*(
        const CoordinateSystemTrafo& right) const;

    void operator*=(
        const CoordinateSystemTrafo& right);
/* ========================================================================== */
/* =============================== Getters ================================== */
    ippl::Vector<double, 3> getOrigin() const;
    Quaternion getRotation() const;
    Matrix_t getRotationMatrix() const;
/* ========================================================================== */
/* =============================== Getters ================================== */
    void print(std::ostream&) const;

private:
    // Displacement
    ippl::Vector<double, 3> origin_m;

    //Rotation
    Quaternion orientation_m;
    Matrix_t rotationMatrix_m;
};

/* ======== Print ======== */
inline std::ostream& operator<<(std::ostream& os, const CoordinateSystemTrafo& trafo) {
    trafo.print(os);
    return os;
}

inline Inform& operator<<(Inform& os, const CoordinateSystemTrafo& trafo) {
    trafo.print(os.getStream());
    return os;
}

inline void CoordinateSystemTrafo::print(std::ostream& os) const {
    os << "Origin: " << origin_m << "\n"
       << "z-axis: " << orientation_m.conjugate().rotate(ippl::Vector<double, 3>(0, 0, 1)) << "\n"
       << "x-axis: " << orientation_m.conjugate().rotate(ippl::Vector<double, 3>(1, 0, 0));
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::transformTo(
    const ippl::Vector<double, 3>& r) const {
    const ippl::Vector<double, 3> delta = r - origin_m;
    return prod_vector(rotationMatrix_m,delta);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::transformFrom(
    const ippl::Vector<double, 3>& r) const {
    return rotateFrom(r) + origin_m;
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateTo(
    const ippl::Vector<double, 3>& r) const {
    return prod_vector(rotationMatrix_m,r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateFrom(
    const ippl::Vector<double, 3>& r) const {
    return prod_vector_transpose(rotationMatrix_m,r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::getOrigin() const {
    return origin_m;
}

inline Quaternion CoordinateSystemTrafo::getRotation() const {
    return orientation_m;
}


inline Matrix_t CoordinateSystemTrafo::getRotationMatrix() const {
    return rotationMatrix_m;
}

inline CoordinateSystemTrafo CoordinateSystemTrafo::inverted() const {
    CoordinateSystemTrafo result(*this);
    result.invert();

    return result;
}

inline void CoordinateSystemTrafo::transformBunchTo(auto Rview) const{
    Matrix_t rot = rotationMatrix_m;
    ippl::Vector<double, 3> ori = origin_m;
    Kokkos::parallel_for("transformBunchTo", Rview.extent(0), 
        KOKKOS_LAMBDA(const size_t i) {
            ippl::Vector<double, 3> delta = Rview(i) - ori;
            Rview(i) = prod_vector(rot, delta);
        });
}

inline void CoordinateSystemTrafo::transformBunchFrom(auto Rview) const{
    Matrix_t rot = rotationMatrix_m;
    ippl::Vector<double, 3> ori = origin_m;
    Kokkos::parallel_for("transformBunchFrom", Rview.extent(0), 
        KOKKOS_LAMBDA(const size_t i) {
            Rview(i) = prod_vector_transpose(rot, Rview(i)) + ori;
        });
}

inline void CoordinateSystemTrafo::rotateBunchTo(auto Pview) const{
    Matrix_t rot = rotationMatrix_m;
    Kokkos::parallel_for("rotateBunchTo", Pview.extent(0), 
        KOKKOS_LAMBDA(const size_t i) {
            Pview(i) = prod_vector(rot, Pview(i));
        });
}

inline void CoordinateSystemTrafo::rotateBunchFrom(auto Pview) const{
    Matrix_t rot = rotationMatrix_m;
    Kokkos::parallel_for("rotateBunchFrom", Pview.extent(0), 
        KOKKOS_LAMBDA(const size_t i) {
            Pview(i) = prod_vector_transpose(rot, Pview(i));
        });
}

#endif
