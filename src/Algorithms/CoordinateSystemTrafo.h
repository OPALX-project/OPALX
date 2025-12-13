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

    void transformBunchTo(auto Rview);

    void transformBunchFrom(auto Rview);

    void rotateBunchTo(auto Pview);

    void rotateBunchFrom(auto Pview);
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

#endif
