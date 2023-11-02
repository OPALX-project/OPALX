#ifndef COORDINATESYSTEMTRAFO
#define COORDINATESYSTEMTRAFO

#include "Algorithms/Quaternion.hpp"

class CoordinateSystemTrafo {
public:
    CoordinateSystemTrafo();

    CoordinateSystemTrafo(const CoordinateSystemTrafo& right);

    CoordinateSystemTrafo(const ippl::Vector<double, 3>& origin, const Quaternion& orientation);

    ippl::Vector<double, 3> transformTo(const ippl::Vector<double, 3>& r) const;
    ippl::Vector<double, 3> transformFrom(const ippl::Vector<double, 3>& r) const;

    ippl::Vector<double, 3> rotateTo(const ippl::Vector<double, 3>& r) const;
    ippl::Vector<double, 3> rotateFrom(const ippl::Vector<double, 3>& r) const;

    void invert();
    CoordinateSystemTrafo inverted() const;

    CoordinateSystemTrafo& operator=(const CoordinateSystemTrafo& right) = default;
    CoordinateSystemTrafo operator*(const CoordinateSystemTrafo& right) const;
    void operator*=(const CoordinateSystemTrafo& right);

    ippl::Vector<double, 3> getOrigin() const;
    Quaternion getRotation() const;

    void print(std::ostream&) const;

private:
    ippl::Vector<double, 3> origin_m;
    Quaternion orientation_m;
    matrix_t rotationMatrix_m;
};

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

inline ippl::Vector<double, 3> CoordinateSystemTrafo::getOrigin() const {
    return origin_m;
}

inline Quaternion CoordinateSystemTrafo::getRotation() const {
    return orientation_m;
}

inline CoordinateSystemTrafo CoordinateSystemTrafo::inverted() const {
    CoordinateSystemTrafo result(*this);
    result.invert();

    return result;
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::transformTo(
    const ippl::Vector<double, 3>& r) const {
    const ippl::Vector<double, 3> delta = r - origin_m;
    return prod_boost_vector(rotationMatrix_m, delta);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::transformFrom(
    const ippl::Vector<double, 3>& r) const {
    return rotateFrom(r) + origin_m;
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateTo(
    const ippl::Vector<double, 3>& r) const {
    return prod_boost_vector(rotationMatrix_m, r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateFrom(
    const ippl::Vector<double, 3>& r) const {
    return prod_boost_vector(boost::numeric::ublas::trans(rotationMatrix_m), r);
}

#endif
