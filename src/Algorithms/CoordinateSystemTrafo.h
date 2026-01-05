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

    matrix_t getRotationMatrix() const;

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


inline matrix_t CoordinateSystemTrafo::getRotationMatrix() const {
    return rotationMatrix_m;
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
    // Manual computation to verify:
    //ippl::Vector<double, 3> manual;
    //manual[0] = -1.0 * r[0] + 0.0 * r[1] + 0.0 * r[2];
    //manual[1] = 0.0 * r[0] -1.0 * r[1] + 0.0 * r[2];  
    //manual[2] = 0.0 * r[0] + 0.0 * r[1] + 1.0 * r[2];
    //std::cout << "Manual: " << manual << endl;
    //std::cout << "Computed: " << prod_boost_vector(rotationMatrix_m, r) << endl;
    return prod_boost_vector(rotationMatrix_m, r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateFrom(
    const ippl::Vector<double, 3>& r) const {
    //std::cout << "rotateFrom called" << std::endl;
    //std::cout << "rotationMatrix_m: " << getRotationMatrix() << std::endl;
    //std::cout << "r: " << r << std::endl;
    // matrix_t tmpMat = boost::numeric::ublas::trans(rotationMatrix_m);
    //std::cout << "transposed rotationMatrix_m: " << tmpMat << std::endl;
    /*
    std::cout << "transposed rotationMatrix_m: [";
    for (size_t i = 0; i < tmpMat.size1(); ++i) {
        for (size_t j = 0; j < tmpMat.size2(); ++j) {
            std::cout << tmpMat(i, j);
            if (j < tmpMat.size2() - 1) std::cout << " ";
        }
        if (i < tmpMat.size1() - 1) std::cout << "; ";
    }
    std::cout << "]" << std::endl; */

    // Also explicitly print out rotationMatrix_m
    /*std::cout << "rotationMatrix_m: [";
    for (size_t i = 0; i < rotationMatrix_m.size1(); ++i) {
        for (size_t j = 0; j < rotationMatrix_m.size2(); ++j) {
            std::cout << rotationMatrix_m(i, j);
            if (j < rotationMatrix_m.size2() - 1) std::cout << " ";
        }
        if (i < rotationMatrix_m.size1() - 1) std::cout << "; ";
    }
    std::cout << "]" << std::endl;*/
    
    //std::cout << "result: " << prod_boost_vector(boost::numeric::ublas::trans(rotationMatrix_m), r) << std::endl;
    
    return prod_boost_vector(boost::numeric::ublas::trans(rotationMatrix_m), r);
}

#endif
