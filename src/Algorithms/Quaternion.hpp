#ifndef OPAL_QUATERNION_H
#define OPAL_QUATERNION_H

// #include <cstdlib>
// #include <ctime>
#include "Algorithms/BoostMatrix.h"
#include "Ippl.h"

class Quaternion : public ippl::Vector<double, 4> {
public:
    Quaternion();
    Quaternion(const Quaternion&);
    Quaternion(const double&, const double&, const double&, const double&);
    Quaternion(const ippl::Vector<double, 3>&);
    Quaternion(const double&, const ippl::Vector<double, 3>&);
    Quaternion(const matrix_t&);

    Quaternion operator*(const double&) const;
    Quaternion operator*(const Quaternion&) const;
    Quaternion& operator=(const Quaternion&) = default;
    Quaternion& operator*=(const Quaternion&);
    Quaternion operator/(const double&) const;

    double Norm() const;
    double length() const;
    Quaternion& normalize();

    bool isUnit() const;
    bool isPure() const;
    bool isPureUnit() const;

    Quaternion inverse() const;
    Quaternion conjugate() const;

    double real() const;
    ippl::Vector<double, 3> imag() const;

    ippl::Vector<double, 3> rotate(const ippl::Vector<double, 3>&) const;

    matrix_t getRotationMatrix() const;
};

typedef Quaternion Quaternion_t;

Quaternion getQuaternion(ippl::Vector<double, 3> vec, ippl::Vector<double, 3> reference);

inline Quaternion::Quaternion() : Vector<double, 4>(1.0, 0.0, 0.0, 0.0) {
    srand(time(0));
}

inline Quaternion::Quaternion(const Quaternion& quat) : Vector<double, 4>(quat) {
    srand(time(0));
}

inline Quaternion::Quaternion(
    const double& x0, const double& x1, const double& x2, const double& x3)
    : Vector<double, 4>(x0, x1, x2, x3) {
    srand(time(0));
}

inline Quaternion::Quaternion(const ippl::Vector<double, 3>& vec)
    : Quaternion(0.0, vec(0), vec(1), vec(2)) {
    srand(time(0));
}

inline Quaternion::Quaternion(const double& realPart, const ippl::Vector<double, 3>& vec)
    : Quaternion(realPart, vec(0), vec(1), vec(2)) {
    srand(time(0));
}

inline double Quaternion::Norm() const {
    return dot(*this);
}

inline double Quaternion::length() const {
    return std::sqrt(this->Norm());
}

inline bool Quaternion::isUnit() const {
    return (std::abs(this->Norm() - 1.0) < 1e-12);
}

inline bool Quaternion::isPure() const {
    return (std::abs((*this)(0)) < 1e-12);
}

inline bool Quaternion::isPureUnit() const {
    return (this->isPure() && this->isUnit());
}

inline Quaternion Quaternion::conjugate() const {
    Quaternion quat(this->real(), -this->imag());

    return quat;
}

inline double Quaternion::real() const {
    return (*this)(0);
}

inline ippl::Vector<double, 3> Quaternion::imag() const {
    ippl::Vector<double, 3> vec((*this)(1), (*this)(2), (*this)(3));

    return vec;
}

#endif
