#include "Algorithms/Quaternion.hpp"
#include <cmath>
#include <iostream>
#include "Algorithms/BoostMatrix.h"
#include "OPALTypes.h"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

namespace {
    ippl::Vector<double, 3> normalize(const ippl::Vector<double, 3>& vec) {
        double length = std::sqrt(dot(vec, vec));

#ifndef NOPAssert
        if (length < 1e-12)
            throw GeneralClassicException("normalize()", "length of vector less than 1e-12");
#endif
        length = (length == 0) ? 1 : length; // just to avoid division by zero and nan results
        return vec / length;
    }
}  // namespace

Quaternion::Quaternion(const matrix_t& M) : ippl::Vector<double, 4>(0.0) {
    (*this)(0) = std::sqrt(std::max(0.0, 1 + M(0, 0) + M(1, 1) + M(2, 2))) / 2;
    (*this)(1) = std::sqrt(std::max(0.0, 1 + M(0, 0) - M(1, 1) - M(2, 2))) / 2;
    (*this)(2) = std::sqrt(std::max(0.0, 1 - M(0, 0) + M(1, 1) - M(2, 2))) / 2;
    (*this)(3) = std::sqrt(std::max(0.0, 1 - M(0, 0) - M(1, 1) + M(2, 2))) / 2;
    (*this)(1) =
        std::abs(M(2, 1) - M(1, 2)) > 0 ? std::copysign((*this)(1), M(2, 1) - M(1, 2)) : 0.0;
    (*this)(2) =
        std::abs(M(0, 2) - M(2, 0)) > 0 ? std::copysign((*this)(2), M(0, 2) - M(2, 0)) : 0.0;
    (*this)(3) =
        std::abs(M(1, 0) - M(0, 1)) > 0 ? std::copysign((*this)(3), M(1, 0) - M(0, 1)) : 0.0;
}

Quaternion getQuaternion(ippl::Vector<double, 3> u, ippl::Vector<double, 3> ref) {
    const double tol = 1e-12;

    u   = normalize(u);
    ref = normalize(ref);

    ippl::Vector<double, 3> axis = cross(u, ref);
    double normAxis              = std::sqrt(dot(axis, axis));

    if (normAxis < tol) {
        if (std::abs(dot(u, ref) - 1.0) < tol) {
            return Quaternion(1.0, ippl::Vector<double, 3>(0.0));
        }
        // vectors are parallel or antiparallel
        do {  // find any vector in plane with ref as normal
            double u = rand() / RAND_MAX;
            double v = 2 * Physics::pi * rand() / RAND_MAX;
            axis(0)  = std::sqrt(1 - u * u) * std::cos(v);
            axis(1)  = std::sqrt(1 - u * u) * std::sin(v);
            axis(2)  = u;
        } while (std::abs(dot(axis, ref)) > 0.9);

        axis -= dot(axis, ref) * ref;
        axis = normalize(axis);

        return Quaternion(0, axis);
    }

    axis /= normAxis;

    double cosAngle = sqrt(0.5 * (1 + dot(u, ref)));
    double sinAngle = sqrt(1 - cosAngle * cosAngle);

    return Quaternion(cosAngle, sinAngle * axis);
}

Quaternion Quaternion::operator*(const double& d) const {
    Quaternion result(d * this->real(), d * this->imag());

    return result;
}

Quaternion Quaternion::operator*(const Quaternion& other) const {
    Quaternion result(*this);
    return result *= other;
}

Quaternion& Quaternion::operator*=(const Quaternion& other) {
    ippl::Vector<double, 3> imagThis  = this->imag();
    ippl::Vector<double, 3> imagOther = other.imag();

    /// \todo there is a scope issue w.r.t dot and ippl::Vector
    double res = 0.0;
    for (unsigned i = 0; i < 3; i++)
        res += imagThis(i) * imagOther(i);
    double dp = std::sqrt(res);

    *this = Quaternion(
        (*this)(0) * other(0) - dp,
        (*this)(0) * imagOther + other(0) * imagThis + cross(imagThis, imagOther));
    return *this;
}

Quaternion Quaternion::operator/(const double& d) const {
    Quaternion result(this->real() / d, this->imag() / d);

    return result;
}

Quaternion& Quaternion::normalize() {
#ifndef NOPAssert
    if (this->Norm() < 1e-12)
        throw GeneralClassicException(
            "Quaternion::normalize()", "length of quaternion less than 1e-12");
#endif

    (*this) /= this->length();

    return (*this);
}

Quaternion Quaternion::inverse() const {
    Quaternion returnValue = conjugate();

    return returnValue.normalize();
}

ippl::Vector<double, 3> Quaternion::rotate(const ippl::Vector<double, 3>& vec) const {
#ifndef NOPAssert
    if (!this->isUnit())
        throw GeneralClassicException(
            "Quaternion::rotate()",
            "quaternion isn't unit quaternion. Norm: " + std::to_string(this->Norm()));
#endif

    Quaternion quat(vec);
    std::cout << "quat (nan?) = " << quat << std::endl;
    std::cout << "quat conj (nan?) = " << (*this).conjugate() << std::endl;

    return ((*this) * (quat * (*this).conjugate())).imag();
}

matrix_t Quaternion::getRotationMatrix() const {
    Quaternion rot(*this);
    matrix_t mat(3, 3);
    double rotNorm = rot.Norm();
    std::cout << "HHeyyyy sooo = " << rotNorm << std::endl;

    // Avoid nan rotation --> nan positions after pushing in ParallelTracker
    if (std::isnan(rotNorm) || rotNorm < 1e-12) {
        std::cout << "yep, rotNorm is nan" << std::endl;
        mat(0, 0) = 1.0;
        mat(0, 1) = 0.0;
        mat(0, 2) = 0.0;
        mat(1, 0) = 0.0;
        mat(1, 1) = 1.0;
        mat(1, 2) = 0.0;
        mat(2, 0) = 0.0;
        mat(2, 1) = 0.0;
        mat(2, 2) = 1.0;

        return mat;
    }

    rot.normalize();

    mat(0, 0) = 1 - 2 * (rot(2) * rot(2) + rot(3) * rot(3));
    mat(0, 1) = 2 * (-rot(0) * rot(3) + rot(1) * rot(2));
    mat(0, 2) = 2 * (rot(0) * rot(2) + rot(1) * rot(3));
    mat(1, 0) = 2 * (rot(0) * rot(3) + rot(1) * rot(2));
    mat(1, 1) = 1 - 2 * (rot(1) * rot(1) + rot(3) * rot(3));
    mat(1, 2) = 2 * (-rot(0) * rot(1) + rot(2) * rot(3));
    mat(2, 0) = 2 * (-rot(0) * rot(2) + rot(1) * rot(3));
    mat(2, 1) = 2 * (rot(0) * rot(1) + rot(2) * rot(3));
    mat(2, 2) = 1 - 2 * (rot(1) * rot(1) + rot(2) * rot(2));

    return mat;
}
