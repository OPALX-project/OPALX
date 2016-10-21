#include "Algorithms/Vektor.h"
#include "AppTypes/Tenzor.h"

void normalize(Vector_t & vec)
{
    double length = sqrt(dot(vec, vec));

    PAssert(length > 1e-12);

    vec /= length;
}

Quaternion getQuaternion(Vector_t u, Vector_t v)
{
    const double tol1 = 1e-5;
    const double tol2 = 1e-8;
    Vector_t imaginaryPart(0.0);
    normalize(u);
    normalize(v);

    double k = 1.0;
    double cosTheta = dot(u, v);

    if (std::abs(cosTheta + 1.0) < tol1) {
        Vector_t xAxis(1.0, 0, 0);
        Vector_t zAxis(0, 0, 1.0);

        if (std::abs(1.0 - dot(u, xAxis)) > tol2)
            imaginaryPart = cross(u, xAxis);
        else
            imaginaryPart = cross(u, zAxis);

        k = 0.0;
        cosTheta = 0.0;
    } else {
        imaginaryPart = cross(u, v);
    }

    Quaternion result(cosTheta + k, imaginaryPart);
    result.normalize();

    return result;
}

const Quaternion Quaternion::operator*(const double & d) const
{
    Quaternion result(d * this->real(), d * this->imag());

    return result;
}

const Quaternion Quaternion::operator*(const Quaternion & other) const
{
    Quaternion result(*this);
    return result *= other;
}

Quaternion & Quaternion::operator*=(const Quaternion & other)
{
    Vector_t imagThis = this->imag();
    Vector_t imagOther = other.imag();

    *this = Quaternion((*this)(0) * other(0) - dot(imagThis, imagOther),
                       (*this)(0) * imagOther + other(0) * imagThis + cross(imagThis, imagOther));

    return *this;
}

const Quaternion Quaternion::operator/(const double & d) const
{
    Quaternion result(this->real() / d, this->imag() / d);

    return result;
}

Quaternion & Quaternion::normalize()
{
    PAssert(this->Norm() > 1e-12);

    (*this) /= this->length();

    return (*this);
}

const Quaternion Quaternion::inverse() const
{
    PAssert(this->Norm() > 1e-12);

    Quaternion returnValue = conjugate();

    return returnValue.normalize();
}

const Vector_t Quaternion::rotate(const Vector_t & vec) const
{
    PAssert(this->isUnit());

    Quaternion quat(vec);

    return ((*this) * (quat * (*this).conjugate())).imag();
}

const Tenzor<double, 3> Quaternion::getRotationMatrix() const
{
    Tenzor<double, 3> mat(1 - 2 * ((*this)(2)*(*this)(2) + (*this)(3)*(*this)(3)),
                          2 * ((*this)(0) * (*this)(3) + (*this)(1) * (*this)(2)),
                          2 * (-(*this)(0) * (*this)(2) + (*this)(1) * (*this)(3)),
                          2 * (-(*this)(0) * (*this)(3) + (*this)(1) * (*this)(2)),
                          1 - 2 * ((*this)(1)*(*this)(1) + (*this)(3)*(*this)(3)),
                          2 * ((*this)(0) * (*this)(1) + (*this)(2) * (*this)(3)),
                          2 * ((*this)(0) * (*this)(2) + (*this)(1) * (*this)(3)),
                          2 * (-(*this)(0) * (*this)(1) + (*this)(2) * (*this)(3)),
                          1 - 2 * ((*this)(1)*(*this)(1) + (*this)(2)*(*this)(2)));

    return mat;
}
