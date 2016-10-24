#ifndef OPAL_VEKTOR_HH
#define OPAL_VEKTOR_HH

#include <limits>

#include "AppTypes/Vektor.h"

template <class, unsigned>
class Tenzor;

typedef Vektor<double, 3> Vector_t;

void normalize(Vector_t & vec);

class Quaternion: public Vektor<double, 4> {
public:
    Quaternion();
    Quaternion(const Quaternion &);
    Quaternion(const double &, const double &, const double &, const double &);
    Quaternion(const Vector_t &);
    Quaternion(const double &, const Vector_t &);

    const Quaternion operator*(const double &) const;
    const Quaternion operator*(const Quaternion &) const;
    Quaternion& operator*=(const Quaternion &);
    const Quaternion operator/(const double &) const;

    double Norm() const;
    double length() const;
    Quaternion & normalize();

    bool isUnit() const;
    bool isPure() const;
    bool isPureUnit() const;

    const Quaternion inverse() const;
    const Quaternion conjugate() const;

    double real() const;
    const Vector_t imag() const;

    const Vector_t rotate(const Vector_t &) const;

    const Tenzor<double, 3> getRotationMatrix() const;
};

typedef Quaternion Quaternion_t;
//typedef Vektor<double, 4> Quaternion_t;  // used in ParallelCyclotron.{cpp,h}

Quaternion getQuaternion(Vector_t u, Vector_t v);


inline
Quaternion::Quaternion():
    Vektor<double, 4>(0.0)
{}

inline
Quaternion::Quaternion(const Quaternion & quat):
    Vektor<double, 4>(quat)
{}

inline
Quaternion::Quaternion(const double & x0, const double & x1, const double & x2, const double & x3):
    Vektor<double, 4>(x0, x1, x2, x3)
{}

inline
Quaternion::Quaternion(const Vector_t & vec):
    Quaternion(0.0, vec(0), vec(1), vec(2))
{}

inline
Quaternion::Quaternion(const double & realPart, const Vector_t & vec):
    Quaternion(realPart, vec(0), vec(1), vec(2))
{}

inline
double Quaternion::Norm() const
{
    return dot(*this, *this);
}

inline
double Quaternion::length() const
{
    return sqrt(this->Norm());
}

inline
bool Quaternion::isUnit() const
{
    return (std::abs(this->Norm() - 1.0) < 1e-12);
}

inline
bool Quaternion::isPure() const
{
    return (std::abs((*this)(0)) < 1e-12);
}

inline
bool Quaternion::isPureUnit() const
{
    return (this->isPure() && this->isUnit());
}

inline
const Quaternion Quaternion::conjugate() const
{
    Quaternion quat(this->real(), -this->imag());

    return quat;
}

inline
double Quaternion::real() const
{
    return (*this)(0);
}

inline
const Vector_t Quaternion::imag() const
{
    Vector_t vec((*this)(1), (*this)(2), (*this)(3));

    return vec;
}


#endif
