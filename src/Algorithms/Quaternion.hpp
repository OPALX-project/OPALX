#ifndef OPAL_QUATERNION_H
#define OPAL_QUATERNION_H

#include "Algorithms/Matrix.h"
#include "Ippl.h"

/**
 * @class Quaternion
 * @brief Quaternion storage and rotation algebra used by OPALX geometry code.
 *
 * The quaternion is stored as
 * \f[
 * q = (w, x, y, z) = \left(w, \mathbf{v}\right)
 * \f]
 * with scalar part \f$w\f$ and vector part \f$\mathbf{v} = (x, y, z)\f$.
 *
 * For a unit quaternion, the active rotation of a spatial vector
 * \f$\mathbf{r}\f$ is computed as
 * \f[
 * \mathbf{r}' = \left(q \, (0, \mathbf{r}) \, q^\ast\right)_{\mathrm{imag}},
 * \f]
 * where \f$q^\ast\f$ is the quaternion conjugate and the imaginary part is
 * mapped back to \f$\mathbb{R}^3\f$.
 *
 * OPALX uses quaternions as the rotational kernel of
 * CoordinateSystemTrafo. The conversion to a rotation matrix is defined for
 * unit quaternions; non-unit inputs are normalized before the matrix is built.
 */
class Quaternion : public ippl::Vector<double, 4> {
public:
    /* ======== Constructors ======== */
    /// Construct the identity rotation quaternion \f$(1, 0, 0, 0)\f$.
    Quaternion();
    Quaternion(const Quaternion&);
    /// Construct a quaternion from its four components \f$(w, x, y, z)\f$.
    Quaternion(const double&, const double&, const double&, const double&);
    /// Construct a pure quaternion \f$(0, \mathbf{v})\f$ from a 3-vector.
    Quaternion(const ippl::Vector<double, 3>&);
    /// Construct a quaternion from scalar part \f$w\f$ and vector part \f$\mathbf{v}\f$.
    Quaternion(const double&, const ippl::Vector<double, 3>&);
    /// Construct a quaternion from a 3x3 rotation matrix.
    Quaternion(const matrix3x3_t&);

    /**
     * @brief Return the multiplicative inverse.
     *
     * The inverse is
     * \f[
     * q^{-1} = \frac{q^\ast}{\| q \|^2}.
     * \f]
     */
    Quaternion inverse() const;

    /// Return the quaternion conjugate \f$q^\ast = (w, -\mathbf{v})\f$.
    Quaternion conjugate() const;

    /* ======== Operators ======== */
    Quaternion operator*(const double&) const;

    /**
     * @brief Hamilton product of two quaternions.
     *
     * For \f$q_1 = (w_1, \mathbf{v}_1)\f$ and
     * \f$q_2 = (w_2, \mathbf{v}_2)\f$,
     * \f[
     * q_1 q_2 =
     * \left(w_1 w_2 - \mathbf{v}_1 \cdot \mathbf{v}_2,
     * w_1 \mathbf{v}_2 + w_2 \mathbf{v}_1 + \mathbf{v}_1 \times \mathbf{v}_2\right).
     * \f]
     */
    Quaternion operator*(const Quaternion&) const;
    Quaternion& operator=(const Quaternion&) = default;
    Quaternion& operator*=(const Quaternion&);
    Quaternion operator/(const double&) const;

    /// Return \f$\| q \|^2 = w^2 + x^2 + y^2 + z^2\f$.
    double Norm() const;

    /// Return the Euclidean norm \f$\| q \|\f$.
    double length() const;

    /**
     * @brief Normalize the quaternion in place.
     *
     * After normalization,
     * \f[
     * \| q \| = 1.
     * \f]
     */
    Quaternion& normalize();

    /// Return true if \f$\| q \|^2 \approx 1\f$ within a fixed tolerance.
    bool isUnit() const;

    /// Return true if the scalar part is approximately zero.
    bool isPure() const;

    /// Return true if the quaternion is both pure and unit length.
    bool isPureUnit() const;

    /// Return the scalar part \f$w\f$.
    double real() const;

    /// Return the vector part \f$\mathbf{v} = (x, y, z)\f$.
    ippl::Vector<double, 3> imag() const;

    /**
     * @brief Rotate a spatial vector by a unit quaternion.
     *
     * This method requires a unit quaternion and evaluates
     * \f[
     * \mathbf{r}' = \left(q \, (0, \mathbf{r}) \, q^\ast\right)_{\mathrm{imag}}.
     * \f]
     */
    ippl::Vector<double, 3> rotate(const ippl::Vector<double, 3>&) const;

    /**
     * @brief Convert the quaternion to a 3x3 rotation matrix.
     *
     * The quaternion is normalized first, so the returned matrix represents the
     * same rotation for any non-zero scalar multiple of \f$q\f$.
     */
    matrix3x3_t getRotationMatrix() const;
};

typedef Quaternion Quaternion_t;

/**
 * @brief Return a unit quaternion that rotates @p vec onto @p reference.
 *
 * Both input vectors are normalized internally. For the generic case, the
 * quaternion is built from the rotation axis
 * \f$\mathbf{a} \propto \mathbf{vec} \times \mathbf{reference}\f$ and the
 * half-angle formula. For antiparallel vectors the current implementation
 * chooses an arbitrary orthogonal axis, so the 180-degree rotation is valid but
 * not deterministic.
 */
Quaternion getQuaternion(ippl::Vector<double, 3> vec, ippl::Vector<double, 3> reference);

inline Quaternion::Quaternion() : Vector<double, 4>(1.0, 0.0, 0.0, 0.0) { srand(time(0)); }

inline Quaternion::Quaternion(const Quaternion& quat) : Vector<double, 4>(quat) { srand(time(0)); }

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

inline double Quaternion::Norm() const { return dot(*this); }

inline double Quaternion::length() const { return std::sqrt(this->Norm()); }

inline bool Quaternion::isUnit() const { return (std::abs(this->Norm() - 1.0) < 1e-12); }

inline bool Quaternion::isPure() const { return (std::abs((*this)(0)) < 1e-12); }

inline bool Quaternion::isPureUnit() const { return (this->isPure() && this->isUnit()); }

inline Quaternion Quaternion::conjugate() const {
    Quaternion quat(this->real(), -this->imag());

    return quat;
}

inline double Quaternion::real() const { return (*this)(0); }

inline ippl::Vector<double, 3> Quaternion::imag() const {
    ippl::Vector<double, 3> vec((*this)(1), (*this)(2), (*this)(3));

    return vec;
}

#endif
