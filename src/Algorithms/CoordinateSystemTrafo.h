#ifndef COORDINATESYSTEMTRAFO
#define COORDINATESYSTEMTRAFO

#include <Kokkos_Core.hpp>
#include "Algorithms/Matrix.h"
#include "Algorithms/Quaternion.hpp"

/**
 * @class CoordinateSystemTrafo
 * @brief Rigid spatial transform between a parent frame and a local frame.
 *
 * OPALX uses this class as a rigid transform with translation and rotation.
 * The stored origin \f$\mathbf{o}\f$ is the location of the local origin
 * expressed in the parent frame. The stored rotation \f$R\f$ maps parent-frame
 * vectors into local-frame vectors.
 *
 * For a point \f$\mathbf{r}\f$ in the parent frame, the local coordinates are
 * computed as
 * \f[
 * \mathbf{r}_{\mathrm{local}} = R \left(\mathbf{r}_{\mathrm{parent}} - \mathbf{o}\right).
 * \f]
 *
 * The inverse mapping is
 * \f[
 * \mathbf{r}_{\mathrm{parent}} = R^T \mathbf{r}_{\mathrm{local}} + \mathbf{o}.
 * \f]
 *
 * The same convention is used for pure vector rotations without translation:
 * \f[
 * \mathbf{v}_{\mathrm{local}} = R \mathbf{v}_{\mathrm{parent}}, \qquad
 * \mathbf{v}_{\mathrm{parent}} = R^T \mathbf{v}_{\mathrm{local}}.
 * \f]
 *
 * Composition follows function composition of transformTo():
 * \f[
 * (A * B)(\mathbf{r}) = A(B(\mathbf{r})).
 * \f]
 * In other words, `A * B` applies `B` first and then `A`.
 */
class CoordinateSystemTrafo {
public:
    /* ============================== Constructors ============================== */
    /// Construct the identity transform.
    CoordinateSystemTrafo();

    CoordinateSystemTrafo(const CoordinateSystemTrafo& right);

    /**
     * @brief Construct a parent-to-local transform.
     *
     * @param origin Position of the local origin in parent coordinates.
     * @param orientation Unit quaternion whose rotation matrix maps
     *                    parent-frame vectors into local-frame vectors.
     */
    CoordinateSystemTrafo(const ippl::Vector<double, 3>& origin, const Quaternion& orientation);

    /// Invert the transform in place.
    void invert();

    /// Return the inverse transform.
    CoordinateSystemTrafo inverted() const;
    /* ========================================================================== */
    /* ======================== Transformation Functions ======================== */
    /// Map a point from the parent frame to the local frame.
    ippl::Vector<double, 3> transformTo(const ippl::Vector<double, 3>& r) const;

    /// Map a point from the local frame back to the parent frame.
    ippl::Vector<double, 3> transformFrom(const ippl::Vector<double, 3>& r) const;

    /// Rotate a vector from the parent frame into the local frame.
    ippl::Vector<double, 3> rotateTo(const ippl::Vector<double, 3>& r) const;

    /// Rotate a vector from the local frame back into the parent frame.
    ippl::Vector<double, 3> rotateFrom(const ippl::Vector<double, 3>& r) const;

    /**
     * @brief Apply transformTo() to a Kokkos view of particle positions.
     *
     * The operation is pointwise and uses the same rigid-transform convention as
     * transformTo().
     */
    template <typename ViewType>
    void transformBunchTo(ViewType Rview, size_t nLocal) const;

    /// Apply transformFrom() to a Kokkos view of particle positions.
    template <typename ViewType>
    void transformBunchFrom(ViewType Rview, size_t nLocal) const;

    /// Apply rotateTo() to a Kokkos view of vectors such as momenta.
    template <typename ViewType>
    void rotateBunchTo(ViewType Pview, size_t nLocal) const;

    /// Apply rotateFrom() to a Kokkos view of vectors such as momenta.
    template <typename ViewType>
    void rotateBunchFrom(ViewType Pview, size_t nLocal) const;
    /* ========================================================================== */
    /* =============================== Operators ================================ */
    CoordinateSystemTrafo& operator=(const CoordinateSystemTrafo& right) = default;

    /**
     * @brief Compose two transforms.
     *
     * For a point \f$\mathbf{r}\f$ in the parent frame,
     * \f[
     * (A * B).transformTo(\mathbf{r}) = A.transformTo(B.transformTo(\mathbf{r})).
     * \f]
     */
    CoordinateSystemTrafo operator*(const CoordinateSystemTrafo& right) const;

    void operator*=(const CoordinateSystemTrafo& right);
    /* ========================================================================== */
    /* =============================== Getters ================================== */
    ippl::Vector<double, 3> getOrigin() const;
    Quaternion getRotation() const;
    matrix3x3_t getRotationMatrix() const;
    /* ========================================================================== */
    /* =============================== Print ==================================== */
    void print(std::ostream&) const;

private:
    /// Local origin expressed in the parent frame.
    ippl::Vector<double, 3> origin_m;

    /// Quaternion whose rotation maps parent-frame vectors into local-frame vectors.
    Quaternion orientation_m;

    /// Cached rotation matrix for the same parent-to-local mapping.
    matrix3x3_t rotationMatrix_m;
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
    return prod_vector(rotationMatrix_m, delta);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::transformFrom(
        const ippl::Vector<double, 3>& r) const {
    return rotateFrom(r) + origin_m;
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateTo(
        const ippl::Vector<double, 3>& r) const {
    return prod_vector(rotationMatrix_m, r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::rotateFrom(
        const ippl::Vector<double, 3>& r) const {
    return prod_vector_transpose(rotationMatrix_m, r);
}

inline ippl::Vector<double, 3> CoordinateSystemTrafo::getOrigin() const { return origin_m; }

inline Quaternion CoordinateSystemTrafo::getRotation() const { return orientation_m; }

inline matrix3x3_t CoordinateSystemTrafo::getRotationMatrix() const { return rotationMatrix_m; }

inline CoordinateSystemTrafo CoordinateSystemTrafo::inverted() const {
    CoordinateSystemTrafo result(*this);
    result.invert();

    return result;
}

template <typename ViewType>
inline void CoordinateSystemTrafo::transformBunchTo(ViewType Rview, size_t nLocal) const {
    matrix3x3_t rot             = rotationMatrix_m;
    ippl::Vector<double, 3> ori = origin_m;
    Kokkos::parallel_for(
            "transformBunchTo", nLocal, KOKKOS_LAMBDA(const size_t i) {
                ippl::Vector<double, 3> delta = Rview(i) - ori;
                Rview(i)                      = prod_vector(rot, delta);
            });
}

template <typename ViewType>
inline void CoordinateSystemTrafo::transformBunchFrom(ViewType Rview, size_t nLocal) const {
    matrix3x3_t rot             = rotationMatrix_m;
    ippl::Vector<double, 3> ori = origin_m;
    Kokkos::parallel_for(
            "transformBunchFrom", nLocal, KOKKOS_LAMBDA(const size_t i) {
                Rview(i) = prod_vector_transpose(rot, Rview(i)) + ori;
            });
}

template <typename ViewType>
inline void CoordinateSystemTrafo::rotateBunchTo(ViewType Pview, size_t nLocal) const {
    matrix3x3_t rot = rotationMatrix_m;
    Kokkos::parallel_for(
            "rotateBunchTo", nLocal,
            KOKKOS_LAMBDA(const size_t i) { Pview(i) = prod_vector(rot, Pview(i)); });
}

template <typename ViewType>
inline void CoordinateSystemTrafo::rotateBunchFrom(ViewType Pview, size_t nLocal) const {
    matrix3x3_t rot = rotationMatrix_m;
    Kokkos::parallel_for(
            "rotateBunchFrom", nLocal,
            KOKKOS_LAMBDA(const size_t i) { Pview(i) = prod_vector_transpose(rot, Pview(i)); });
}

#endif
