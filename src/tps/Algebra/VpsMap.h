#ifndef TPS_VpsMap_HH
#define TPS_VpsMap_HH
// Copyright: see Copyright.readme

// Class: VpsMap

#include "tps/Algebra/Matrix.h"
#include "tps/Algebra/Tps.h"
#include "tps/Algebra/TpsData.h"
#include "tps/Algebra/TpsMonomial.h"
#include "tps/Algebra/TpsSubstitution.h"
#include "tps/Algebra/Vector.h"
#include "tps/Algebra/Vps.h"
#include "tps/Error.hpp"
#include <iosfwd>

namespace tps {

// Template class VpsMap<T>, declaration.
/** @brief Vector-valued polynomial map with component-wise TPS storage.
 * @details Input and output dimensions may differ. F.substitute(G,p) computes
 * the degree-p projection of F composed with G: G acts first.
 * @see tps_mathematics
 */

template <class T> class VpsMap : public Vps<T> {

  public:
    /// Constructor.
    ///  Construct `nDim` series in `nVar` variables each.
    VpsMap(int nDim, int nVar);

    /// Convert.
    VpsMap(const Vps<T>& rhs);

    /// Convert.
    ///  The constant terms are zero.
    ///  The linear terms are taken from `M`.
    VpsMap(const Matrix<T>& M);

    /// Convert.
    ///  The constant terms are taken from `V`.
    ///  The linear terms are zero.
    VpsMap(const Vector<T>& V);

    VpsMap();
    VpsMap(const VpsMap<T>& rhs);
    ~VpsMap();
    VpsMap<T>& operator=(const VpsMap<T>& y);

    /// Substitute.
    VpsMap<T> substitute(const VpsMap<T>& vv) const;

    /// Substitute.
    VpsMap<T> substitute(const Matrix<T>& M) const;

    /// Substitute and truncate.
    VpsMap<T> substitute(const VpsMap<T>& y, int trunc) const;

    /// Substitute.
    VpsMap<T> substituteInto(const Matrix<T>& x) const;

    /// Evaluate map at point.
    ///  Return the point `y` mapped by the map.
    Vector<T> constantTerm(const Vector<T>& y) const;

    /// Evaluate map at origin.
    ///  Extract constant part.
    Vector<T> constantTerm() const;

    /// Extract linear terms at point.
    ///  Return the linear terms of the map around the point `y`,
    ///  i.e. the linear part of the partial differentials at `y`.
    Matrix<T> linearTerms(const Vector<T>& y) const;

    /// Extract linear terms at origin.
    ///  Return the linear part.
    Matrix<T> linearTerms() const;

    /// Derivative with respect to variable `var`.
    VpsMap<T> derivative(int var) const;

    /// Integral with respect to variable `var`.
    VpsMap<T> integral(int var) const;
};

// Global Operators on VpsMap<T>

/// Multiply.
template <class T> VpsMap<T> operator*(const VpsMap<T>& x, const Tps<T>& y);

/// Multiply.
template <class T> VpsMap<T> operator*(const Tps<T>& x, const VpsMap<T>& y);

/// Multiply.
template <class T> VpsMap<T> operator*(const VpsMap<T>& x, const T& y);

/// Multiply.
template <class T> VpsMap<T> operator*(const T& x, const VpsMap<T>& y);

/// Divide.
///  Throw DivideError if constant part of `y` is zero.
template <class T> VpsMap<T> operator/(const VpsMap<T>& x, const Tps<T>& y);

/// Divide.
///  Throw DivideError if `y` is zero.
template <class T> VpsMap<T> operator/(const VpsMap<T>& x, const T& y);

/// Add.
template <class T> VpsMap<T> operator+(const VpsMap<T>& x, const VpsMap<T>& y);

/// Subtract.
template <class T> VpsMap<T> operator-(const VpsMap<T>& x, const VpsMap<T>& y);

/// Add.
template <class T> VpsMap<T> operator+(const VpsMap<T>& x, const Vector<T>& y);

/// Subtract.
template <class T> VpsMap<T> operator-(const VpsMap<T>& x, const Vector<T>& y);

/// Add.
template <class T> VpsMap<T> operator+(const Vector<T>& x, const VpsMap<T>& y);

/// Subtract.
template <class T> VpsMap<T> operator-(const Vector<T>& x, const VpsMap<T>& y);

/// Extract from stream.
template <class T> std::istream& operator>>(std::istream&, VpsMap<T>& x);

/// Insert to stream.
template <class T> std::ostream& operator<<(std::ostream&, const VpsMap<T>& x);

/// Substitute.
template <class T> inline VpsMap<T> operator*(const Matrix<T>& x, const VpsMap<T>& y) {
    return y.substituteInto(x);
}

// Template class VpsMap<T>, implementation.

template <class T> inline VpsMap<T>::VpsMap() : Vps<T>() {}

template <class T> inline VpsMap<T>::VpsMap(int nDim, int nVar) : Vps<T>(nDim, nVar) {}

template <class T> inline VpsMap<T>::VpsMap(const VpsMap<T>& rhs) : Vps<T>(rhs) {}

template <class T> inline VpsMap<T>::VpsMap(const Vps<T>& rhs) : Vps<T>(rhs) {}

template <class T> inline VpsMap<T>::VpsMap(const Matrix<T>& x) : Vps<T>(x) {}

template <class T> inline VpsMap<T>::VpsMap(const Vector<T>& x) : Vps<T>(x) {}

template <class T> inline VpsMap<T>::~VpsMap() {}

template <class T> inline VpsMap<T>& VpsMap<T>::operator=(const VpsMap<T>& y) {
    Vps<T>::operator=(y);
    return *this;
}

template <class T> inline VpsMap<T> VpsMap<T>::substitute(const VpsMap<T>& y) const {
    return substitute(y, this->getTruncOrder());
}

template <class T> inline VpsMap<T> VpsMap<T>::substitute(const Matrix<T>& y) const {
    return substitute(VpsMap<T>(y));
}

template <class T> Vector<T> VpsMap<T>::constantTerm(const Vector<T>& y) const {
    if (this->getVariables() != y.size()) {
        throw SizeError("VpsMap::constantTerm()", "Inconsistent dimensions.");
    }

    // This statement makes sure that the substitution table exists.
    const Array1D<TpsSubstitution>& table =
        TpsData::getTpsData(this->getTopOrder(), this->getVariables())->getSubTable();

    // Copy the constant term.
    Vector<T> z = constantTerm();
    Array1D<T> product(table.size());
    product[0] = T(1);

    for (int next = 1; next < table.size();) {
        const TpsSubstitution& s = table[next];
        product[s.order] = product[s.order - 1] * y[s.variable];

        for (int v = 0; v < this->getDimension(); v++) {
            if (s.order <= this->data[v].getMaxOrder() && this->data[v][s.index] != T(0)) {
                z[v] += this->data[v][s.index] * product[s.order];
            }
        }

        next = (s.order < this->getTopOrder()) ? next + 1 : s.skip;
    }

    return z;
}

template <class T> Vector<T> VpsMap<T>::constantTerm() const {
    Vector<T> z(this->getDimension());
    for (int i = 0; i < this->getDimension(); i++)
        z[i] = this->data[i][0];
    return z;
}

template <class T> Matrix<T> VpsMap<T>::linearTerms(const Vector<T>& y) const {
    const int variables = this->getVariables();
    if (variables && variables != y.size())
        throw SizeError("VpsMap::linearTerms()", "Point dimension mismatch.");
    // All-scalar maps have no stored input dimension: the point supplies it.
    Matrix<T> z(this->getDimension(), y.size());
    for (int i = 0; i < z.nrows(); i++) {
        for (int j = 0; j < z.ncols(); j++) {
            z[i][j] = (this->data[i].getVariables() == 0 || this->data[i].getMaxOrder() == 0)
                          ? T(0)
                          : this->data[i].derivative(j).evaluate(y);
        }
    }

    return z;
}

template <class T> Matrix<T> VpsMap<T>::linearTerms() const {
    Matrix<T> M(this->getDimension(), this->getVariables());

    for (int i = 0; i < this->getDimension(); i++) {
        for (int j = 0; j < this->getVariables(); j++) {
            M(i, j) = (this->data[i].getVariables() != 0 && this->data[i].getMaxOrder() >= 1)
                          ? this->data[i][j + 1]
                          : T(0);
        }
    }

    return M;
}

template <class T> VpsMap<T> VpsMap<T>::derivative(int var) const {
    VpsMap<T> z(this->getDimension(), this->getVariables());

    for (int i = 0; i < this->getDimension(); i++) {
        z[i] = this->data[i].derivative(var);
    }

    return z;
}

template <class T> VpsMap<T> VpsMap<T>::integral(int var) const {
    VpsMap<T> z(this->getDimension(), this->getVariables());

    for (int i = 0; i < this->getDimension(); i++) {
        z[i] = this->data[i].integral(var);
    }

    return z;
}

template <class T> VpsMap<T> VpsMap<T>::substitute(const VpsMap<T>& y, int trunc) const {
    // Consistency checks.
    this->check();
    y.check();

    if (this->variables != y.getDimension()) {
        throw SizeError("VpsMap::substitute()", "Inconsistent dimensions.");
    }

    // This statement makes sure that the substitution table exists.
    const Array1D<TpsSubstitution>& table =
        TpsData::getTpsData(trunc, this->getVariables())->getSubTable();

    // Highest order present in this Vps.
    int topOrder = this->getTopOrder();
    VpsMap<T> z = constantTerm();

    if (topOrder) {
        Array1D<Tps<T>> product(topOrder + 1);
        product[0] = Tps<T>(T(1));

        for (int next = 1; next < table.size();) {
            const TpsSubstitution& s = table[next];
            product[s.order] = product[s.order - 1].multiply(y[s.variable], trunc);

            for (int v = 0; v < this->getDimension(); v++) {
                int maxOrder = this->data[v].getMaxOrder();
                int cutOrder = std::min(this->data[v].getTruncOrder(), trunc);

                if (s.order <= maxOrder && this->data[v][s.index] != T(0)) {
                    z.data[v] += this->data[v][s.index] * product[s.order].truncate(cutOrder);
                }
            }

            next = (s.order < topOrder) ? next + 1 : s.skip;
        }
    }

    return z;
}

template <class T> VpsMap<T> VpsMap<T>::substituteInto(const Matrix<T>& x) const {
    if (x.ncols() != this->getDimension()) {
        throw SizeError("VpsMap::substituteInto()", "Inconsistent dimensions.");
    }

    VpsMap<T> z(x.nrows(), this->getVariables());

    for (int i = 0; i < x.nrows(); i++) {
        for (int j = 0; j < this->data.size(); j++) {
            z[i] += this->data[j] * x(i, j);
        }
    }

    return z;
}

// Implementation of global functions for class Vps<T>.

template <class T> inline VpsMap<T> operator*(const VpsMap<T>& x, const Tps<T>& y) {
    VpsMap<T> z(x);
    return z *= y;
}

template <class T> inline VpsMap<T> operator*(const Tps<T>& x, const VpsMap<T>& y) {
    VpsMap<T> z(x);
    return z *= y;
}

template <class T> inline VpsMap<T> operator*(const VpsMap<T>& x, const T& y) {
    VpsMap<T> z(x);
    return z *= y;
}

template <class T> inline VpsMap<T> operator*(const T& x, const VpsMap<T>& y) {
    VpsMap<T> z(x);
    return z *= y;
}

template <class T> inline VpsMap<T> operator/(const VpsMap<T>& x, const Tps<T>& y) {
    VpsMap<T> z(x);
    return z /= y;
}

template <class T> inline VpsMap<T> operator/(const VpsMap<T>& x, const T& y) {
    VpsMap<T> z(x);
    return z /= y;
}

template <class T> inline VpsMap<T> operator+(const VpsMap<T>& x, const VpsMap<T>& y) {
    VpsMap<T> z(x);
    return z += y;
}

template <class T> inline VpsMap<T> operator-(const VpsMap<T>& x, const VpsMap<T>& y) {
    VpsMap<T> z(x);
    return z -= y;
}

template <class T> inline VpsMap<T> operator+(const VpsMap<T>& x, const Vector<T>& y) {
    VpsMap<T> z(x);
    return z += y;
}

template <class T> inline VpsMap<T> operator-(const VpsMap<T>& x, const Vector<T>& y) {
    VpsMap<T> z(x);
    return z -= y;
}

template <class T> inline VpsMap<T> operator+(const Vector<T>& x, const VpsMap<T>& y) {
    VpsMap<T> z(y);
    return z += x;
}

template <class T> inline VpsMap<T> operator-(const Vector<T>& x, const VpsMap<T>& y) {
    VpsMap<T> z(-y);
    return z += x;
}

template <class T> inline std::istream& operator>>(std::istream& is, VpsMap<T>& x) {
    return x.get(is);
}

template <class T> inline std::ostream& operator<<(std::ostream& os, const VpsMap<T>& x) {
    return x.put(os);
}

} // namespace tps

#endif // TPS_VpsMap_HH
