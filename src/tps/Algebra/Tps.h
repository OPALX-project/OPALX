#ifndef TPS_Tps_HH
#define TPS_Tps_HH
#include "Utilities/CLRangeError.h"
#include "Utilities/DivideError.h"
#include "Utilities/LogicalError.h"
// Copyright: see NOTICE and LICENSE

// Template Class: Tps<class T>
//   Truncated power series in n variables of type T.

#include <climits>
#include <iosfwd>
#include <limits>
#include <span>

namespace tps {

    template <class T>
    class Tps;
    template <class T>
    class TpsRep;
    template <class T>
    class Matrix;
    template <class T>
    class Vector;
    template <class T>
    class VpsMap;

    class TpsMonomial;

    // Class Tps
    /** @brief Dense multivariate polynomial or finite Taylor jet.
     * @tparam T Coefficient type supporting operations used by the instantiated algorithm.
     * @details Stores ordinary coefficients \f$c_\alpha\f$ of \f$x^\alpha\f$;
     * derivatives are divided by \f$\alpha!\f$. Storage degree and truncation
     * precision differ. Copy-on-write coefficients are owned and correctly aligned.
     * @see tps_mathematics tps_numerics
     * @warning Global limits and the shared index cache require external synchronization.
     */

    template <class T>
    class Tps {
    public:
        /// Constructor.
        ///  Define maximum order and number of variables, value = 0.
        Tps(int maxOrder, int nVar);

        /// Conversion.
        ///  From constant `y`.
        Tps(const T& y);

        /// Conversion.
        ///  From int `y`.
        Tps(int y);

        Tps();
        Tps(const Tps<T>& x);
        ~Tps();
        Tps<T>& operator=(const Tps<T>& y);

        /// Convert and assign.
        Tps<T>& operator=(const T& y);

        /// Extract orders.
        ///  Return the orders `lowOrder` through `highOrder`,
        ///  both included.
        Tps<T> filter(int lowOrder, int highOrder) const;

        /// Truncate.
        ///  Change the maximum order to `trunc`; may also reserve more space.
        Tps<T> truncate(int trunc);

        /// Get coefficient.
        ///  Return the coefficient denoted by the Giorgilli index.
        ///  Throw ::CLRangeError, if index is out of range.
        const T getCoefficient(int index) const;

        /// Set coefficient.
        ///  Assign the coefficient denoted by the Giorgilli index.
        ///  Expand size if index is out of range.
        void setCoefficient(int index, const T& value);

        /// Get coefficient.
        ///  Return the coefficient denoted by monomial exponents.
        ///  Throw ::CLRangeError, if index is out of range.
        const T getCoefficient(const TpsMonomial& monomial) const;

        /// Set coefficient.
        ///  Assign the coefficient denoted by monomial exponents.
        ///  Throw ::CLRangeError, if index is out of range.
        void setCoefficient(const TpsMonomial& monomial, const T& value);

        /// Get coefficient.
        ///  Return the coefficient denoted by the Giorgilli index.
        ///  Throws ::CLRangeError outside allocated coefficients; does not expand storage.
        const T operator[](int index) const;

        /// Get coefficient.
        ///  Return a reference to the coefficient denoted by the Giorgilli index.
        ///  Throws ::CLRangeError outside allocated coefficients; does not expand storage.
        T& operator[](int index);

        /// Get coefficient.
        ///  Return the coefficient denoted by monomial exponents.
        ///  Throws ::CLRangeError outside allocated coefficients; does not expand storage.
        ///  Throws ::LogicalError for mismatched variable counts.
        const T operator[](const TpsMonomial& monomial) const;

        /// Get coefficient.
        ///  Return reference to the coefficient denoted by monomial exponents.
        ///  Throws ::CLRangeError outside allocated coefficients; does not expand storage.
        ///  Throws ::LogicalError for mismatched variable counts.
        T& operator[](const TpsMonomial& monomial);

        /// Make variable.
        ///  Construct the variable identified by the index `var`,
        ///  with total of `nVar` variables.
        static Tps<T> makeVariable(int nVar, int var);

        /// Make power.
        ///  Construct power of variable `var`, with total of `nVar`
        ///  variables.
        static Tps<T> makeVarPower(int nVar, int var, int power);

        /// Make monomial.
        ///  Construct the monomial with the exponents in `m` and coefficient
        ///  `t`.
        static Tps<T> makeMonomial(const TpsMonomial& m, const T& t);

        /// Unary plus.
        Tps<T> operator+() const;

        /// Unary minus.
        Tps<T> operator-() const;

        /// Add and assign.
        Tps<T>& operator+=(const Tps<T>& y);

        /// Subtract and assign.
        Tps<T>& operator-=(const Tps<T>& y);

        /// Multiply and assign.
        Tps<T>& operator*=(const Tps<T>& y);

        /// Divide and assign.
        Tps<T>& operator/=(const Tps<T>& y);

        /// Add constant and assign.
        Tps<T>& operator+=(const T& y);

        /// Subtract constant and assign.
        Tps<T>& operator-=(const T& y);

        /// Multiply by constant and assign.
        Tps<T>& operator*=(const T& y);

        /// Divide by constant and assign.
        Tps<T>& operator/=(const T& y);

        /// Equality operator.
        bool operator==(const Tps<T>& y) const;

        /// Equlatity with constant.
        bool operator==(const T& y) const;

        /// Inequality operator.
        bool operator!=(const Tps<T>& y) const;

        /// Inequality with constant.
        bool operator!=(const T& y) const;

        /// Substitute.
        ///  Substitute matrix `M` in Tps, giving new Tps.
        Tps<T> substitute(const Matrix<T>& M) const;

        /// Substitute.
        ///  Substitute map `m` in Tps, giving new Tps.
        Tps<T> substitute(const VpsMap<T>& m) const;

        /// Substitute.
        ///  Substitute vector `v` in Tps, giving a constant.
        T evaluate(const Vector<T>& v) const;

        /// Set to zero.
        void clear();

        /** @brief Finite Taylor jet of the reciprocal.
         * @param order Requested order, constrained by representation/global precision.
         * @return Reciprocal jet; EXACT resolves to a finite stored order.
         * @throws ::DivideError If the constant coefficient is zero.
         * @throws ::CLRangeError If order is negative.
         * @see tps_mathematics
         */
        Tps<T> inverse(int order = truncOrder) const;

        /// Truncated multiplication.
        Tps<T> multiply(const Tps<T>& y, int trunc) const;

        /// Partial derivative.
        ///  Return partial derivative with respect to variable `var`.
        ///  Return zero for a constant.
        Tps<T> derivative(int var) const;

        /// Partial integral.
        ///  Return partial integral with respect to variable `var`.
        ///  Throw ::LogicalError for a constant.
        Tps<T> integral(int var) const;

        /// Multiply by variable `var`.
        ///  Throw ::LogicalError for a constant.
        Tps<T> multiplyVariable(int var) const;

        /// Multiply monomial-wise.
        ///  Construct new Tps<T> by multipying corresponding monomials pair-wise.
        Tps<T> scaleMonomials(const Tps<T>& y) const;

        /** @brief Evaluate a Taylor polynomial in the zero-constant part of this series.
         * @param coefficients Coefficients of successive powers of
         * \f$u=f-f(0)\f$, including the constant coefficient.
         * @return The Horner-evaluated polynomial, subject to truncation limits.
         * @throws ::CLRangeError For an empty or unrepresentably large coefficient span.
         * @note Coefficients must already include factorial normalization.
         */
        [[nodiscard]] Tps<T> taylor(std::span<const T> coefficients) const;

        /// Get maximal order.
        [[nodiscard]] int getMaxOrder() const;

        /// Get truncation order.
        [[nodiscard]] int getTruncOrder() const;

        /// Get number of variables.
        [[nodiscard]] int getVariables() const;

        /// Get number of coefficients.
        int getSize() const;

        /// Get Tps from the stream is.
        std::istream& get(std::istream& is);

        /// Put Tps to the stream is.
        std::ostream& put(std::ostream& os) const;

        /// Get exponents.
        ///  Return the exponent array for the given Giorgilli index.
        const TpsMonomial& getExponents(int index) const;

        /// Get order.
        ///  Return the order of the monomial with the given Giorgilli index.
        int getOrder(int index) const;

        /// Size.
        ///  Return the number of monomials required for a Tps of order `order`.
        int getSize(int order) const;

        /// Test for constant.
        bool isConstant() const;

        /// Get global truncation order.
        static int getGlobalTruncOrder();

        /// Finite order for Taylor functions: explicit limits, or stored order for EXACT.
        /// Constants need only the scalar (order-zero) term. Negative orders throw.
        [[nodiscard]] int getTaylorOrder(int order = EXACT) const;

        /// Set global truncation order (nonnegative, or EXACT).
        static void setGlobalTruncOrder(int order);

        /// Representation of infinite precision.
        static constexpr int EXACT = std::numeric_limits<int>::max();

    private:
        // Validate dimensions and stored degree before monomial cache lookup.
        int checkedCoefficientIndex(const TpsMonomial& monomial) const;

        // Internal constructor.
        Tps(TpsRep<T>*);

        // Set maximum order.
        void setMaxOrder(int);

        // Make representation unique.
        void unique();

        // Pointer to representation.
        TpsRep<T>* rep;

        // Current maximum truncation order.
        inline static int truncOrder = EXACT;
    };

    // Global functions acting on Tps<T> objects.

    /// Add.
    template <class T>
    Tps<T> operator+(const Tps<T>& x, const Tps<T>& y);

    /// Subtract.
    template <class T>
    Tps<T> operator-(const Tps<T>& x, const Tps<T>& y);

    /// Add.
    template <class T>
    Tps<T> operator+(const Tps<T>& x, const T& y);

    /// Subtract.
    template <class T>
    Tps<T> operator-(const Tps<T>& x, const T& y);

    /// Add.
    template <class T>
    Tps<T> operator+(const T& x, const Tps<T>& y);

    /// Subtract.
    template <class T>
    Tps<T> operator-(const T& x, const Tps<T>& y);

    /// Multiply.
    template <class T>
    Tps<T> operator*(const Tps<T>& x, const Tps<T>& y);

    /// Divide.
    template <class T>
    Tps<T> operator/(const Tps<T>& x, const Tps<T>& y);

    /// Multiply.
    template <class T>
    Tps<T> operator*(const Tps<T>& x, const T& y);

    /// Divide.
    template <class T>
    Tps<T> operator/(const Tps<T>& x, const T& y);

    /// Multiply.
    template <class T>
    Tps<T> operator*(const T& x, const Tps<T>& y);

    /// Divide.
    template <class T>
    Tps<T> operator/(const T& x, const Tps<T>& y);

    /// Test for equality.
    template <class T>
    bool operator==(const T& x, const Tps<T>& y);

    /// Extract from stream.
    template <class T>
    std::istream& operator>>(std::istream&, Tps<T>& x);

    /// Insert to stream.
    template <class T>
    std::ostream& operator<<(std::ostream&, const Tps<T>& x);

    // Inline global functions.

    template <class T>
    Tps<T> operator+(const Tps<T>& x, const Tps<T>& y) {
        Tps<T> z(x);
        return z += y;
    }

    template <class T>
    Tps<T> operator-(const Tps<T>& x, const Tps<T>& y) {
        Tps<T> z(x);
        return z -= y;
    }

    template <class T>
    inline Tps<T> operator+(const Tps<T>& x, const T& y) {
        Tps<T> z(x);
        return z += y;
    }

    template <class T>
    inline Tps<T> operator-(const Tps<T>& x, const T& y) {
        Tps<T> z(x);
        return z -= y;
    }

    template <class T>
    inline Tps<T> operator+(const T& x, const Tps<T>& y) {
        Tps<T> z(y);
        return z += x;
    }

    template <class T>
    Tps<T> operator-(const T& x, const Tps<T>& y) {
        Tps<T> z(-y);
        return z += x;
    }

    template <class T>
    Tps<T> operator*(const Tps<T>& x, const Tps<T>& y) {
        Tps<T> z(x);
        return z *= y;
    }

    template <class T>
    Tps<T> operator/(const Tps<T>& x, const Tps<T>& y) {
        Tps<T> z(x);
        return z /= y;
    }

    template <class T>
    Tps<T> operator*(const Tps<T>& x, const T& y) {
        Tps<T> z(x);
        return z *= y;
    }

    template <class T>
    Tps<T> operator/(const Tps<T>& x, const T& y) {
        Tps<T> z(x);
        return z /= y;
    }

    template <class T>
    Tps<T> operator*(const T& x, const Tps<T>& y) {
        Tps<T> z(y);
        return z *= x;
    }

    template <class T>
    Tps<T> operator/(const T& x, const Tps<T>& y) {
        Tps<T> z(y.inverse(Tps<T>::getGlobalTruncOrder()));
        return z *= x;
    }

    template <class T>
    bool operator==(const T& x, const Tps<T>& y) {
        return y == x;
    }

    template <class T>
    std::istream& operator>>(std::istream& is, Tps<T>& x) {
        return x.get(is);
    }

    template <class T>
    std::ostream& operator<<(std::ostream& os, const Tps<T>& x) {
        return x.put(os);
    }

}  // namespace tps

#endif  // TPS_Tps_HH
