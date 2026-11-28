#ifndef TPS_TpsMonomial_HH
#define TPS_TpsMonomial_HH
// Copyright: see LICENSE

// Class: TpsMonomial

#include "tps/Algebra/Array1D.h"

namespace tps {

    // Class TpsMonomial
    /// Exponent array for Tps<T>.
    ///  A representation of the exponents for a Tps monomial.

    class TpsMonomial {
    public:
        /// Constructor.
        ///  Array of zeros with length `nVar`, defines `nVar` variables
        ///  with zero exponent each.
        TpsMonomial(int nVar);

        /// Constructor.
        ///  Array with length `nVar`, define variable with index `var`.
        TpsMonomial(int nVar, int var);

        TpsMonomial()                              = default;
        TpsMonomial(const TpsMonomial&)            = default;
        ~TpsMonomial()                             = default;
        TpsMonomial& operator=(const TpsMonomial&) = default;

        /// Get exponent.
        ///  Return a reference to the exponent of variable `index`.
        int& operator[](int index);

        /// Get exponent.
        ///  Return value of the exponent of variable `index`.
        int operator[](int index) const;

        /// Product.
        ///  Return the exponent set for the product of the monomials denoted by
        ///  `this` and `rhs`, i.e. the sum of the exponents for each
        ///  variable.
        TpsMonomial operator*(const TpsMonomial& rhs) const;

        /// Convert.
        ///  Convert the monomial to the corresponding Giorgilli index.
        int getIndex() const;

        /// Get order.
        ///  Compute the monomial's order, i.e. the sum of all exponents.
        int getOrder() const;

        /// Get variables.
        ///  Return the monomial's number of variables, i.e. the number of exponents.
        int getVariables() const;

    private:
        // The monomial's exponents.
        Array1D<int> array;
    };

}  // namespace tps

#endif  // TPS_TpsMonomial_HH
