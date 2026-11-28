#ifndef TPS_TpsData_H
#define TPS_TpsData_H
#include "Utilities/CLRangeError.h"
// Copyright: see LICENSE

// Class: TpsData

#include "tps/Algebra/Array1D.h"
#include "tps/Algebra/Array2D.h"
#include "tps/Algebra/TpsMonomial.h"
#include "tps/Algebra/TpsSubstitution.h"

namespace tps {

    // Class TpsData
    /** @brief Internal shared monomial indexing and substitution tables for Tps.
     * @details One cache exists per variable count. Growing its order replaces the
     * underlying tables and invalidates borrowed element references/product pointers.
     * Access and growth require external synchronization; rebuilding is not transactional
     * on allocation failure. Prefer the scalar/map APIs over direct cache access.
     */

    class TpsData {
    public:
        /// Fixed cache supports scalar (0) and 1..99 variables.
        static const int MAX_VARIABLES = 99;
        static int validateVariables(int variables);
        TpsData();
        ~TpsData();

        /// Obtain/grow shared tables through nOrd; zero variables returns nullptr.
        /// Invalid orders/counts and nonrepresentable table sizes throw ::CLRangeError.
        static TpsData* getTpsData(int nOrd, int nVar);

        /// Power array for monomial with index "index"
        const TpsMonomial& getExponents(int index) const;

        /// Order of monomial with index "index".
        int getOrder(int index) const;

        /// Borrow product indices for one monomial; valid until cache growth.
        /// The row covers getSize(built order - getOrder(index)) factors.
        const int* getProductArray(int index) const;

        /// Number of terms through total degree order; getSize(-1) returns zero.
        /// Orders outside [-1, built order] throw ::CLRangeError.
        int getSize(int order) const;

        /// Depth-first traversal records, with subtree skip indices.
        const Array1D<TpsSubstitution>& getSubTable() const;

        /// Number of variables
        int getVariables() const;

        /// Index for a monomial from its exponents
        int indexMonomial(const TpsMonomial&) const;

    private:
        // Assignment/copy is deliberately disabled.
        TpsData(const TpsData&)        = delete;
        void operator=(const TpsData&) = delete;

        // TpsData builder.
        void build(int max, int nVar);

        // Initialise TpsData to empty.
        void clear();

        // Initialise Substitution tables.
        void fillSubst(int var, int order, TpsMonomial& pow, int& next);

        // Number of variables.
        int variables;

        // Maximum order to which the data are expanded.
        int topOrder;

        // Maximum Tps size.
        int maxSize;

        // Exponents: expon[index] = set of exponents in monomial "index".
        typedef Array1D<TpsMonomial> ExponentTable;
        ExponentTable expon;

        // Indexing information, used for Tps multiplication.
        typedef Array1D<int> ProductRow;
        typedef Array1D<ProductRow> ProductTable;
        ProductTable prod;

        // Binomial coefficients: binom[i-1][k]=(i+k) ! / i! / k!
        // used in the indexMonomial() and getSize(int) methods.
        typedef Array2D<int> BinomialTable;
        BinomialTable binom;

        // The substitution table.
        Array1D<TpsSubstitution> subTable;
    };

    // Class TpsData, inline methods.

    inline const TpsMonomial& TpsData::getExponents(int index) const {
        if (index < 0 || index >= maxSize)
            throw ::CLRangeError("TpsData::getExponents()", "Index out of range.");
        return expon[index];
    }

    inline const int* TpsData::getProductArray(int index) const {
        if (index < 0 || index >= maxSize)
            throw ::CLRangeError("TpsData::getProductArray()", "Index out of range.");
        return prod[index].begin();
    }

    inline int TpsData::getOrder(int index) const { return getExponents(index).getOrder(); }

    inline int TpsData::getSize(int order) const {
        if (order < -1 || order > topOrder)
            throw ::CLRangeError("TpsData::getSize()", "Order outside built cache.");
        return order == -1 ? 0 : binom[0][order + 1];
    }

    inline const Array1D<TpsSubstitution>& TpsData::getSubTable() const { return subTable; }

    inline int TpsData::getVariables() const { return variables; }

}  // namespace tps

#endif  // TPS_TpsData_H
