#ifndef TPS_Tps_CC
#define TPS_Tps_CC 1
#include "Utilities/CLRangeError.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include "Utilities/LogicalError.h"
// Copyright: see LICENSE

// Template Class: Tps<class T>
//   Truncated power series in n variables of type T.
//   This file contains the template implementation only.

#include <algorithm>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <new>
#include "tps/Algebra/Array1D.h"
#include "tps/Algebra/Matrix.h"
#include "tps/Algebra/Tps.h"
#include "tps/Algebra/TpsData.h"
#include "tps/Algebra/TpsMonomial.h"
#include "tps/Algebra/TpsSubstitution.h"
#include "tps/Algebra/VpsMap.h"

namespace tps {

    // Helper functions.

    template <class T>
    inline void Tps<T>::setMaxOrder(int order) {
        rep->maxOrd = order;
    }

    template <class T>
    inline void Tps<T>::unique() {
        if (rep->ref > 1) {
            auto* copy = rep->clone();
            rep->ref--;
            rep = copy;
        }
    }

    // Template class TpsRep<T>.

    // Coefficients own their storage and construct/destruct each T correctly.
    // Reference-counted copy-on-write avoids coefficient copies on reads.

    template <class T>
    class TpsRep {
        friend class Tps<T>;

        // Memory management.
        explicit TpsRep(int size) : coefficients(size, T(0)) {}
        Array1D<T> coefficients;
        static TpsRep<T>* create(int maxOrder, int trcOrder, int variables);
        static TpsRep<T>* zero();
        static void release(TpsRep<T>*);

        // Make a copy of the representation.
        TpsRep<T>* clone();

        // Grab a new reference.
        TpsRep<T>* grab();

        // Return the monomial array.
        T* data();

        // The reference count.
        int ref;

        // Order and size of this object.
        int maxOrd;
        int trcOrd;

        // The length of the monomial table.
        int len;

        // The data structure for bookkeeping.
        TpsData* help;

        // Assignment/copy is deliberately disabled.
        TpsRep<T>& operator=(const TpsRep<T>&) = delete;
    };

    template <class T>
    inline T* TpsRep<T>::data() {
        return coefficients.begin();
    }

    template <class T>
    inline TpsRep<T>* TpsRep<T>::create(int maxOrder, int trcOrder, int variables) {
        TpsData::validateVariables(variables);
        if (maxOrder < 0 || maxOrder > INT_MAX - 2 || trcOrder < 0)
            throw ::CLRangeError("TpsRep::create()", "Invalid order.");
        // Construct descriptor and size.
        TpsData* d = 0;
        int s      = 1;
        if (variables) {
            d = TpsData::getTpsData(maxOrder, variables);
            s = d->getSize(maxOrder);
        }

        // Allocate representation and fill in data.
        auto owned = std::unique_ptr<TpsRep<T>>(new TpsRep<T>(s));
        auto* p    = owned.get();
        p->ref     = 1;
        p->maxOrd  = maxOrder;
        p->trcOrd  = trcOrder;
        p->len     = s;
        p->help    = d;

        return owned.release();
    }

    template <class T>
    inline TpsRep<T>* TpsRep<T>::zero() {
        // Allocate representation and fill in data.
        auto owned = std::unique_ptr<TpsRep<T>>(new TpsRep<T>(1));
        auto* p    = owned.get();
        p->ref     = 1;
        p->maxOrd  = 0;
        p->trcOrd  = Tps<T>::EXACT;
        p->len     = 1;
        p->help    = 0;

        return owned.release();
    }

    template <class T>
    inline TpsRep<T>* TpsRep<T>::clone() {
        // Allocate copy and copy monomial coefficients.
        auto owned      = std::unique_ptr<TpsRep<T>>(new TpsRep<T>(len));
        auto* p         = owned.get();
        p->coefficients = coefficients;

        // Copy limits and descriptor.
        p->ref    = 1;
        p->maxOrd = maxOrd;
        p->trcOrd = trcOrd;
        p->len    = len;
        p->help   = help;
        return owned.release();
    }

    template <class T>
    inline TpsRep<T>* TpsRep<T>::grab() {
        ++ref;
        return this;
    }

    template <class T>
    inline void TpsRep<T>::release(TpsRep<T>* p) {
        if (--(p->ref) <= 0) delete p;
    }

    // Template class Tps<T>.

    template <class T>
    Tps<T>::Tps() : rep(TpsRep<T>::zero()) {}

    template <class T>
    Tps<T>::Tps(int maxOrder, int nVar) : rep(TpsRep<T>::create(maxOrder, EXACT, nVar)) {}

    template <class T>
    Tps<T>::Tps(const Tps<T>& rhs) : rep(rhs.rep->grab()) {}

    template <class T>
    Tps<T>::Tps(const T& rhs) : Tps() {
        rep->data()[0] = rhs;
    }

    template <class T>
    Tps<T>::Tps(int rhs) : Tps(T(rhs)) {}

    template <class T>
    Tps<T>::~Tps() {
        TpsRep<T>::release(rep);
    }

    template <class T>
    Tps<T>& Tps<T>::operator=(const Tps<T>& rhs) {
        if (rep != rhs.rep) {
            TpsRep<T>::release(rep);
            rep = rhs.rep->grab();
        }

        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator=(const T& rhs) {
        auto replacement       = std::unique_ptr<TpsRep<T>>(TpsRep<T>::zero());
        replacement->data()[0] = rhs;
        TpsRep<T>::release(rep);
        rep = replacement.release();
        return *this;
    }

    template <class T>
    Tps<T> Tps<T>::filter(int minOrder, int maxOrder) const {
        // Compute order limits.
        maxOrder      = std::min(maxOrder, getMaxOrder());
        int trcOrder  = getTruncOrder();
        int variables = getVariables();

        // Construct filtered TpsRep.
        TpsRep<T>* p = TpsRep<T>::create(maxOrder, trcOrder, variables);
        int bot      = getSize(minOrder - 1);
        int top      = getSize(maxOrder);
        std::copy(rep->data() + bot, rep->data() + top, p->data() + bot);
        return Tps<T>(p);
    }

    template <class T>
    Tps<T> Tps<T>::truncate(int trunc) {
        return filter(0, trunc);
    }

    template <class T>
    const T Tps<T>::operator[](int index) const {
        return getCoefficient(index);
    }

    template <class T>
    T& Tps<T>::operator[](int index) {
        if (index < 0 || index >= rep->len)
            throw ::CLRangeError("Tps::operator[]", "Coefficient index out of range.");
        unique();
        return rep->data()[index];
    }

    template <class T>
    const T Tps<T>::operator[](const TpsMonomial& monomial) const {
        return (*this)[checkedCoefficientIndex(monomial)];
    }

    template <class T>
    T& Tps<T>::operator[](const TpsMonomial& monomial) {
        return (*this)[checkedCoefficientIndex(monomial)];
    }

    template <class T>
    const T Tps<T>::getCoefficient(int index) const {
        if (index < 0 || index >= rep->len) {
            throw ::CLRangeError("Tps::getCoefficients()", "Monomial index out of range.");
        }

        return rep->data()[index];
    }

    template <class T>
    void Tps<T>::setCoefficient(int index, const T& value) {
        if (index < 0 || (!rep->help && index != 0))
            throw ::CLRangeError("Tps::setCoefficient()", "Index out of range.");
        int order = getOrder(index);
        if (order > rep->trcOrd) return;

        if (order > rep->maxOrd) {
            TpsRep<T>* p = TpsRep<T>::create(order, rep->trcOrd, getVariables());
            std::copy(rep->data(), rep->data() + rep->len, p->data());
            TpsRep<T>::release(rep);
            rep = p;
        } else {
            unique();
        }

        rep->data()[index] = value;
    }

    template <class T>
    int Tps<T>::checkedCoefficientIndex(const TpsMonomial& monomial) const {
        if (monomial.getVariables() != getVariables())
            throw ::LogicalError("Tps::coefficientIndex()", "Inconsistent number of variables.");
        if (monomial.getOrder() > getMaxOrder())
            throw ::CLRangeError("Tps::coefficientIndex()", "Monomial exceeds stored order.");
        return monomial.getIndex();
    }

    template <class T>
    const T Tps<T>::getCoefficient(const TpsMonomial& monomial) const {
        return getCoefficient(checkedCoefficientIndex(monomial));
    }

    template <class T>
    void Tps<T>::setCoefficient(const TpsMonomial& monomial, const T& value) {
        int v1 = monomial.getVariables();
        int v2 = getVariables();

        if (v1 != v2) {
            throw ::LogicalError("Tps::getCoefficients()", "Inconsistent number of variables.");
        }

        setCoefficient(monomial.getIndex(), value);
    }

    template <class T>
    Tps<T> Tps<T>::makeVariable(int nVar, int var) {
        TpsData::validateVariables(nVar);
        if (var < 0 || var >= nVar)
            throw ::CLRangeError("Tps::makeVariable()", "Variable index out of range.");
        TpsRep<T>* p       = TpsRep<T>::create(1, EXACT, nVar);
        p->data()[var + 1] = T(1);
        return Tps<T>(p);
    }

    template <class T>
    Tps<T> Tps<T>::makeVarPower(int nVar, int var, int order) {
        TpsMonomial monomial(nVar);
        monomial[var] = order;
        Tps<T> z      = Tps<T>(TpsRep<T>::create(order, EXACT, nVar));
        z[monomial]   = T(1);
        return z;
    }

    template <class T>
    Tps<T> Tps<T>::makeMonomial(const TpsMonomial& monomial, const T& cc) {
        int order   = monomial.getOrder();
        int nVar    = monomial.getVariables();
        Tps<T> z    = Tps<T>(TpsRep<T>::create(order, EXACT, nVar));
        z[monomial] = cc;
        return z;
    }

    template <class T>
    Tps<T> Tps<T>::operator+() const {
        return *this;
    }

    template <class T>
    Tps<T> Tps<T>::operator-() const {
        TpsRep<T>* p = TpsRep<T>::create(getMaxOrder(), getTruncOrder(), getVariables());
        std::transform(rep->data(), rep->data() + rep->len, p->data(), std::negate<T>());
        return Tps<T>(p);
    }

    template <class T>
    Tps<T>& Tps<T>::operator+=(const Tps<T>& rhs) {
        if (int v1 = getVariables()) {
            if (int v2 = rhs.getVariables()) {
                if (v1 != v2) {
                    throw ::LogicalError(
                            "TpsRep::operator+=()", "Number of variables inconsistent.");
                }

                int trunc    = std::min(getTruncOrder(), rhs.getTruncOrder());
                int xOrder   = std::min(getMaxOrder(), trunc);
                int yOrder   = std::min(rhs.getMaxOrder(), trunc);
                int xLength  = getSize(xOrder);
                int yLength  = getSize(yOrder);
                int xyLength = std::min(xLength, yLength);
                TpsRep<T>* p = TpsRep<T>::create(std::max(xOrder, yOrder), trunc, v1);
                const T* x   = rep->data();
                const T* y   = rhs.rep->data();
                T* z         = p->data();
                std::transform(x, x + xyLength, y, z, std::plus<T>());
                std::copy(x + xyLength, x + xLength, z + xyLength);
                std::copy(y + xyLength, y + yLength, z + xyLength);
                TpsRep<T>::release(rep);
                rep = p;
            } else {
                unique();
                rep->data()[0] += rhs[0];
            }
        } else {
            *this = rhs + rep->data()[0];
        }
        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator-=(const Tps<T>& rhs) {
        if (int v1 = getVariables()) {
            if (int v2 = rhs.getVariables()) {
                if (v1 != v2) {
                    throw ::LogicalError(
                            "TpsRep::operator-=()", "Number of variables inconsistent.");
                }

                int trunc    = std::min(getTruncOrder(), rhs.getTruncOrder());
                int xOrder   = std::min(getMaxOrder(), trunc);
                int yOrder   = std::min(rhs.getMaxOrder(), trunc);
                int xLength  = getSize(xOrder);
                int yLength  = getSize(yOrder);
                int xyLength = std::min(xLength, yLength);

                TpsRep<T>* p = TpsRep<T>::create(std::max(xOrder, yOrder), trunc, v1);
                const T* x   = rep->data();
                const T* y   = rhs.rep->data();
                T* z         = p->data();
                std::transform(x, x + xyLength, y, z, std::minus<T>());
                std::copy(x + xyLength, x + xLength, z + xyLength);
                std::transform(y + xyLength, y + yLength, z + xyLength, std::negate<T>());
                TpsRep<T>::release(rep);
                rep = p;
            } else {
                unique();
                rep->data()[0] -= rhs[0];
            }
        } else {
            *this = -rhs + rep->data()[0];
        }
        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator*=(const Tps<T>& rhs) {
        return *this = multiply(rhs, truncOrder);
    }

    template <class T>
    Tps<T>& Tps<T>::operator/=(const Tps<T>& rhs) {
        return *this = multiply(rhs.inverse(truncOrder), truncOrder);
    }

    template <class T>
    Tps<T>& Tps<T>::operator+=(const T& rhs) {
        unique();
        rep->data()[0] += rhs;
        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator-=(const T& rhs) {
        unique();
        rep->data()[0] -= rhs;
        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator*=(const T& rhs) {
        unique();
        T* x = rep->data();
        std::transform(x, x + rep->len, x, [&rhs](const T& value) {
            return value * rhs;
        });
        return *this;
    }

    template <class T>
    Tps<T>& Tps<T>::operator/=(const T& rhs) {
        if (rhs == T(0)) throw ::DivideError("Tps::operator/()");
        T* x = rep->data();
        std::transform(x, x + rep->len, x, [&rhs](const T& value) {
            return value / rhs;
        });
        return *this;
    }

    template <class T>
    bool Tps<T>::operator==(const Tps<T>& rhs) const {
        if (int v1 = getVariables()) {
            if (int v2 = rhs.getVariables()) {
                if (v1 == v2) {
                    int trunc    = std::min(getTruncOrder(), rhs.getTruncOrder());
                    int xOrder   = std::min(getMaxOrder(), trunc);
                    int yOrder   = std::min(rhs.getMaxOrder(), trunc);
                    int xLength  = getSize(xOrder);
                    int yLength  = getSize(yOrder);
                    int xyLength = getSize(std::min(xOrder, yOrder));
                    const T* x   = rep->data();
                    const T* y   = rhs.rep->data();

                    for (int i = 0; i < xyLength; i++) {
                        if (x[i] != y[i]) return false;
                    }

                    for (int i = xyLength; i < xLength; i++) {
                        if (x[i] != T(0)) return false;
                    }

                    for (int i = xyLength; i < yLength; i++) {
                        if (y[i] != T(0)) return false;
                    }

                    return true;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        } else {
            if (rhs.getVariables()) {
                return false;
            } else {
                return rep->data()[0] == rhs.rep->data()[0];
            }
        }
    }

    template <class T>
    bool Tps<T>::operator==(const T& rhs) const {
        const T* x = rep->data();

        if (x[0] != rhs) return false;

        for (int i = 1; i < getSize(); ++i) {
            if (x[i] != T(0)) return false;
        }

        return true;
    }

    template <class T>
    bool Tps<T>::operator!=(const Tps<T>& rhs) const {
        return !(*this == rhs);
    }

    template <class T>
    bool Tps<T>::operator!=(const T& rhs) const {
        return !(*this == rhs);
    }

    template <class T>
    Tps<T> Tps<T>::substitute(const Matrix<T>& M) const {
        if (getVariables()) {
            int v1 = getVariables();
            int v2 = M.nrows();
            if (v1 != v2) {
                throw ::LogicalError("Tps::substitute()", "Matrix not consistent with Tps.");
            }
            int nRow = M.nrows();
            int nCol = M.ncols();

            // Define the nRow linear transformations.
            Array1D<Tps<T>> y(nRow);

            for (int i = 0; i < nRow; ++i) {
                y[i] = Tps<T>(1, nCol);
                for (int j = 0; j < nCol; ++j)
                    y[i][j + 1] = M[i][j];
            }

            // Evaluate the substitution.
            const T* x = rep->data();
            Tps<T> z(x[0]);

            if (int maxOrd = getMaxOrder()) {
                // Own the traversal: coefficient arithmetic can grow/rebuild the shared cache.
                const Array1D<TpsSubstitution> table = rep->help->getSubTable();
                Array1D<Tps<T>> product(maxOrd + 1);
                product[0] = Tps<T>(T(1));

                for (int next = 1; next < table.size();) {
                    const TpsSubstitution& s = table[next];
                    product[s.order]         = product[s.order - 1] * y[s.variable];
                    z += x[s.index] * product[s.order];
                    next = (s.order < maxOrd) ? next + 1 : s.skip;
                }
            }

            return z;
        } else {
            return *this;
        }
    }

    template <class T>
    Tps<T> Tps<T>::substitute(const VpsMap<T>& rhs) const {
        int v1 = getVariables();
        int v2 = rhs.getDimension();
        if (v1 != v2) {
            throw ::LogicalError("Tps::substitute()", "VpsMap is inconsistent with Tps.");
        }

        const T* x = rep->data();
        Tps<T> z(x[0]);

        if (int maxOrd = getMaxOrder()) {
            // Own the traversal: coefficient arithmetic can grow/rebuild the shared cache.
            const Array1D<TpsSubstitution> table = rep->help->getSubTable();
            Array1D<Tps<T>> product(maxOrd + 1);
            product[0] = Tps<T>(T(1));
            int trunc  = getTruncOrder();

            for (int next = 1; next < table.size();) {
                const TpsSubstitution& s = table[next];
                product[s.order]         = product[s.order - 1].multiply(rhs[s.variable], trunc);
                z += x[s.index] * product[s.order];
                next = (s.order < maxOrd) ? next + 1 : s.skip;
            }
        }

        return z;
    }

    template <class T>
    T Tps<T>::evaluate(const Vector<T>& rhs) const {
        int v1 = getVariables();
        int v2 = rhs.size();
        if (v1 != v2) {
            throw ::LogicalError("Tps::evaluate()", "Vector is inconsistent with Tps.");
        }

        const T* x = rep->data();
        T z        = x[0];

        if (int maxOrd = getMaxOrder()) {
            // Own the traversal: coefficient arithmetic can grow/rebuild the shared cache.
            const Array1D<TpsSubstitution> table = rep->help->getSubTable();
            Array1D<T> product(maxOrd + 1);
            product[0] = T(1);

            for (int next = 1; next < table.size();) {
                const TpsSubstitution& s = table[next];
                product[s.order]         = product[s.order - 1] * rhs[s.variable];
                z += x[s.index] * product[s.order];
                next = (s.order < maxOrd) ? next + 1 : s.skip;
            }
        }

        return z;
    }

    template <class T>
    void Tps<T>::clear() {
        auto* replacement = TpsRep<T>::zero();
        TpsRep<T>::release(rep);
        rep = replacement;
    }

    template <class T>
    std::istream& Tps<T>::get(std::istream& is) {
        is.flags(std::ios::skipws);
        char head[4];
        is.get(head, 4);
        if (strcmp(head, "Tps") != 0) {
            throw ::FormatError("Tps::get()", "Flag word \"Tps\" missing.");
        }

        int maxOrder, truncOrder, nVar;
        is >> maxOrder >> truncOrder >> nVar;
        Tps<T> z(TpsRep<T>::create(maxOrder, truncOrder, nVar));
        T coeff;

        if (nVar <= 0 || truncOrder == 0) {
            is >> coeff;
            z[0] = coeff;

            if (coeff != T(0)) {
                z[0] = coeff;
                is >> coeff;
            }
        } else {
            TpsMonomial monomial(nVar);
            maxOrder  = 0;
            bool done = false;
            bool fail = false;

            while (true) {
                is >> coeff;
                fail = is.fail();

                int order = 0;
                for (int var = 0; var < nVar; var++) {
                    int p;
                    is >> p;
                    fail |= is.fail();
                    if (p < 0) done = true;
                    monomial[var] = p;
                    order += monomial[var];
                }

                if (done) break;
                if (fail) throw ::FormatError("Tps::get()", "File read error");
                int index = monomial.getIndex();

                if (coeff != T(0)) {
                    maxOrder = order;
                } else if (index == 0) {
                    break;
                }

                z[index] = coeff;
            }

            z.setMaxOrder(maxOrder);
            *this = z;
        }

        return is;
    }

    template <class T>
    std::ostream& Tps<T>::put(std::ostream& os) const {
        std::streamsize old_prec = os.precision(14);
        os.setf(std::ios::scientific, std::ios::floatfield);

        int nVar = getVariables();
        os << "Tps " << getMaxOrder() << ' ' << getTruncOrder() << ' ' << nVar << std::endl;

        if (nVar == 0) {
            os << std::setw(24) << rep->data()[0] << std::endl;
        } else {
            for (int i = 0; i < getSize(); ++i) {
                if (rep->data()[i] != T(0)) {
                    os << std::setw(24) << rep->data()[i];

                    for (int var = 0; var < nVar; var++) {
                        os << std::setw(3) << getExponents(i)[var];
                    }

                    os << std::endl;
                }
            }

            os << std::setw(24) << T(0);

            for (int var = 0; var < nVar; var++) {
                os << std::setw(3) << (-1);
            }
        }

        os << std::endl;

        os.precision(old_prec);
        os.setf(std::ios::fixed, std::ios::floatfield);
        return os;
    }

    template <class T>
    Tps<T> Tps<T>::multiply(const Tps<T>& rhs, int trunc) const {
        int v1 = getVariables();
        int v2 = rhs.getVariables();
        if (v1) {
            if (v2) {
                if (v1 != v2) {
                    throw ::LogicalError("TpsRep::multiply()", "Number of variables inconsistent.");
                }

                if (getTruncOrder() != EXACT) {
                    int cut = getTruncOrder();
                    if (rhs[0] == 0.0) ++cut;
                    trunc = std::min(trunc, cut);
                }

                if (rhs.getTruncOrder() != EXACT) {
                    int cut = rhs.getTruncOrder();
                    if ((*this)[0] == 0.0) ++cut;
                    trunc = std::min(trunc, cut);
                }

                int maxOrder = std::min(getMaxOrder() + rhs.getMaxOrder(), trunc);

                TpsRep<T>* p = TpsRep<T>::create(maxOrder, trunc, v1);
                const T* x   = rep->data();
                T* z         = p->data();
                int yBot     = 0;
                int yHig     = std::min(rhs.getMaxOrder(), trunc);

                for (int yOrd = 0; yOrd <= yHig; yOrd++) {
                    int xOrd = std::min(getMaxOrder(), trunc - yOrd);
                    int xTop = getSize(xOrd);
                    int yTop = getSize(yOrd);

                    for (int yInd = yBot; yInd < yTop; yInd++) {
                        T y = rhs.rep->data()[yInd];
                        if (y != T(0)) {
                            const int* prod = rep->help->getProductArray(yInd);

                            for (int xInd = 0; xInd < xTop; xInd++) {
                                z[prod[xInd]] += x[xInd] * y;
                            }
                        }
                    }

                    yBot = yTop;
                }

                return Tps<T>(p);
            } else {
                Tps<T> result(TpsRep<T>::create(getMaxOrder(), getTruncOrder(), v1));
                const T* x = rep->data();
                const T y  = rhs.rep->data()[0];
                T* z       = result.rep->data();
                std::transform(x, x + rep->len, z, [&y](const T& value) {
                    return value * y;
                });
                return result;
            }
        } else {
            Tps<T> result(TpsRep<T>::create(rhs.getMaxOrder(), rhs.getTruncOrder(), v2));
            const T x  = rep->data()[0];
            const T* y = rhs.rep->data();
            T* z       = result.rep->data();
            std::transform(y, y + rhs.rep->len, z, [&x](const T& value) {
                return value * x;
            });
            return result;
        }
    }

    template <class T>
    Tps<T> Tps<T>::inverse(int trunc) const {
        if (trunc < 0) throw ::CLRangeError("Tps::inverse()", "Negative truncation order.");
        T aZero = rep->data()[0];
        if (aZero == T(0)) throw ::DivideError("Tps::inverse()");

        if (isConstant()) {
            return Tps<T>(T(1) / aZero);
        } else {
            int cut = getTaylorOrder(trunc);
            Array1D<T> series(cut + 1);
            series[0] = T(1) / aZero;

            for (int i = 1; i <= cut; i++) {
                series[i] = -series[i - 1] / aZero;
            }

            Tps<T> z = taylor(series);
            return z;
        }
    }

    template <class T>
    Tps<T> Tps<T>::derivative(int var) const {
        if (var < 0 || (getVariables() && var >= getVariables()))
            throw ::CLRangeError("Tps::derivative()", "Variable index out of range.");
        if (getVariables() && getMaxOrder() > 0) {
            int maxOrder = getMaxOrder() - 1;
            int trcOrder = getTruncOrder();

            TpsRep<T>* p = TpsRep<T>::create(maxOrder, trcOrder, getVariables());
            const T* x   = rep->data();
            T* z         = p->data();

            const int* product = rep->help->getProductArray(var + 1);

            for (int i = getSize(maxOrder); i-- > 0;) {
                int k = product[i];
                z[i]  = x[k] * double(getExponents(k)[var]);
            }

            Tps<T> result(p);
            return result;
        } else {
            Tps<T> result(TpsRep<T>::zero());
            return result;
        }
    }

    template <class T>
    Tps<T> Tps<T>::integral(int var) const {
        if (var < 0 || (getVariables() && var >= getVariables()))
            throw ::CLRangeError("Tps::integral()", "Variable index out of range.");
        if (getVariables() == 0) {
            throw ::LogicalError("TpsRep::integral()", "Cannot integrate a constant.");
        }

        int trcO     = std::min(rep->trcOrd == EXACT ? EXACT : rep->trcOrd + 1, truncOrder);
        int maxO     = std::min(rep->maxOrd + 1, trcO);
        TpsRep<T>* p = TpsRep<T>::create(maxO, trcO, getVariables());

        const T* x         = rep->data();
        T* z               = p->data();
        const int* product = rep->help->getProductArray(var + 1);

        for (int i = getSize(maxO - 1); i-- > 0;) {
            int k = product[i];
            z[k]  = x[i] / double(getExponents(k)[var]);
        }

        return Tps<T>(p);
    }

    template <class T>
    Tps<T> Tps<T>::multiplyVariable(int var) const {
        if (var < 0 || (getVariables() && var >= getVariables()))
            throw ::CLRangeError("Tps::multiplyVariable()", "Variable index out of range.");
        if (getVariables() == 0) {
            throw ::LogicalError(
                    "TpsRep::multiplyVariable()",
                    "Cannot multiply a constant by a numbered variable.");
        }

        int trcO     = std::min(rep->trcOrd == EXACT ? EXACT : rep->trcOrd + 1, truncOrder);
        int maxO     = std::min(rep->maxOrd + 1, trcO);
        TpsRep<T>* p = TpsRep<T>::create(maxO, trcO, getVariables());

        const T* x         = rep->data();
        T* z               = p->data();
        const int* product = rep->help->getProductArray(var + 1);

        for (int i = getSize(maxO - 1); i-- > 0;) {
            z[product[i]] = x[i];
        }

        return Tps<T>(p);
    }

    template <class T>
    Tps<T> Tps<T>::scaleMonomials(const Tps<T>& rhs) const {
        int v1 = getVariables();
        int v2 = rhs.getVariables();

        if (v1 != v2) {
            throw ::LogicalError("TpsRep::scaleMonomials()", "Number of variables inconsistent.");
        }

        int order = std::min(getMaxOrder(), rhs.getMaxOrder());
        int trunc = std::min(getTruncOrder(), rhs.getTruncOrder());

        TpsRep<T>* p = TpsRep<T>::create(std::min(order, trunc), trunc, v1);
        const T* x   = rep->data();
        const T* y   = rhs.rep->data();
        T* z         = p->data();
        std::transform(x, x + p->len, y, z, std::multiplies<T>());
        return Tps<T>(p);
    }

    template <class T>
    Tps<T> Tps<T>::taylor(std::span<const T> series) const {
        if (series.empty() || series.size() > static_cast<std::size_t>(INT_MAX))
            throw ::CLRangeError("Tps::taylor()", "Invalid coefficient span size.");
        const int order = static_cast<int>(series.size()) - 1;
        if (isConstant()) {
            return Tps<T>(series[0]);
        } else {
            Tps<T> x(*this);
            x[0] = T(0);
            Tps<T> z(series[order]);
            for (int maxOrder = 1; maxOrder <= order; maxOrder++) {
                z    = x.multiply(z, maxOrder);
                z[0] = series[order - maxOrder];
            }
            return z;
        }
    }

    template <class T>
    int Tps<T>::getMaxOrder() const {
        return rep->maxOrd;
    }

    template <class T>
    int Tps<T>::getTruncOrder() const {
        return std::min(rep->trcOrd, truncOrder);
    }

    template <class T>
    int Tps<T>::getVariables() const {
        return (rep->help != 0) ? rep->help->getVariables() : 0;
    }

    template <class T>
    int Tps<T>::getSize() const {
        return rep->len;
    }

    template <class T>
    bool Tps<T>::isConstant() const {
        return rep->help == 0;
    }

    template <class T>
    int Tps<T>::getTaylorOrder(int order) const {
        if (order < 0) throw ::CLRangeError("Tps::getTaylorOrder()", "Negative truncation order.");
        if (isConstant()) return 0;
        const int cut = std::min(order, getTruncOrder());
        return cut == EXACT ? getMaxOrder() : cut;
    }

    template <class T>
    int Tps<T>::getGlobalTruncOrder() {
        return truncOrder;
    }

    template <class T>
    void Tps<T>::setGlobalTruncOrder(int order) {
        if (order < 0)
            throw ::CLRangeError("Tps::setGlobalTruncOrder()", "Negative truncation order.");
        truncOrder = order;
    }

    template <class T>
    const TpsMonomial& Tps<T>::getExponents(int index) const {
        if (rep->help == 0) {
            throw ::LogicalError("Tps::getExponents()", "Cannot get exponents of a constant.");
        }

        return rep->help->getExponents(index);
    }

    template <class T>
    int Tps<T>::getOrder(int index) const {
        return rep->help ? rep->help->getOrder(index) : 0;
    }

    template <class T>
    int Tps<T>::getSize(int order) const {
        return rep->help ? rep->help->getSize(order) : 1;
    }

    template <class T>
    inline Tps<T>::Tps(TpsRep<T>* d) : rep(d) {}

}  // namespace tps

#endif  // TPS_Tps_CC
