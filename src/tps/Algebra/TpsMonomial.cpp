// Copyright: see Copyright.readme

// Class: TpsMonomial
//   Index set for a Tps monomial.

#include "tps/Algebra/TpsMonomial.h"
#include "tps/Algebra/TpsData.h"
#include <climits>

namespace tps {

// Class TpsMonomial

TpsMonomial::TpsMonomial(int nVar) : array(TpsData::validateVariables(nVar), int(0)) {}

TpsMonomial::TpsMonomial(int nVar, int var) : array(TpsData::validateVariables(nVar), int(0)) {
    (*this)[var] = 1;
}

int& TpsMonomial::operator[](int index) {
    if (index < 0 || index >= array.size())
        throw CLRangeError("TpsMonomial::operator[]", "Variable index out of range.");
    return array[index];
}

int TpsMonomial::operator[](int index) const {
    if (index < 0 || index >= array.size())
        throw CLRangeError("TpsMonomial::operator[]", "Variable index out of range.");
    return array[index];
}

TpsMonomial TpsMonomial::operator*(const TpsMonomial& rhs) const {
    int n = rhs.getVariables();
    if (n != getVariables())
        throw SizeError("TpsMonomial::operator*()", "Variable count mismatch.");
    const int leftOrder = getOrder(), rightOrder = rhs.getOrder();
    if (leftOrder > INT_MAX - rightOrder)
        throw CLRangeError("TpsMonomial::operator*()", "Order overflow.");
    TpsMonomial z(*this);
    for (int i = 0; i < n; i++)
        z.array[i] += rhs.array[i];
    return z;
}

int TpsMonomial::getIndex() const {
    // The empty monomial is the scalar constant and has no cache.
    if (array.size() == 0)
        return 0;
    TpsData* data = TpsData::getTpsData(getOrder(), array.size());
    return data->indexMonomial(*this);
}

int TpsMonomial::getOrder() const {
    int order = 0;

    for (int i = 0; i < array.size(); i++) {
        if (array[i] < 0 || array[i] > INT_MAX - order)
            throw CLRangeError("TpsMonomial::getOrder()", "Negative exponent or order overflow.");
        order += array[i];
    }

    return order;
}

int TpsMonomial::getVariables() const {
    return array.size();
}

} // namespace tps
