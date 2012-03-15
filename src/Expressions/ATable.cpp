// ------------------------------------------------------------------------
// $RCSfile: ATable.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ATable
//   Expression class. Generates a list of values from a TABLE() function.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/01/22 15:16:02 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Expressions/ATable.h"

#include <iostream>

using namespace std;

// Class ATable
// ------------------------------------------------------------------------

namespace Expressions {

    ATable::ATable(const ATable &rhs):
        Array<double>(rhs), itsExpr(rhs.itsExpr),
        itsBegin(rhs.itsBegin), itsEnd(rhs.itsEnd), itsStep(rhs.itsEnd)
    {}


    ATable::ATable(int n1, int n2, int n3):
        Array<double>(), itsExpr(), itsBegin(n1), itsEnd(n2), itsStep(n3)
    {}


    ATable::~ATable()
    {}


    Array<double> *ATable::clone() const {
        return new ATable(*this);
    }


    void ATable::defineExpression(PtrToScalar<double> expr) {
        itsExpr = expr;
    }


    std::vector<double> ATable::evaluate() const {
        std::vector<double> result(itsEnd, 0.0);

        for(itsHash = itsBegin; itsHash <= itsEnd; itsHash += itsStep) {
            result[itsHash-1] = itsExpr->evaluate();
        }

        return result;
    }


    double ATable::getHash() const {
        return itsHash;
    }


    std::ostream &ATable::print(std::ostream &os, int) const {
        os << "TABLE(" << itsBegin << ':' << itsEnd << ':' << itsStep << ',';
        itsExpr->print(os);
        return os << ')';
    }

}
