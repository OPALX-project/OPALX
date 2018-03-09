// ------------------------------------------------------------------------
// $RCSfile: ConcreteFun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConcreteFun
//   This class is the concrete class for a single constrained value
//   for matching.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/ConcreteFun.h"
#include "Attributes/Attributes.h"
#include "Utilities/LogicalError.h"
#include <iomanip>
#include <iostream>


// Character form of relational operators.
static const char *relString[] = { " == ", " > ", " < " };


// Class ConcreteFun
// ------------------------------------------------------------------------

ConcreteFun::ConcreteFun
(Attribute &lhs, int rel, Attribute &rhs, Attribute &wgt):
    itsLhs(lhs),
    itsRhs(rhs),
    relation(rel),
    itsWeight(Attributes::getRealArray(wgt))
{}


ConcreteFun::~ConcreteFun()
{}


int ConcreteFun::countConstraints() const {
    return itsWeight.size();
}


void ConcreteFun::evaluate(Vector<double> &fun, int &index) const {
    lValue = Attributes::getRealArray(itsLhs);
    rValue = Attributes::getRealArray(itsRhs);
    std::vector<double>::size_type size = itsWeight.size();

    if(lValue.size() != size || rValue.size() != size) {
        throw LogicalError("ConcreteFun::evaluate()",
                           "Matching constraint is inconsistent, size_of(LHS), "
                           "size_of(RHS), size_of(WGT) should all be equal.");
    }

    value.erase(value.begin(), value.end());
    value.resize(size, 0.0);

    for(std::vector<double>::size_type i = 0; i < size; ++i) {
        switch(relation) {

            case 0:
                value[i] = lValue[i] - rValue[i];
                break;

            case 1:
                if(lValue[i] < rValue[i]) value[i] = lValue[i] - rValue[i];
                break;

            case 2:
                if(lValue[i] > rValue[i]) value[i] = lValue[i] - rValue[i];
        }

        fun[index++] = value[i] * itsWeight[i];
    }
}


void ConcreteFun::print(std::ostream &os) const {
    // Print the symbolic representation of the constraint.
    os << "\nCondition: " << itsLhs << relString[relation] << itsRhs
       << std::endl;

    // Print the value(s) of the constraint.
    std::streamsize old_prec = os.precision(9);
    os.setf(std::ios::scientific, std::ios::floatfield);

    for(std::vector<double>::size_type i = 0; i < lValue.size(); ++i) {
        os << "  lhs = " << std::setw(16) << lValue[i]
           << ", rhs = " << std::setw(16) << rValue[i]
           << ", wgt = " << std::setw(16) << itsWeight[i]
           << ", val = " << std::setw(16) << value[i]
           << std::endl;
    }

    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(old_prec);
}
