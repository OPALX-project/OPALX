// ------------------------------------------------------------------------
// $RCSfile: ConcreteVar.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.5 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConcreteVar
//   This class is the concrete class for a single variable to be varied
//   during matching.  It implements the setting and retrieving of the
//   value in the system to be adjusted.
//
//   JMJ 7/4/2000:  changed to make output from MATCH like OPAL input.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/03/28 21:27:54 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Match/ConcreteVar.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>


// Class ConcreteVar
// ------------------------------------------------------------------------


ConcreteVar::ConcreteVar
(const string &name, Attribute &attr, int limits, double pars[4]):
    AbstractVar(name),
    itsAttr(attr),
    itsStep(std::abs(pars[1])),
    itsLimits(limits),
    itsMin(pars[2]),
    itsMax(pars[3]) {
    // Builtin constant.
    const double stplim = 0.2;

    // Convert the step sizes.
    switch(limits) {
        case UPPER_LIMIT:
        case LOWER_LIMIT:
            if(pars[0] != 0.0) itsStep /= pars[0];
            break;
        case BOTH_LIMITS:
            itsStep /= sqrt(pars[0] * (pars[0] - itsMin - itsMax) + itsMin * itsMax);
            break;
        default:
            break;
    }

    // Bring the step in range.
    if(itsStep > stplim) itsStep = stplim;
}


ConcreteVar::~ConcreteVar()
{}


double ConcreteVar::getInternalValue() const {
    double value = Attributes::getReal(itsAttr);

    switch(itsLimits) {
        case LOWER_LIMIT:
            if(value < itsMin) value = itsMin;
            return sqrt(value - itsMin);
        case UPPER_LIMIT:
            if(value > itsMax) value = itsMax;
            return sqrt(itsMax - value);
        case BOTH_LIMITS:
            if(value < itsMin) value = itsMin;
            if(value > itsMax) value = itsMax;
            return asin((value + value - itsMin - itsMax) / (itsMax - itsMin));
        default:
            return value;
    }
}


void ConcreteVar::setInternalValue(double value) {
    // Store internal value to chache.
    switch(itsLimits) {
        case LOWER_LIMIT:
            value = itsMin + value * value;
            break;
        case UPPER_LIMIT:
            value = itsMax - value * value;
            break;
        case BOTH_LIMITS:
            value = ((itsMin + itsMax) + (itsMax - itsMin) * sin(value)) / 2.0;
            break;
        default:
            break;
    }

    // Transmit to data structure.
    setExternalValue(value);
}


double ConcreteVar::getExternalValue() const {
    return Attributes::getReal(itsAttr);
}


void ConcreteVar::setExternalValue(double value) {
    // Transmit external value to cache and to data structure.
    Attributes::setReal(itsAttr, value);
    OPAL.makeDirty(0);
}


void ConcreteVar::print(std::ostream &os) const {
    std::streamsize old_prec = os.precision(10);
    os  << itsName << " := "
        << std::setw(16) << getExternalValue() << ";" << std::endl;
    os.precision(old_prec);
}
//  Above changed by JMJ, 7/4/2000, increased precision to 10, swapped columns
//  and slipped in the ":=".
//  Doesn't seem to be happening in version of July 2000 ??
