// ------------------------------------------------------------------------
// $RCSfile: RealVector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RealVector
//   Implements a OPAL REAL_VECTOR definition.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:49 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "ValueDefinitions/RealVector.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include <iostream>
#include <vector>


// Class RealVector
// ------------------------------------------------------------------------

RealVector::RealVector():
    ValueDefinition(1, "REAL_VECTOR",
                    "The \"REAL VECTOR\" statement defines a global "
                    "real vector.\n"
                    "\tREAL VECTOR<name>=<real-vector-expression>;\n") {
    itsAttr[0] = Attributes::makeRealArray("VALUE", "The vector value");
}


RealVector::RealVector(const std::string &name, RealVector *parent):
    ValueDefinition(name, parent)
{}


RealVector::~RealVector()
{}


bool RealVector::canReplaceBy(Object *object) {
    // Replace only by another vector.
    return (dynamic_cast<RealVector *>(object) != 0);
}


RealVector *RealVector::clone(const std::string &name) {
    return new RealVector(name, this);
}


void RealVector::print(std::ostream &os) const {
    // WARNING: Cannot print in OPAL-8 format.
    os << "REAL VECTOR " << getOpalName() << ":="
       << itsAttr[0] << ';' << std::endl;
}

void RealVector::printValue(std::ostream &os) const {
    os << itsAttr[0];
}

double RealVector::getRealComponent(int index) const {
    std::vector<double> array = Attributes::getRealArray(itsAttr[0]);
    return array[index-1];
}