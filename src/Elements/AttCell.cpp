// ------------------------------------------------------------------------
// $RCSfile: AttCell.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AttCell
//   The representation of a table cell for OPAL ATTLIST commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/28 11:50:57 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Elements/AttCell.h"
#include <iostream>


//  Class AttCell.
// ------------------------------------------------------------------------

AttCell::AttCell()
{}


AttCell::~AttCell()
{}


void AttCell::setReal(double)
{}


void AttCell::setString(const std::string &)
{}


//  Class AttString.
// ------------------------------------------------------------------------

AttString::AttString():
    itsValue()
{}


AttString::~AttString()
{}


void AttString::clearValue() {
    itsValue.erase();
}


void AttString::printFormat(std::ostream &os) const {
    os << "%s";
}


void AttString::printValue(std::ostream &os) const {
    if(itsValue.empty()) {
        os << "Null";
    } else {
        os << itsValue;
    }
}


void AttString::setString(const std::string &value) {
    itsValue = value;
}


//  Class AttReal.
// ------------------------------------------------------------------------

AttReal::AttReal():
    itsValue(0.0)
{}


AttReal::~AttReal()
{}


void AttReal::clearValue() {
    itsValue = 0.0;
}


void AttReal::printFormat(std::ostream &os) const {
    os << "%le";
}


void AttReal::printValue(std::ostream &os) const {
    os << itsValue;
}


void AttReal::setReal(double value) {
    itsValue = value;
}
