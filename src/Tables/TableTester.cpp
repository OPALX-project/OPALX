// ------------------------------------------------------------------------
// $RCSfile: TableTester.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TableTester
//   Ancillary class for testing for dependency of a Table on a name.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:45 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Tables/TableTester.h"
#include "AbsBeamline/ElementBase.h"


// Class TableTester
// ------------------------------------------------------------------------

TableTester::TableTester(const Beamline &beamline, const std::string &name):
    DefaultVisitor(beamline, false, false), itsName(name)
{}


TableTester::~TableTester()
{}


void TableTester::applyDefault(const ElementBase &elem) {
    if(elem.getName() == itsName) throw 0;
}
