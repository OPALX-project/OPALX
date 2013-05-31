// ------------------------------------------------------------------------
// $RCSfile: AbstractVar.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractVar
//   This abstract class defines the interface for a variable to be
//   adjusted during matching.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/AbstractVar.h"


// Class AbstractVar
// ------------------------------------------------------------------------


AbstractVar::AbstractVar(const string &name):
    itsName(name)
{}


AbstractVar::~AbstractVar()
{}


const string &AbstractVar::getName() const {
    return itsName;
}

