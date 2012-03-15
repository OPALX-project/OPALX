// ------------------------------------------------------------------------
// $RCSfile: MPRemover.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPRemover
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:41 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/MPRemover.h"
#include "Fields/BMultipoleField.h"


// Class MPRemover
// ------------------------------------------------------------------------

MPRemover::MPRemover()
{}


MPRemover::~MPRemover()
{}


void MPRemover::fieldError(const string &, int,
                           const BMultipoleField &,
                           BMultipoleField &errorField) {
    errorField = BMultipoleField();
}
