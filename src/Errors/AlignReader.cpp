// ------------------------------------------------------------------------
// $RCSfile: AlignReader.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignReader
//   Ancillary class for reading align errors from DOOM data base.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignReader.h"
#include "AbsBeamline/AlignWrapper.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class AlignReader
// ------------------------------------------------------------------------


AlignReader::AlignReader(const string &name) {

}


AlignReader::~AlignReader()
{}


void AlignReader::misalignment(const AlignWrapper &wrap, int occur) {
    Euclid3D offset;

}
