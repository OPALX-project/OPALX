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
#include "AbstractObjects/DoomDB.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class AlignReader
// ------------------------------------------------------------------------


AlignReader::AlignReader(const string &name) {
    // Set DOOM environment for align error.
    DOOM_DB.setEnvironment("USE", name);
}


AlignReader::~AlignReader()
{}


void AlignReader::misalignment(const AlignWrapper &wrap, int occur) {
    Euclid3D offset;

    if(DOOM_DB.readAlign(wrap.getElement()->getName(), occur, offset)) {
        wrap.offset() = offset;
    }
}
