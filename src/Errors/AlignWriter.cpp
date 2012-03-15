// ------------------------------------------------------------------------
// $RCSfile: AlignWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignWriter
//   Ancillary class for writing align errors to DOOM data base.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignWriter.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/Table.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class AlignWriter
// ------------------------------------------------------------------------


AlignWriter::AlignWriter(const string &name) {
    // Set DOOM environment for align error.
    DOOM_DB.setEnvironment("USE", name);
}


AlignWriter::~AlignWriter()
{}


void AlignWriter::misalignment(const AlignWrapper &wrap, int occur) {
    const Euclid3D &offset = wrap.offset();

    if(! offset.isIdentity()) {
        DOOM_DB.writeAlign(wrap.getElement()->getName(), occur, offset);
    }
}
