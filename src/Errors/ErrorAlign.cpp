// ------------------------------------------------------------------------
// $RCSfile: ErrorAlign.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ErrorAlign
//   Class for OPAL EALIGN error command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/14 07:07:44 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Errors/ErrorAlign.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbstractObjects/BeamSequence.h"
#include "Attributes/Attributes.h"
#include "Errors/AlignHandler.h"
#include "Errors/Error.h"
#include "Utilities/Options.h"
#include <cmath>
#include <iostream>

using std::cerr;
using std::endl;


// Class ErrorAlign
// ------------------------------------------------------------------------

// The attributes of class ErrorAlign.
namespace {
    enum {
        DELTAX,     // The horizontal displacement in m.
        DELTAY,     // The vertical displacement in m.
        DELTAS,     // The longitudinal displacement in m.
        DTHETA, // The rotation around y axis in rad.
        DPHI,   // The rotation around x axis in rad.
        DPSI,   // The rotation around s axis in rad.
        MREX,   // The monitor x read error.
        MREY,   // The monitor y read error.
        MREDX,  // The monitor D_x read error.
        MREDY,  // The monitor D_y read error.
        SIZE
    };
}


ErrorAlign::ErrorAlign():
    Action(SIZE, "EALIGN",
           "The \"EALIGN\" sub-command assigns misalignments to the selected "
           "elements.") {
    itsAttr[DELTAX] = Attributes::makeReal
                      ("DELTAX", "Horizontal displacement in m");
    itsAttr[DELTAX].setDeferred(true);

    itsAttr[DELTAY] = Attributes::makeReal
                      ("DELTAY", "Vertical displacement in m");
    itsAttr[DELTAY].setDeferred(true);

    itsAttr[DELTAS] = Attributes::makeReal
                      ("DELTAS", "Longitudinal displacement in m");
    itsAttr[DELTAS].setDeferred(true);

    itsAttr[DTHETA] = Attributes::makeReal
                      ("DTHETA", "Rotation around y axis in rad");
    itsAttr[DTHETA].setDeferred(true);

    itsAttr[DPHI] = Attributes::makeReal
                    ("DPHI", "Rotation around x axis in rad");
    itsAttr[DPHI].setDeferred(true);

    itsAttr[DPSI] = Attributes::makeReal
                    ("DPSI", "Rotation around s axis in rad");
    itsAttr[DPSI].setDeferred(true);

    itsAttr[MREX] = Attributes::makeReal
                    ("MREX", "Monitor x read error in m");
    itsAttr[MREX].setDeferred(true);

    itsAttr[MREY] = Attributes::makeReal
                    ("MREY", "Monitor y read error in m");
    itsAttr[MREY].setDeferred(true);

    itsAttr[MREDX] = Attributes::makeReal
                     ("MREDX", "Monitor D_x read error in m");
    itsAttr[MREDX].setDeferred(true);

    itsAttr[MREDY] = Attributes::makeReal
                     ("MREDY", "Monitor D_y read error in m");
    itsAttr[MREDY].setDeferred(true);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


ErrorAlign::ErrorAlign(const std::string &name, ErrorAlign *parent):
    Action(name, parent)
{}


ErrorAlign::~ErrorAlign()
{}


ErrorAlign *ErrorAlign::clone(const std::string &name) {
    return new ErrorAlign(name, this);
}


void ErrorAlign::execute() {
    oldError = newError = 0;
    AlignHandler inserter(*Error::block->itsLine, *this, false);
    inserter.execute();
    cerr << '\n';

    if(Options::info) {
        if(oldError + newError == 0) {
            cerr << "No alignment errors assigned.\n";
        } else {
            if(oldError != 0) {
                cerr << "Alignment errors "
                     << (Error::block->addError ? "superposed on " : "replaced on ")
                     << oldError << " element";
                if(oldError > 1) cerr << "s";
                cerr << ".\n";
            }

            if(newError != 0) {
                cerr << "Alignment errors assigned to " << newError << " element";
                if(newError > 1) cerr << "s";
                cerr << ".\n";
            }
        }
    }

    cerr << endl;
}


void ErrorAlign::misalignment(const AlignWrapper &wrap, int) {
    double dx = Attributes::getReal(itsAttr[DELTAX]);
    double dy = Attributes::getReal(itsAttr[DELTAY]);
    double ds = Attributes::getReal(itsAttr[DELTAS]);
    double vx = Attributes::getReal(itsAttr[DTHETA]);
    double vy = Attributes::getReal(itsAttr[DPHI]);
    double vz = Attributes::getReal(itsAttr[DPSI]);
    Euclid3D newOffset(dx, dy, ds, vx, vy, vz);
    Euclid3D oldOffset = wrap.offset();

    if(oldOffset.isIdentity()) {
        newError++;
    } else {
        oldError++;
        if(Error::block->addError) newOffset *= oldOffset;
    }

    wrap.offset() = newOffset;
}