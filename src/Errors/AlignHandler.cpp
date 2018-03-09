// ------------------------------------------------------------------------
// $RCSfile: AlignHandler.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: AlignHandler
//   This Visitor class provides access to the misalignment errors.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignHandler.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/Marker.h"
#include "AbstractObjects/Table.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Errors/AlignBase.h"
#include "Utilities/Options.h"


// Class AlignHandler
// ------------------------------------------------------------------------

AlignHandler::AlignHandler(const Beamline &beamline,
                           AlignBase &cmd,
                           bool full):
    DefaultVisitor(beamline, false, false),
    command(cmd),
    sHandler(),
    applyToAll(full)
{}


AlignHandler::~AlignHandler()
{}


void AlignHandler::visitAlignWrapper(const AlignWrapper &wrap) {
    // Must be called before generating misalignments.
    DefaultVisitor::visitAlignWrapper(wrap);

    // Generate misalignment.
    if(applyToAll || isSelected) {
        command.misalignment(wrap, occurCount);
    }
}


void AlignHandler::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    // Set selection flag for use in the embedded element.
    isSelected = fep.getSelectionFlag();
    occurCount = fep.getCounter();

    // Forward call to element.
    DefaultVisitor::visitFlaggedElmPtr(fep);
}


void AlignHandler::applyDefault(const ElementBase &elem) {
    // Update the length for the S-function.
    sHandler.update(elem.getElementLength());
}
