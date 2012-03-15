// ------------------------------------------------------------------------
// $RCSfile: RangeSelector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RangeSelector
//   This abstract class runs through a beam line and calls the pure
//   virtual methods RangeSelector::handleXXX() for each element or
//   beamline in a range.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:35:16 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Tables/RangeSelector.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

class Element;


// Class RangeSelector
// ------------------------------------------------------------------------

RangeSelector::RangeSelector(const Beamline &beamline, const RangeRep &range):
    DefaultVisitor(beamline, false, false), itsRange(range)
{}


RangeSelector::~RangeSelector()
{}


void RangeSelector::execute() {
    itsRange.initialize();
    DefaultVisitor::execute();
}


void RangeSelector::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    // Are we past the beginning of the range ?
    itsRange.enter(fep);

    // Do the required operations on the beamline or element.
    ElementBase *base = fep.getElement()->removeWrappers();

    if(dynamic_cast<Beamline *>(base)) {
        handleBeamline(fep);
    } else {
        handleElement(fep);
    }

    // Are we past the end of the range ?
    itsRange.leave(fep);
}


void RangeSelector::handleBeamline(const FlaggedElmPtr &fep) {
    DefaultVisitor::visitFlaggedElmPtr(fep);
}


void RangeSelector::handleElement(const FlaggedElmPtr &fep) {
    // Default: delegate algorithm to the element, if in range.
    if(itsRange.isActive()) {
        DefaultVisitor::visitFlaggedElmPtr(fep);
    }
}
