// ------------------------------------------------------------------------
// $RCSfile: Selector.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Selector
//   This class sets selected selection flags in a USE object.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:25:22 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Tables/Selector.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/Table.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "Elements/OpalElement.h"
#include "Utilities/Options.h"
#include "Utilities/RegularExpression.h"
#include <iostream>


// Class Selector
// ------------------------------------------------------------------------

Selector::Selector(const Beamline &bl, const RangeRep &range,
                   const string &clsName, const string &typName,
                   const string &pattern):
    RangeSelector(bl, range),
    itsClass(0), itsType(typName), itsPattern(0), itsCount(0) {
    if(! clsName.empty()  && (itsClass = Element::find(clsName)) == 0) {
        if(Options::warn) {
            std::cerr << "\n### Warning ### Unknown class name \""
                      << clsName << "\"; will select all classes.\n" << std::endl;
        }
    }

    if(! pattern.empty()) {
        itsPattern = new RegularExpression(pattern);
    }
}


Selector::~Selector() {
    if(itsPattern) delete itsPattern;
}


void Selector::execute() {
    itsCount = 0;
    RangeSelector::execute();
}


void Selector::handleElement(const FlaggedElmPtr &fep) {
    // Skip elements which are not in range.
    if(itsRange.isActive()) {
        const string &name = fep.getElement()->getName();
        if(name[0] != '[') {
            bool set = true;
            OpalElement &elem = dynamic_cast<OpalElement &>(*Element::find(name));

            // If class exists and element is not class member, then skip.
            if(itsClass != 0  &&  ! elem.isTreeMember(itsClass)) set = false;

            // If pattern does exists and element name does not match, then skip.
            if(itsPattern != 0  &&  ! itsPattern->match(name)) set = false;

            // If type name is not blank and element type is different, then skip.
            if(! itsType.empty() && itsType != elem.getTypeName()) set = false;

            // The current element matches all conditions.
            if(set) {
                fep.setSelectionFlag(true);
                ++itsCount;
            }
        }
    }
}


int Selector::getCount() const {
    return itsCount;
}
