// ------------------------------------------------------------------------
// $RCSfile: SBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/SBend.h"
#include "Algorithms/PartPusher.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
    Bend()
{
    setMessageHeader("SBend ");
}

SBend::SBend(const SBend &right):
    Bend(right)
{
    setMessageHeader("SBend ");
}

SBend::SBend(const std::string &name):
    Bend(name)
{
    setMessageHeader("SBend ");
}

SBend::~SBend() {
}

void SBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitSBend(*this);
}

/*
 * OPAL-MAP methods
 * ================
 */
double SBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}

double SBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}

void SBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}

void SBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}


/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */

ElementBase::ElementType SBend::getType() const {
    return SBEND;
}


bool SBend::FindChordLength(Inform &msg,
                            double &chordLength,
                            bool &chordLengthFromMap) {

    /*
     * Find bend chord length. If this was not set by the user using the
     * L (length) attribute, infer it from the field map.
     */
    chordLength = getLength();
    if(chordLength > 0.0) {
        chordLengthFromMap = false;
        return true;
    } else {

        if(chordLength == 0.0)
            chordLength = getLength();

        chordLengthFromMap = true;

        if(chordLength <= 0.0) {
            ERRORMSG(level2 << "Magnet length inferred from field map is less than or equal"
                     << " to zero. Check your bend magnet input." << endl);
            return false;
        } else
            return true;

    }
}