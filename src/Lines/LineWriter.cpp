// ------------------------------------------------------------------------
// $RCSfile: LineWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: LineWriter
//   The class for the OPAL LineWriter command.
//
// ------------------------------------------------------------------------
//
// Revision History:
// $Date: 2001/08/13 15:16:16 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Lines/LineWriter.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Algorithms/MPSplitIntegrator.h"
#include "Attributes/Attributes.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"
#include "doomex.h"
#include <vector>


// Class LineWriter.
// ------------------------------------------------------------------------
//: A visitor to write out a line or sequence definition.

LineWriter::LineWriter(const Beamline &bl, const string &oName,
                       const string &pName):
    DefaultVisitor(bl, false, false),
    itsWriter(), itsCount(0), itsPosition(0.0) {
    static int keyList[2] = { 1, SEQU_CODE };
    itsWriter.setObjectName(oName, keyList);
    itsWriter.setParentName(pName);
    itsWriter.setTypeName("EX_SEQUENCE");
}


LineWriter::~LineWriter()
{}


void LineWriter::visitDrift(const Drift &drf) {
    // All drifts are skipped, but their length is accumulated.
    itsPosition += drf.getElementLength();
}


void LineWriter::visitMapIntegrator(const MapIntegrator &i) {
    ElementBase *base = i.getElement();

    if(const MPSplitIntegrator *mpi =
           dynamic_cast<const MPSplitIntegrator *>(&i)) {
        // For a MPSliceIntegrator store the equivalent thin lenses.
        double length = base->getElementLength();
        std::vector<double> itsSlices;
        mpi->getSlices(itsSlices);
        int slices = itsSlices.size() - 1;

        if(slices > 1  &&  length != 0.0) {
            // Construct internal element name.
            const string name = base->getName();
            char extension[4] = "$00";
            extension[1] = char(slices / 10 + '0');
            extension[2] = char(slices % 10 + '0');
            const string splitName = name + extension;

            // Find or construct internal element.
            if(OPAL.find(splitName) == 0) {
                Object *split = Element::find(name)->clone(splitName);
                Attribute *attr = split->findAttribute("L");
                Attributes::setReal(*attr, length / double(slices));
                OPAL.define(split);
                DOOM_DB.writeObject(split);
            }

            // Handle the slices.
            for(int k = 0; k < slices; ++k) {
                itsWriter.putString(itsCount, splitName);
                itsWriter.putInt(itsCount, occurCount);
                itsWriter.putReal(itsCount, itsPosition + itsSlices[k] * length);
                ++itsCount;
            }

            itsPosition += length;
        } else {
            LineWriter::applyDefault(*base);
        }
    } else {
        LineWriter::applyDefault(*base);
    }
}


void LineWriter::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    // Save occurrence counter for use in the embedded element.
    occurCount = fep.getCounter();

    // Enter the embedded element.
    DefaultVisitor::visitFlaggedElmPtr(fep);
}


void LineWriter::applyDefault(const ElementBase &elem) {
    const string name = elem.getName();

    if(name[0] != '#') {
        double length = elem.getElementLength();
        itsWriter.putString(itsCount, name);
        itsWriter.putInt(itsCount, occurCount);
        itsWriter.putReal(itsCount, itsPosition + length / 2.0);
        itsPosition += length;
        ++itsCount;
    }
}
