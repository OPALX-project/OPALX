// ------------------------------------------------------------------------
// $RCSfile: OpalMultipole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalMultipole
//   The class of OPAL general multipoles.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalMultipole.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/DoomReader.h"
#include "AbstractObjects/DoomWriter.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/MultipoleRep.h"
#include "ComponentWrappers/MultipoleWrapper.h"
#include "Expressions/SValue.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/Options.h"
#include "Attributes/Attributes.h"  // JMJ added 20/12/2000
#include <iostream>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif
#include <vector>


// Class OpalMultipole
// ------------------------------------------------------------------------

OpalMultipole::OpalMultipole():
    OpalElement(SIZE, "MULTIPOLE",
                "The \"MULTIPOLE\" element defines a thick multipole.\n"
                "* If the length is non-zero, the strengths are per unit "
                "length.\n* If the length is zero, the strengths are the "
                "values integrated over the length.\n"
                "* With zero length no synchrotron radiation can be calculated.") {
    itsAttr[KN] = Attributes::makeRealArray
                  ("KN", "Normalised multipole strengths (normal) in m^(-k)");
    itsAttr[KS] = Attributes::makeRealArray
                  ("KS", "Normalised multipole strengths (skew) in m^(-k)");

    setElement((new MultipoleRep("MULTIPOLE"))->makeWrappers());
}


OpalMultipole::OpalMultipole(const string &name, OpalMultipole *parent):
    OpalElement(name, parent) {
    setElement((new MultipoleRep(name))->makeWrappers());
}


OpalMultipole::~OpalMultipole()
{}


OpalMultipole *OpalMultipole::clone(const string &name) {
    return new OpalMultipole(name, this);
}


void OpalMultipole::doomGet(const DoomReader &reader) {
    // Read the multipole length.
    itsAttr[LENGTH].doomGet(reader, DoomDB::getAttributeIndex("L"));

    // Read the multipole components.
    int index = DoomDB::getAttributeIndex("K0");
    int size = (reader.getRealSize() - index + 1) / 2;
    std::vector<double> norm(size, 0.0);
    std::vector<double> skew(size, 0.0);

    for(int i = 0; i < size; ++i) {
        norm[i] = reader.getReal(index++);
        skew[i] = reader.getReal(index++);
    }

    Attributes::setRealArray(itsAttr[KN], norm);
    Attributes::setRealArray(itsAttr[KS], skew);
}


void OpalMultipole::doomPut(DoomWriter &writer) const {
    // Set the object type name.
    writer.setTypeName("ELEMENT");

    // Save the multipole length.
    itsAttr[LENGTH].doomPut(writer, DoomDB::getAttributeIndex("L"));

    // Write the multipole components.
    const std::vector<double> norm = Attributes::getRealArray(itsAttr[KN]);
    const std::vector<double> skew = Attributes::getRealArray(itsAttr[KS]);
    int normSize = norm.size();
    int skewSize = skew.size();
    int size = normSize > skewSize ? normSize : skewSize;
    int index = DoomDB::getAttributeIndex("K0");

    for(int i = 0; i < size; ++i) {
        if(i < normSize) writer.putReal(index, norm[i]);
        index++;
        if(i < skewSize) writer.putReal(index, skew[i]);
        index++;
    }

    // Store the deflection angle and curvature.
    double length = Attributes::getReal(itsAttr[LENGTH]);
    double k0 = (normSize > 0) ? norm[0] : 0.0;
    index = DoomDB::getAttributeIndex("ANGLE");
    writer.putReal(index, k0 * length);
    writer.putInt(index, 1);

    index = DoomDB::getAttributeIndex("RHOINV");
    writer.putReal(index, k0);
    writer.putInt(index, 1);
}


void OpalMultipole::print(std::ostream &os) const {
        OpalElement::print(os);
}


void OpalMultipole::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
    const MultipoleWrapper *mult =
        dynamic_cast<const MultipoleWrapper *>(base.removeAlignWrapper());
    BMultipoleField field;

    // Get the desired field.
    if(flag == ERROR_FLAG) {
        field = mult->errorField();
    } else if(flag == ACTUAL_FLAG) {
        field = mult->getField();
    } else if(flag == IDEAL_FLAG) {
        field = mult->getDesign().getField();
    }

    double length = getLength();
    double scale = Physics::c / OPAL.getP0();
    if(length != 0.0) scale *= length;

    for(int order = 1; order <= field.order(); ++order) {
#if defined(__GNUC__) && __GNUC__ < 3
        char buffer[10];
        std::ostrstream ss(buffer, 10);
#else
        std::ostringstream ss;
#endif
        ss << (order - 1) << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
        string orderString(buffer);
#else
        std::string orderString = ss.str();
#endif

        string normName = "K" + orderString + "L";
        registerRealAttribute(normName)->setReal(scale * field.normal(order));

        string skewName = "K" + orderString + "SL";
        registerRealAttribute(skewName)->setReal(scale * field.skew(order));

        scale *= double(order);
    }
}


void OpalMultipole::update() {
    // Magnet length.
    MultipoleRep *mult =
        dynamic_cast<MultipoleRep *>(getElement()->removeWrappers());
    double length = getLength();
    mult->setElementLength(length);

    // Field components.
    BMultipoleField field;

    const std::vector<double> norm = Attributes::getRealArray(itsAttr[KN]);
    const std::vector<double> skew = Attributes::getRealArray(itsAttr[KS]);
    int normSize = norm.size();
    int skewSize = skew.size();
    double factor = OPAL.getP0() / Physics::c;
    int top = (normSize > skewSize) ? normSize : skewSize;

    for(int comp = 1; comp <= top; comp++) {
        factor /= double(comp);
        if(comp <= normSize) {
            field.setNormalComponent(comp, norm[comp-1] * factor);
	    mult->setNormalComponent(comp, norm[comp-1]);
        }
        if(comp <= skewSize) {
            field.setSkewComponent(comp, skew[comp-1] * factor);
	    mult->setSkewComponent(comp,skew[comp-1]);
        }
    }

    mult->setField(field);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(mult);
}
