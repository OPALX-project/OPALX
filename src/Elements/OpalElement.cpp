// ------------------------------------------------------------------------
// $RCSfile: OpalElement.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalElement
//   The base class for all OPAL beamline elements.

//   and for printing in OPAL-8 format.
//
// ------------------------------------------------------------------------
//
// $Date: 2002/12/09 15:06:07 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"
#include "AbsBeamline/AlignWrapper.h"
#include "AbsBeamline/ElementImage.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Parser/Statement.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include "Utilities/Round.h"
#include <cmath>
#include <cctype>
#if defined(__GNUC__) && __GNUC__ < 3
#include <strstream>
#else
#include <sstream>
#endif
#include <vector>


// Class OpalElement
// ------------------------------------------------------------------------

std::map < std::string, OwnPtr<AttCell> > OpalElement::attributeRegistry;


OpalElement::~OpalElement()
{}


void OpalElement::
fillRegisteredAttributes(const ElementBase &base, ValueFlag) {
    // Fill in the common data for all elements.
    attributeRegistry["NAME"]->setString(getOpalName());
    attributeRegistry["TYPE"]->setString(getTypeName());
    attributeRegistry["CLASS"]->setString(getParent()->getOpalName());
    attributeRegistry["KEYWORD"]->setString(getBaseObject()->getOpalName());
    attributeRegistry["L"]->setReal(base.getElementLength());

    // Misalignments.
    const AlignWrapper *wrap = dynamic_cast<const AlignWrapper *>(&base);
    if(wrap) {
        double dx, dy, ds, dphi, dtheta, dpsi;
        wrap->offset().getAll(dx, dy, ds, dtheta, dphi, dpsi);
        attributeRegistry["DELTAX"]->setReal(dx);
        attributeRegistry["DELTAY"]->setReal(dy);
        attributeRegistry["DELTAS"]->setReal(ds);
        attributeRegistry["DTHETA"]->setReal(dtheta);
        attributeRegistry["DPHI"]->setReal(dphi);
        attributeRegistry["DPSI"]->setReal(dpsi);
    }

    // Fill in the "unknown" attributes.
    ElementImage *image = base.ElementBase::getImage();
    AttributeSet::const_iterator cur = image->begin();
    AttributeSet::const_iterator end = image->end();
    for(; cur != end; ++cur) {
        attributeRegistry[cur->first]->setReal(cur->second);
    }
}


AttCell *OpalElement::findRegisteredAttribute(const std::string &name) {
    AttCell *cell = &*attributeRegistry[name];

    if(cell == 0) {
        std::string::size_type i = 0;

        if(name[i] == 'K') {
            ++i;
            while(isdigit(name[i])) ++i;
            if(name[i] == 'S') ++i;

            if(name[i] == 'L'  &&  ++i == name.length()) {
                attributeRegistry[name] = cell = new AttReal();
            } else {
                throw OpalException("OpalElement::findRegisteredAttribute()",
                                    "There is no element which has an attribute "
                                    "called \"" + name + "\".");
            }
        }
    }

    return cell;
}

std::vector<double> OpalElement::getApert() const {
    return Attributes::getRealArray(itsAttr[APERT]);
}

double OpalElement::getLength() const {
    return Attributes::getReal(itsAttr[LENGTH]);
}


const std::string OpalElement::getTypeName() const {
    const Attribute *attr = findAttribute("TYPE");
    return attr ? Attributes::getString(*attr) : std::string();
}

/**
   Functions to get the wake field parametes
*/

const std::string OpalElement::getWakeF() const {
    const Attribute *attr = findAttribute("WAKEF");
    return attr ? Attributes::getString(*attr) : std::string();
}

const std::string OpalElement::getSurfacePhysics() const {
    const Attribute *attr = findAttribute("SURFACEPHYSICS");
    return attr ? Attributes::getString(*attr) : std::string();
}

void OpalElement::parse(Statement &stat) {
    while(stat.delimiter(',')) {
        std::string name = Expressions::parseString(stat, "Attribute name expected.");
        Attribute *attr = findAttribute(name);

        if(attr == 0) {
            // Append unknown attribute to the attribute list. Keep the original
            // size in itsSize, so as to recognize the added attributes.
            itsAttr.push_back(Attributes::makeReal(name, "user-defined"));
            attr = &itsAttr.back();
            registerRealAttribute(name);
        }

        if(stat.delimiter('[')) {
            int index = int(Round(Expressions::parseRealConst(stat)));
            Expressions::parseDelimiter(stat, ']');

            if(stat.delimiter('=')) {
                attr->parseComponent(stat, true, index);
            } else if(stat.delimiter(":=")) {
                attr->parseComponent(stat, false, index);
            } else {
                throw ParseError("OpalElement::parse()",
                                 "Delimiter \"=\" or \":=\" expected.");
            }
        } else {
            if(stat.delimiter('=')) {
                attr->parse(stat, true);
            } else if(stat.delimiter(":=")) {
                attr->parse(stat, false);
            } else {
                attr->setDefault();
            }
        }
    }
}


void OpalElement::print(std::ostream &os) const {
    std::string head = getOpalName();

    Object *parent = getParent();
    if(parent != 0  &&  ! parent->getOpalName().empty()) {
        if(! getOpalName().empty()) head += ':';
        head += parent->getOpalName();
    }

    os << head;
    os << ';'; // << "JMJdebug OPALElement.cc" ;
    os << std::endl;
}


void OpalElement::setRegisteredAttribute
(const std::string &name, double value) {
    attributeRegistry[name]->setReal(value);
}


void OpalElement::setRegisteredAttribute
(const std::string &name, const std::string &value) {
    attributeRegistry[name]->setString(value);
}


void OpalElement::printMultipoleStrength
(std::ostream &os, int order, int &len,
 const std::string &sName, const std::string &tName,
 const Attribute &length, const Attribute &sNorm, const Attribute &sSkew) {
    // Find out which type of output is required.
    int flag = 0;
    if(sNorm) {
        if(sNorm.getBase().isExpression()) {
            flag += 2;
        } else if(Attributes::getReal(sNorm) != 0.0) {
            flag += 1;
        }
    }

    if(sSkew) {
        if(sSkew.getBase().isExpression()) {
            flag += 6;
        } else if(Attributes::getReal(sSkew) != 0.0) {
            flag += 3;
        }
    }
    //  cout << "JMJdebug, OpalElement.cc: flag=" << flag << endl ;
    // Now do the output.
    int div = 2 * (order + 1);

    switch(flag) {

        case 0:
            // No component at all.
            break;

        case 1:
        case 2:
            // Pure normal component.
        {
            std::string normImage = sNorm.getImage();
            if(length) {
                normImage = "(" + normImage + ")*(" + length.getImage() + ")";
            }
            printAttribute(os, sName, normImage, len);
        }
        break;

        case 3:
        case 6:
            // Pure skew component.
        {
            std::string skewImage = sSkew.getImage();
            if(length) {
                skewImage = "(" + skewImage + ")*(" + length.getImage() + ")";
            }
            printAttribute(os, sName, skewImage, len);
            double tilt = Physics::pi / double(div);
            printAttribute(os, tName, tilt, len);
        }
        break;

        case 4:
            // Both components are non-zero constants.
        {
            double sn = Attributes::getReal(sNorm);
            double ss = Attributes::getReal(sSkew);
            double strength = sqrt(sn * sn + ss * ss);
            if(strength) {
#if defined(__GNUC__) && __GNUC__ < 3
                char buffer[80];
                std::ostrstream ts(buffer, 80);
#else
                std::ostringstream ts;
#endif
                ts << strength;
#if defined(__GNUC__) && __GNUC__ < 3
                std::string image(buffer);
#else
                std::string image = ts.str();
#endif
                if(length) {
                    image = "(" + image + ")*(" + length.getImage() + ")";
                }
                printAttribute(os, sName, image, len);
                double tilt = - atan2(ss, sn) / double(div);
                if(tilt) printAttribute(os, tName, tilt, len);
            }
        }
        break;

        case 5:
        case 7:
        case 8:
            // One or both components is/are expressions.
        {
            std::string normImage = sNorm.getImage();
            std::string skewImage = sSkew.getImage();
            std::string image =
                "SQRT((" + normImage + ")^2+(" + skewImage + ")^2)";
            printAttribute(os, sName, image, len);
            if(length) {
                image = "(" + image + ")*(" + length.getImage() + ")";
            }
            std::string divisor;
            if(div < 9) {
                divisor = "0";
                divisor[0] += div;
            } else {
                divisor = "00";
                divisor[0] += div / 10;
                divisor[1] += div % 10;
            }
            image = "-ATAN2(" + skewImage + ',' + normImage + ")/" + divisor;
            printAttribute(os, tName, image, len);
            break;
        }
    }
}


void OpalElement::updateUnknown(ElementBase *base) {
    for(std::vector<Attribute>::size_type i = itsSize;
        i < itsAttr.size(); ++i) {
        Attribute &attr = itsAttr[i];
        base->setAttribute(attr.getName(), Attributes::getReal(attr));

    }
}


void OpalElement::printAttribute
(std::ostream &os, const std::string &name, const std::string &image, int &len) {
    len += name.length() + image.length() + 2;
    if(len > 74) {
        os << ",&\n  ";
        len = name.length() + image.length() + 3;
    } else {
        os << ',';
    }
    os << name << '=' << image;
}

void OpalElement::printAttribute
(std::ostream &os, const std::string &name, double value, int &len) {
#if defined(__GNUC__) && __GNUC__ < 3
    char buffer[80];
    std::ostrstream ss(buffer, sizeof(buffer));
#else
    std::ostringstream ss;
#endif
    ss << value << std::ends;
#if defined(__GNUC__) && __GNUC__ < 3
    printAttribute(os, name, std::string(buffer), len);
#else
    printAttribute(os, name, ss.str(), len);
#endif
}


OpalElement::OpalElement(int size, const char *name, const char *help):
    Element(size, name, help), itsSize(size) {
    itsAttr[TYPE]   = Attributes::makeString
                      ("TYPE", "The element design type (the project name)");
    itsAttr[LENGTH] = Attributes::makeReal
                      ("L", "The element length in m");
    itsAttr[APERT]  = Attributes::makeRealArray
                      ("APERTURE", "The element aperture");
    itsAttr[WAKEF]   = Attributes::makeString
                       ("WAKEF", "Defines the wake function");
    itsAttr[SURFACEPHYSICS]   = Attributes::makeString
                                ("SURFACEPHYSICS", "Defines the surface physics handler");
    static bool first = true;
    if(first) {
        registerStringAttribute("NAME");
        registerStringAttribute("TYPE");
        registerStringAttribute("CLASS");
        registerStringAttribute("KEYWORD");
        registerRealAttribute("L");
        registerRealAttribute("DELTAX");
        registerRealAttribute("DELTAY");
        registerRealAttribute("DELTAS");
        registerRealAttribute("DTHETA");
        registerRealAttribute("DPHI");
        registerRealAttribute("DPSI");
        registerStringAttribute("WAKEF");
        registerStringAttribute("SURFACEPHYSICS");
        first = false;
    }

}


OpalElement::OpalElement(const std::string &name, OpalElement *parent):
    Element(name, parent), itsSize(parent->itsSize)
{}


AttCell *OpalElement::registerRealAttribute(const std::string &name) {
    OwnPtr<AttCell> &cell = attributeRegistry[name];
    if(! cell.isValid()) {
        cell = new AttReal();
    }
    return &*cell;
}


AttCell *OpalElement::registerStringAttribute(const std::string &name) {
    OwnPtr<AttCell> &cell = attributeRegistry[name];
    if(! cell.isValid()) {
        cell = new AttString();
    }
    return &*cell;
}
