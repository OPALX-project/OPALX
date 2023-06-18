//
// Class OpalBeamline
//   :FIXME: add class description
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Elements/OpalBeamline.h"

#include "AbstractObjects/OpalData.h"
#include "Physics/Units.h"
#include "Structure/MeshGenerator.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <fstream>

OpalBeamline::OpalBeamline():
    elements_m(),
    prepared_m(false)
{
}

OpalBeamline::OpalBeamline(const Vector_t& origin,
                           const Quaternion& rotation):
    elements_m(),
    prepared_m(false),
    coordTransformationTo_m(origin, rotation)
{
}

OpalBeamline::~OpalBeamline() {
    elements_m.clear();
}

std::set<std::shared_ptr<Component>> OpalBeamline::getElements(const Vector_t &x) {
    std::set<std::shared_ptr<Component> > elementSet;
    FieldList::iterator it = elements_m.begin();
    const FieldList::iterator end = elements_m.end();
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        Vector_t r = element->getCSTrafoGlobal2Local().transformTo(x);

        if (element->isInside(r)) {
            elementSet.insert(element);
        }
    }

    return elementSet;
}

unsigned long OpalBeamline::getFieldAt(const unsigned int &/*index*/, const Vector_t &/*pos*/, const long &/*sindex*/, const double &/*t*/, Vector_t &/*E*/, Vector_t &/*B*/) {

    unsigned long rtv = 0x00;

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(const Vector_t &position,
                                       const Vector_t &momentum,
                                       const double &t,
                                       Vector_t &Ef,
                                       Vector_t &Bf) {
    unsigned long rtv = 0x00;

    std::set<std::shared_ptr<Component>> elements = getElements(position);

    std::set<std::shared_ptr<Component>>::const_iterator it = elements.begin();
    const std::set<std::shared_ptr<Component>>::const_iterator end = elements.end();

    for (; it != end; ++ it) {
        ElementType type = (*it)->getType();
        if (type == ElementType::MARKER) continue;

        Vector_t localR = transformToLocalCS(*it, position);
        Vector_t localP = rotateToLocalCS(*it, momentum);
        Vector_t localE(0.0), localB(0.0);

        (*it)->applyToReferenceParticle(localR, localP, t, localE, localB);

        Ef += rotateFromLocalCS(*it, localE);
        Bf += rotateFromLocalCS(*it, localB);
    }

    //         if(section.hasWake()) {
    //             rtv |= BEAMLINE_WAKE;
    //         }
    //         if(section.hasParticleMatterInteraction()) {
    //             rtv |= BEAMLINE_PARTICLEMATTERINTERACTION;
    //         }

    return rtv;
}

void OpalBeamline::switchElements(const double &min, const double &max, const double &kineticEnergy, const bool &nomonitors) {

    FieldList::iterator fprev;
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit) {
        // don't set online monitors if the centroid of the bunch is allready inside monitor
        // or if explicitly not desired (eg during auto phasing)
        if(!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
            (*flit).setOn(kineticEnergy);
        }

        fprev = flit;
    }
}

void OpalBeamline::switchElementsOff() {
    for(FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++ flit)
        (*flit).setOff();
}

void OpalBeamline::prepareSections() {
    if (elements_m.empty()) {
        prepared_m = true;
        return;
    }
    elements_m.sort(ClassicField::SortAsc);
    prepared_m = true;
}

void OpalBeamline::print(Inform &/*msg*/) const {
}

void OpalBeamline::swap(OpalBeamline & rhs) {
    std::swap(elements_m, rhs.elements_m);
    std::swap(prepared_m, rhs.prepared_m);
    std::swap(coordTransformationTo_m, rhs.coordTransformationTo_m);
}

void OpalBeamline::merge(OpalBeamline &rhs) {
    elements_m.insert(elements_m.end(),
                      rhs.elements_m.begin(),
                      rhs.elements_m.end());
    prepared_m = false;
}


FieldList OpalBeamline::getElementByType(ElementType type) {
    if (type == ElementType::ANY) {
        return elements_m;
    }

    FieldList elements_of_requested_type;
    for(FieldList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++ fit) {
        if((*fit).getElement()->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}

void OpalBeamline::positionElementRelative(std::shared_ptr<ElementBase> element) {
    if (!element->isPositioned()) {
        return;
    }

    element->releasePosition();
    CoordinateSystemTrafo toElement = element->getCSTrafoGlobal2Local();
    toElement *= coordTransformationTo_m;

    element->setCSTrafoGlobal2Local(toElement);
    element->fixPosition();
}

void OpalBeamline::compute3DLattice() {

}

void OpalBeamline::save3DLattice() {

}

namespace {
    std::string parseInput() {

        std::ifstream in(OpalData::getInstance()->getInputFn());
        std::string source("");
        std::string str;
        char testBit;
        const std::string commentFormat("");
        const boost::regex empty("^[ \t]*$");
        const boost::regex lineEnd(";");
        const std::string lineEndFormat(";\n");
        const boost::regex cppCommentExpr("//.*");
        const boost::regex cCommentExpr("/\\*.*?\\*/"); // "/\\*(?>[^*/]+|\\*[^/]|/[^*])*(?>(?R)(?>[^*/]+|\\*[^/]|/[^*])*)*\\*/"
        bool priorEmpty = true;

        in.get(testBit);
        while (!in.eof()) {
            in.putback(testBit);

            std::getline(in, str);
            str = boost::regex_replace(str, cppCommentExpr, commentFormat);
            str = boost::regex_replace(str, empty, commentFormat);
            if (!str.empty()) {
                source += str;// + '\n';
                priorEmpty = false;
            } else if (!priorEmpty) {
                source += "##EMPTY_LINE##";
                priorEmpty = true;
            }

            in.get(testBit);
        }

        source = boost::regex_replace(source, cCommentExpr, commentFormat);
        source = boost::regex_replace(source, lineEnd, lineEndFormat, boost::match_default | boost::format_all);

        // Since the positions of the elements are absolute in the laboratory coordinate system we have to make
        // sure that the line command doesn't provide an origin and orientation. Everything after the sequence of
        // elements can be deleted and only "LINE = (...);", the first sub-expression (denoted by '\1'), should be kept.
        const boost::regex lineCommand("(LINE[ \t]*=[ \t]*\\([^\\)]*\\))[ \t]*,[^;]*;", boost::regex::icase);
        const std::string firstSubExpression("\\1;");
        source = boost::regex_replace(source, lineCommand, firstSubExpression);

        return source;
    }

    unsigned int getMinimalSignificantDigits(double num, const unsigned int maxDigits) {
        char buf[32];
        snprintf(buf, 32, "%.*f", maxDigits + 1, num);
        std::string numStr(buf);
        unsigned int length = numStr.length();

        unsigned int numDigits = maxDigits;
        unsigned int i = 2;
        while (i < maxDigits + 1 && numStr[length - i] == '0') {
            --numDigits;
            ++ i;
        }

        return numDigits;
    }

    std::string round2string(double num, const unsigned int maxDigits) {
        char buf[64];

        snprintf(buf, 64, "%.*f", getMinimalSignificantDigits(num, maxDigits), num);

        return std::string(buf);
    }
}

void OpalBeamline::save3DInput() {
    if (Ippl::myNode() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    FieldList::iterator it = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::string input = parseInput();
    std::string fname = Util::combineFilePath({
        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
        OpalData::getInstance()->getInputBasename() + "_3D.opal"
    });
    std::ofstream pos(fname);

    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        std::string elementName = element->getName();
        const boost::regex replacePSI("(" + elementName + "\\s*:[^\\n]*)PSI\\s*=[^,;]*,?", boost::regex::icase);
        input = boost::regex_replace(input, replacePSI, "\\1\\2");

        const boost::regex replaceELEMEDGE("(" + elementName + "\\s*:[^\\n]*)ELEMEDGE\\s*=[^,;]*(.)", boost::regex::icase);

        CoordinateSystemTrafo cst = element->getCSTrafoGlobal2Local();
        Vector_t origin = cst.getOrigin();
        Vector_t orient = Util::getTaitBryantAngles(cst.getRotation().conjugate(), elementName);
        for (unsigned int d = 0; d < 3; ++ d)
            orient(d) *= Units::rad2deg;

        std::string x = (std::abs(origin(0)) > 1e-10? "X = " + round2string(origin(0), 10) + ", ": "");
        std::string y = (std::abs(origin(1)) > 1e-10? "Y = " + round2string(origin(1), 10) + ", ": "");
        std::string z = (std::abs(origin(2)) > 1e-10? "Z = " + round2string(origin(2), 10) + ", ": "");

        std::string theta = (orient(0) > 1e-10? "THETA = " + round2string(orient(0), 6) + " * PI / 180, ": "");
        std::string phi = (orient(1) > 1e-10? "PHI = " + round2string(orient(1), 6) + " * PI / 180, ": "");
        std::string psi = (orient(2) > 1e-10? "PSI = " + round2string(orient(2), 6) + " * PI / 180, ": "");
        std::string coordTrafo = x + y + z + theta + phi + psi;
        if (coordTrafo.length() > 2) {
            coordTrafo = coordTrafo.substr(0, coordTrafo.length() - 2); // remove last ', '
        }

        std::string position = ("\\1" + coordTrafo + "\\2");

        input = boost::regex_replace(input, replaceELEMEDGE, position);
    }

    const boost::regex empty("##EMPTY_LINE##");
    const std::string emptyFormat("\n");
    input = boost::regex_replace(input, empty, emptyFormat);

    pos << input << std::endl;
}

void OpalBeamline::activateElements() {
    auto it = elements_m.begin();
    const auto end = elements_m.end();

    double designEnergy = 0.0;
    for (; it != end; ++ it) {
        std::shared_ptr<Component> element = (*it).getElement();
        (*it).setOn(designEnergy);
        // element->goOnline(designEnergy);
    }
}