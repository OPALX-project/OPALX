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

#include "Elements/PlacementResolver.h"

#include "AbstractObjects/OpalData.h"
#include "Physics/Units.h"
#include "Structure/MeshGenerator.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <filesystem>
#include <fstream>
#include <regex>

OpalBeamline::OpalBeamline()
    : elements_m(), prepared_m(false), referencePathPlacementCompiled_m(false) {}

OpalBeamline::OpalBeamline(const Vector_t<double, 3>& origin, const Quaternion& rotation)
    : elements_m(),
      prepared_m(false),
      referencePathPlacementCompiled_m(false),
      coordTransformationTo_m(origin, rotation) {}

OpalBeamline::~OpalBeamline() { elements_m.clear(); }

std::set<std::shared_ptr<ElementBase>> OpalBeamline::getElements(const Vector_t<double, 3>& x) {
    std::set<std::shared_ptr<ElementBase>> elementSet;
    FieldList::iterator it        = elements_m.begin();
    const FieldList::iterator end = elements_m.end();
    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it).getElement();
        Vector_t<double, 3> r                = getCSTrafoLab2Local(element).transformTo(x);

        if (element->isInside(r)) {
            elementSet.insert(element);
        }
    }

    return elementSet;
}

std::set<std::shared_ptr<ElementBase>> OpalBeamline::getElements() {
    std::set<std::shared_ptr<ElementBase>> elementSet;
    for (auto& item : elements_m) {
        elementSet.insert(item.getElement());
    }
    return elementSet;
}

unsigned long OpalBeamline::getFieldAt(
        const unsigned int& /*index*/, const Vector_t<double, 3>& /*pos*/, const long& /*sindex*/,
        const double& /*t*/, Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    unsigned long rtv = 0x00;

    return rtv;
}

unsigned long OpalBeamline::getFieldAt(
        const Vector_t<double, 3>& position, const Vector_t<double, 3>& momentum, const double& t,
        Vector_t<double, 3>& Ef, Vector_t<double, 3>& Bf) {
    unsigned long rtv = 0x00;

    std::set<std::shared_ptr<ElementBase>> elements = getElements(position);

    std::set<std::shared_ptr<ElementBase>>::const_iterator it        = elements.begin();
    const std::set<std::shared_ptr<ElementBase>>::const_iterator end = elements.end();

    for (; it != end; ++it) {
        ElementType type = (*it)->getType();
        if (type == ElementType::MARKER) continue;

        Vector_t<double, 3> localR = transformToLocalCS(*it, position);
        Vector_t<double, 3> localP = rotateToLocalCS(*it, momentum);
        Vector_t<double, 3> localE(0.0), localB(0.0);

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

void OpalBeamline::switchElements(
        const double& min, const double& max, const double& kineticEnergy,
        const bool& /*nomonitors*/) {
    FieldList::iterator fprev;
    for (FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++flit) {
        // don't set online monitors if the centroid of the bunch is allready inside monitor
        // or if explicitly not desired (eg during auto phasing)
        if (!(*flit).isOn() && max > (*flit).getStart() && min < (*flit).getEnd()) {
            (*flit).setOn(kineticEnergy);
        }

        fprev = flit;
    }
}

void OpalBeamline::switchElementsOff() {
    for (FieldList::iterator flit = elements_m.begin(); flit != elements_m.end(); ++flit)
        (*flit).setOff();
}

void OpalBeamline::prepareSections() {
    if (elements_m.empty()) {
        prepared_m = true;
        return;
    }
    elements_m.sort(BeamlineFieldElement::SortAsc);
    placeElementsAlongReferencePath();
    prepared_m = true;
}

void OpalBeamline::print(Inform& /*msg*/) const {}

void OpalBeamline::swap(OpalBeamline& rhs) {
    std::swap(elements_m, rhs.elements_m);
    std::swap(prepared_m, rhs.prepared_m);
    std::swap(referencePathPlacementCompiled_m, rhs.referencePathPlacementCompiled_m);
    std::swap(coordTransformationTo_m, rhs.coordTransformationTo_m);
}

void OpalBeamline::merge(OpalBeamline& rhs) {
    elements_m.insert(elements_m.end(), rhs.elements_m.begin(), rhs.elements_m.end());
    prepared_m                       = false;
    referencePathPlacementCompiled_m = false;
}

FieldList OpalBeamline::getElementByType(ElementType type) {
    if (type == ElementType::ANY) {
        return elements_m;
    }

    FieldList elements_of_requested_type;
    for (FieldList::iterator fit = elements_m.begin(); fit != elements_m.end(); ++fit) {
        if ((*fit).getElement()->getType() == type) {
            elements_of_requested_type.push_back((*fit));
        }
    }
    return elements_of_requested_type;
}

void OpalBeamline::placeElementsAlongReferencePath() {
    if (referencePathPlacementCompiled_m) {
        return;
    }
    // Single PLACE stage: resolve every element's global-to-local transform (6D pose +
    // ELEMEDGE) in one place.
    PlacementResolver::resolve(elements_m, coordTransformationTo_m);
    referencePathPlacementCompiled_m = true;
}

void OpalBeamline::compute3DLattice() { placeElementsAlongReferencePath(); }

void OpalBeamline::save3DLattice() {
    if (ippl::Comm->rank() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    elements_m.sort([](const BeamlineFieldElement& a, const BeamlineFieldElement& b) {
        return a.order_m < b.order_m;
    });

    FieldList::iterator it  = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::ofstream pos;
    std::string fileName = Util::combineFilePath(
            {OpalData::getInstance()->getAuxiliaryOutputDirectory(),
             OpalData::getInstance()->getInputBasename() + "_ElementPositions.txt"});
    if (OpalData::getInstance()->getOpenMode() == OpalData::OpenMode::APPEND
        && std::filesystem::exists(fileName)) {
        pos.open(fileName, std::ios_base::app);
    } else {
        pos.open(fileName);
    }

    MeshGenerator mesh;
    for (auto scan = it; scan != end; ++scan) {
        const std::shared_ptr<ElementBase> scanElement = (*scan).getElement();
        if (scanElement->getType() == ElementType::DRIFT) {
            continue;
        }

        double minor = 0.0;
        double major = 0.0;
        if (MeshGenerator::getTransverseSupport(*scanElement, minor, major)) {
            mesh.setDriftReference(0.5 * minor, 0.5 * major);
            break;
        }
    }

    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it).getElement();
        CoordinateSystemTrafo nominalBody    = getCSTrafoLab2Local(element);
        CoordinateSystemTrafo toBegin        = getNominalEntryTransform(element);
        CoordinateSystemTrafo toEnd          = getNominalExitTransform(element);
        Vector_t<double, 3> entry3D          = toBegin.getOrigin();
        Vector_t<double, 3> exit3D           = toEnd.getOrigin();

        mesh.add(*(element.get()));

        if (element->getType() == ElementType::SBEND || element->getType() == ElementType::RBEND) {
            std::vector<Vector_t<double, 3>> designPath = element->getDesignPath();
            unsigned int size                           = designPath.size();

            unsigned int minNumSteps = std::max(
                    20u, static_cast<unsigned int>(
                                 std::ceil(std::abs(element->getBendAngle() * Units::rad2deg))));

            unsigned int frequency =
                    std::max(1u, static_cast<unsigned int>(std::floor((double)size / minNumSteps)));

            pos << std::setw(30) << std::left
                << std::string("\"ENTRY EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2) << std::setw(18)
                << std::setprecision(10) << entry3D(0) << std::setw(18) << std::setprecision(10)
                << entry3D(1) << "\n";

            Vector_t<double, 3> position = nominalBody.transformFrom(designPath.front());
            pos << std::setw(30) << std::left
                << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << position(2) << std::setw(18)
                << std::setprecision(10) << position(0) << std::setw(18) << std::setprecision(10)
                << position(1) << std::endl;

            for (unsigned int i = frequency; i + 1 < size; i += frequency) {
                position = nominalBody.transformFrom(designPath[i]);
                pos << std::setw(30) << std::left
                    << std::string("\"MID: ") + element->getName() + std::string("\"")
                    << std::setw(18) << std::setprecision(10) << position(2) << std::setw(18)
                    << std::setprecision(10) << position(0) << std::setw(18)
                    << std::setprecision(10) << position(1) << std::endl;
            }

            position = nominalBody.transformFrom(designPath.back());
            pos << std::setw(30) << std::left
                << std::string("\"END: ") + element->getName() + std::string("\"") << std::setw(18)
                << std::setprecision(10) << position(2) << std::setw(18) << std::setprecision(10)
                << position(0) << std::setw(18) << std::setprecision(10) << position(1)
                << std::endl;

            pos << std::setw(30) << std::left
                << std::string("\"EXIT EDGE: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << exit3D(2) << std::setw(18)
                << std::setprecision(10) << exit3D(0) << std::setw(18) << std::setprecision(10)
                << exit3D(1) << std::endl;
        } else {
            pos << std::setw(30) << std::left
                << std::string("\"BEGIN: ") + element->getName() + std::string("\"")
                << std::setw(18) << std::setprecision(10) << entry3D(2) << std::setw(18)
                << std::setprecision(10) << entry3D(0) << std::setw(18) << std::setprecision(10)
                << entry3D(1) << "\n";

            pos << std::setw(30) << std::left
                << std::string("\"END: ") + element->getName() + std::string("\"") << std::setw(18)
                << std::setprecision(10) << exit3D(2) << std::setw(18) << std::setprecision(10)
                << exit3D(0) << std::setw(18) << std::setprecision(10) << exit3D(1) << std::endl;
        }
    }
    elements_m.sort(BeamlineFieldElement::SortAsc);
    mesh.write(OpalData::getInstance()->getInputBasename());
}

namespace {
    std::string parseInput() {
        std::ifstream in(OpalData::getInstance()->getInputFn());
        std::string source("");
        std::string str;
        char testBit;
        const std::string commentFormat("");
        const std::regex empty("^[ \t]*$");
        const std::regex lineEnd(";");
        const std::string lineEndFormat(";\n");
        const std::regex cppCommentExpr("//.*");
        const std::regex cCommentExpr(
                "/\\*.*?\\*/");  // "/\\*(?>[^*/]+|\\*[^/]|/[^*])*(?>(?R)(?>[^*/]+|\\*[^/]|/[^*])*)*\\*/"
        bool priorEmpty = true;

        in.get(testBit);
        while (!in.eof()) {
            in.putback(testBit);

            std::getline(in, str);
            str = std::regex_replace(
                    str, cppCommentExpr, commentFormat, std::regex_constants::format_default);
            str = std::regex_replace(
                    str, empty, commentFormat, std::regex_constants::format_default);
            if (!str.empty()) {
                source += str;  // + '\n';
                priorEmpty = false;
            } else if (!priorEmpty) {
                source += "##EMPTY_LINE##";
                priorEmpty = true;
            }

            in.get(testBit);
        }

        source = std::regex_replace(source, cCommentExpr, commentFormat);
        source = std::regex_replace(
                source, lineEnd, lineEndFormat,
                std::regex_constants::match_default | std::regex_constants::format_default);

        // Since the positions of the elements are absolute in the laboratory coordinate system we
        // have to make sure that the line command doesn't provide an origin and orientation.
        // Everything after the sequence of elements can be deleted and only "LINE = (...);", the
        // first sub-expression (denoted by '\1'), should be kept.
        const std::regex lineCommand(
                "(LINE[ \t]*=[ \t]*\\([^\\)]*\\))[ \t]*,[^;]*;", std::regex::icase);
        const std::string firstSubExpression("\\1;");
        source = std::regex_replace(source, lineCommand, firstSubExpression);

        return source;
    }

    unsigned int getMinimalSignificantDigits(double num, const unsigned int maxDigits) {
        char buf[32];
        snprintf(buf, 32, "%.*f", maxDigits + 1, num);
        std::string numStr(buf);
        unsigned int length = numStr.length();

        unsigned int numDigits = maxDigits;
        unsigned int i         = 2;
        while (i < maxDigits + 1 && numStr[length - i] == '0') {
            --numDigits;
            ++i;
        }

        return numDigits;
    }

    std::string round2string(double num, const unsigned int maxDigits) {
        char buf[64];

        snprintf(buf, 64, "%.*f", getMinimalSignificantDigits(num, maxDigits), num);

        return std::string(buf);
    }
}  // namespace

void OpalBeamline::save3DInput() {
    if (ippl::Comm->rank() != 0 || OpalData::getInstance()->isOptimizerRun()) return;

    FieldList::iterator it  = elements_m.begin();
    FieldList::iterator end = elements_m.end();

    std::string input = parseInput();
    std::string fname = Util::combineFilePath(
            {OpalData::getInstance()->getAuxiliaryOutputDirectory(),
             OpalData::getInstance()->getInputBasename() + "_3D.opal"});
    std::ofstream pos(fname);

    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it).getElement();
        std::string elementName              = element->getName();
        const std::regex replacePSI(
                "(" + elementName + "\\s*:[^\\n]*)PSI\\s*=[^,;]*,?", std::regex::icase);
        input = std::regex_replace(input, replacePSI, "\\1\\2");

        const std::regex replaceELEMEDGE(
                "(" + elementName + "\\s*:[^\\n]*)ELEMEDGE\\s*=[^,;]*(.)", std::regex::icase);

        CoordinateSystemTrafo cst  = getCSTrafoLab2Local(element);
        Vector_t<double, 3> origin = cst.getOrigin();
        Vector_t<double, 3> orient =
                Util::getTaitBryantAngles(cst.getRotation().conjugate(), elementName);
        for (unsigned int d = 0; d < 3; ++d)
            orient(d) *= Units::rad2deg;

        std::string x =
                (std::abs(origin(0)) > 1e-10 ? "X = " + round2string(origin(0), 10) + ", " : "");
        std::string y =
                (std::abs(origin(1)) > 1e-10 ? "Y = " + round2string(origin(1), 10) + ", " : "");
        std::string z =
                (std::abs(origin(2)) > 1e-10 ? "Z = " + round2string(origin(2), 10) + ", " : "");

        std::string theta =
                (orient(0) > 1e-10 ? "THETA = " + round2string(orient(0), 6) + " * PI / 180, "
                                   : "");
        std::string phi =
                (orient(1) > 1e-10 ? "PHI = " + round2string(orient(1), 6) + " * PI / 180, " : "");
        std::string psi =
                (orient(2) > 1e-10 ? "PSI = " + round2string(orient(2), 6) + " * PI / 180, " : "");
        std::string coordTrafo = x + y + z + theta + phi + psi;
        if (coordTrafo.length() > 2) {
            coordTrafo = coordTrafo.substr(0, coordTrafo.length() - 2);  // remove last ', '
        }

        std::string position = ("\\1" + coordTrafo + "\\2");

        input = std::regex_replace(input, replaceELEMEDGE, position);
    }

    const std::regex empty("##EMPTY_LINE##");
    const std::string emptyFormat("\n");
    input = std::regex_replace(input, empty, emptyFormat);

    pos << input << std::endl;
}

void OpalBeamline::activateElements() {
    auto it             = elements_m.begin();
    const auto end      = elements_m.end();
    double designEnergy = 0.0;
    for (; it != end; ++it) {
        std::shared_ptr<ElementBase> element = (*it).getElement();
        (*it).setOn(designEnergy);
        element->goOnline(designEnergy);
    }
}
