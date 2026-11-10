// Copyright (c) 2023, Chris Rogers
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

#include "Elements/OpalOutputPlane.h"
#include "AbsBeamline/OutputPlane.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"

const std::string docstring =
    std::string("The \"OUTPUTPLANE\" element writes out position at which")+
    std::string("trajectories cross a given plane.");

OpalOutputPlane::OpalOutputPlane():
    OpalElement(SIZE, "OUTPUTPLANE", docstring.c_str()) {

    itsAttr[CENTRE] = Attributes::makeRealArray
                      ("CENTRE", "3-vector position of the plane centre [m]");
    itsAttr[NORMAL] = Attributes::makeRealArray
                      ("NORMAL", "3-vector normal to the plane");
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Maximum distance from centre of plane for crossings [m].");
    registerOwnership();
    setElement(new OutputPlane("OUTPUTPLANE"));
}


OpalOutputPlane::OpalOutputPlane(const std::string &name, OpalOutputPlane *parent):
    OpalElement(name, parent) {
    OutputPlane* plane = dynamic_cast<OutputPlane*>(parent->getElement());
    setElement(new OutputPlane(*plane));
}


OpalOutputPlane::~OpalOutputPlane() {
}


OpalOutputPlane *OpalOutputPlane::clone(const std::string &name) {
    return new OpalOutputPlane(name, this);
}

void OpalOutputPlane::update() {
    OpalElement::update();
    std::vector<double> centreVec = Attributes::getRealArray(itsAttr[CENTRE]);
    //if (centreVec.size() != 3) {
    //    throw OpalException("OpalOutputPlane::update()",
    //        "OutputPlane centre should be a 3-vector");
    //}
    Vector_t<double, 3> centre;
    if (centreVec.size() == 3) {
        centre = Vector_t<double, 3>({centreVec[0], centreVec[1], centreVec[2]});
    }
    std::vector<double> normalVec = Attributes::getRealArray(itsAttr[NORMAL]);
    //if (normalVec.size() != 3) {
    //    throw OpalException("OpalOutputPlane::update()",
    //        "OutputPlane normal should be a 3-vector");
    //}
    Vector_t<double, 3> normal;
    if (normalVec.size() == 3) {
        normal = Vector_t<double, 3>({normalVec[0], normalVec[1], normalVec[2]});
    }
    const double radius = Attributes::getReal(itsAttr[RADIUS]);
    //if (radius < 0) {
    //    throw OpalException("OpalOutputPlane::update()",
    //        "OutputPlane had negative radius");
    //}
    auto *output = dynamic_cast<OutputPlane *>(getElement());
    output->setCentre(centre);
    output->setNormal(normal);
    output->setRadialExtent(radius);
    output->setOutputFN(Attributes::getString(itsAttr[OUTFN]));
    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(output);
}
