//
// Class OpalBeamBeam
//   The class of OPAL drift spaces.
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
#include "Elements/OpalBeamBeam.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/BeamBeamRep.h"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <exception>
#include <sstream>
#include <string>

namespace {

    std::string normalizeWitnessContainerSpec(std::string value) {
        value.erase(
                std::remove_if(
                        value.begin(), value.end(),
                        [](unsigned char ch) {
                            return std::isspace(ch) != 0;
                        }),
                value.end());
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
            return static_cast<char>(std::toupper(ch));
        });
        return value;
    }

    double parseWitnessContainerMask(const std::string& rawValue) {
        const std::string normalized = normalizeWitnessContainerSpec(rawValue);
        if (normalized == "NONE") {
            return 0.0;
        }
        if (normalized.empty()) {
            throw OpalException(
                    "OpalBeamBeam::parseWitnessContainerMask",
                    "WITNESS_CONTAINERS must be \"NONE\" or a comma-separated list such as "
                    "\"1,2\".");
        }

        unsigned long long mask = 0ULL;
        std::stringstream stream(normalized);
        std::string token;
        while (std::getline(stream, token, ',')) {
            if (token.empty()) {
                throw OpalException(
                        "OpalBeamBeam::parseWitnessContainerMask",
                        "WITNESS_CONTAINERS contains an empty container index.");
            }
            for (const char ch : token) {
                if (std::isdigit(static_cast<unsigned char>(ch)) == 0) {
                    throw OpalException(
                            "OpalBeamBeam::parseWitnessContainerMask",
                            "WITNESS_CONTAINERS contains non-integer entry \"" + token + "\".");
                }
            }

            std::size_t index = 0;
            try {
                index = std::stoull(token);
            } catch (const std::exception&) {
                throw OpalException(
                        "OpalBeamBeam::parseWitnessContainerMask",
                        "WITNESS_CONTAINERS entry \"" + token + "\" is out of range.");
            }
            if (index == 0) {
                throw OpalException(
                        "OpalBeamBeam::parseWitnessContainerMask",
                        "container[0] is the BeamBeam source and cannot be a witness container.");
            }
            if (index >= 53) {
                throw OpalException(
                        "OpalBeamBeam::parseWitnessContainerMask",
                        "WITNESS_CONTAINERS supports indices below 53 when encoded as a scalar "
                        "mask.");
            }
            mask |= (1ULL << index);
        }

        return static_cast<double>(mask);
    }

}  // namespace

OpalBeamBeam::OpalBeamBeam()
    : OpalElement(
              SIZE, "BEAMBEAM",
              "The \"BEAMBEAM\" element defines a beam-beam interaction point for colliding "
              "beams.") {
    itsAttr[GEOMETRY] =
            Attributes::makeString("GEOMETRY", "BoundaryGeometry for beam-beam elements");
    itsAttr[COPY_TIME] = Attributes::makeReal(
            "COPY_TIME",
            "Start the mirrored-bunch copy model once simulation time reaches this value [s]. "
            "Use 0 to disable copied fields.",
            0.0);
    itsAttr[VISUALIZE] = Attributes::makeBool(
            "VISUALIZE", "If true, emit the ASCII beam-beam-window visualization during tracking.",
            false);
    itsAttr[WITNESS_CONTAINERS] = Attributes::makeString(
            "WITNESS_CONTAINERS",
            "Passive particle-container indices that gather the source BeamBeam field, or NONE.",
            "NONE");
    itsAttr[RETIRE_TIME] = Attributes::makeReal(
            "RETIRE_TIME",
            "Delete container[0] source particles once simulation time reaches this value [s]. "
            "Use 0 to keep the source particles active.",
            0.0);

    registerOwnership();

    setElement(new BeamBeamRep("DRIFT"));
}

OpalBeamBeam::OpalBeamBeam(const std::string& name, OpalBeamBeam* parent)
    : OpalElement(name, parent) {
    setElement(new BeamBeamRep(name));
}

OpalBeamBeam::~OpalBeamBeam() {}

OpalBeamBeam* OpalBeamBeam::clone(const std::string& name) { return new OpalBeamBeam(name, this); }

bool OpalBeamBeam::isBeamBeam() const { return true; }

void OpalBeamBeam::update() {
    OpalElement::update();

    BeamBeamRep* beamBeam = static_cast<BeamBeamRep*>(getElement());
    beamBeam->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    const double copyTime = Attributes::getReal(itsAttr[COPY_TIME]);
    if (copyTime < 0.0) {
        throw OpalException(
                "OpalBeamBeam::update", "COPY_TIME must be non-negative and is specified in s.");
    }
    beamBeam->setAttribute("COPY_TIME", copyTime);
    beamBeam->setAttribute("VISUALIZE", Attributes::getBool(itsAttr[VISUALIZE]) ? 1.0 : 0.0);
    beamBeam->setAttribute(
            "WITNESS_CONTAINERS_MASK",
            parseWitnessContainerMask(Attributes::getString(itsAttr[WITNESS_CONTAINERS])));
    const double retireTime = Attributes::getReal(itsAttr[RETIRE_TIME]);
    if (retireTime < 0.0) {
        throw OpalException(
                "OpalBeamBeam::update", "RETIRE_TIME must be non-negative and is specified in s.");
    }
    beamBeam->setAttribute("RETIRE_TIME", retireTime);
    beamBeam->setAttribute("APERTURE_SET", itsAttr[APERT] ? 1.0 : 0.0);

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(beamBeam);
}
