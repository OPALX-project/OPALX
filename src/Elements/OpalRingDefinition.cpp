/* 
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <limits>

#include "Elements/OpalRingDefinition.h"
#include "Elements/OpalRing.h"

#include "Attributes/Attributes.h"

OpalRingDefinition::OpalRingDefinition() :
    OpalElement(SIZE, "RINGDEFINITION",
              "The \"RINGDEFINITION\" element defines basic ring parameters.") {
    itsAttr[HARMONIC_NUMBER] = Attributes::makeReal("HARMONIC_NUMBER",
        "The assumed harmonic number of the ring (i.e. number of bunches in the ring on a given turn).");
    itsAttr[LATTICE_RINIT] = Attributes::makeReal("LATTICE_RINIT",
        "The initial radius of the first element to be placed in the ring [mm].");
    itsAttr[LATTICE_PHIINIT] = Attributes::makeReal("LATTICE_PHIINIT",
        "The initial angle around the ring of the first element to be placed. [deg]");
    itsAttr[LATTICE_THETAINIT] = Attributes::makeReal("LATTICE_THETAINIT",
        "The angle relative to the tangent of the ring for the first element to be placed [deg].");
    itsAttr[BEAM_PHIINIT] = Attributes::makeReal("BEAM_PHIINIT",
        "The initial angle around the ring of the beam [deg].");
    itsAttr[BEAM_PRINIT] = Attributes::makeReal("BEAM_PRINIT",
        "An initial pr momentum offset of the beam.");
    itsAttr[BEAM_RINIT] = Attributes::makeReal("BEAM_RINIT",
        "The initial radius of the beam [mm].");
    itsAttr[SYMMETRY] = Attributes::makeReal("SYMMETRY",
        "The rotational symmetry of the lattice.");
    // should be in RF cavity definition; this comes from cyclotron definition,
    // but not right
    itsAttr[RFFREQ] = Attributes::makeReal("RFFREQ",
        "The nominal RF frequency of the ring [MHz].");
    // I see also makeBool, but dont know how it works; no registerBoolAttribute
    itsAttr[IS_CLOSED] = Attributes::makeString("IS_CLOSED",
        "Set to 'false' to disable checking for closure of the ring");

    registerRealAttribute("LATTICE_RINIT");
    registerRealAttribute("LATTICE_PHIINIT");
    registerRealAttribute("LATTICE_THETAINIT");
    registerRealAttribute("BEAM_RINIT");
    registerRealAttribute("BEAM_PHIINIT");
    registerRealAttribute("BEAM_PRINIT");
    registerRealAttribute("HARMONIC_NUMBER");
    registerRealAttribute("SYMMETRY");
    registerRealAttribute("RFFREQ");
    registerStringAttribute("IS_CLOSED");

    setElement((new OpalRing("OPALRING"))->makeAlignWrapper());
}

OpalRingDefinition* OpalRingDefinition::clone(const string &name) {
    return new OpalRingDefinition(name, this);
}

void OpalRingDefinition::print(std::ostream& out) const {
    OpalElement::print(out);
}

OpalRingDefinition::OpalRingDefinition(const string &name, OpalRingDefinition *parent):
    OpalElement(name, parent) {
    setElement((new OpalRing(name))->makeAlignWrapper());
}

OpalRingDefinition::~OpalRingDefinition() {}

void OpalRingDefinition::fillRegisteredAttributes
                                     (const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}

void OpalRingDefinition::update() {
    OpalRing *ring = dynamic_cast<OpalRing*>(getElement()->removeWrappers());
    double degree = Physics::pi/180.;
    ring->setBeamPhiInit(Attributes::getReal(itsAttr[BEAM_PHIINIT]));
    ring->setBeamPRInit(Attributes::getReal(itsAttr[BEAM_PRINIT]));
    ring->setBeamRInit(Attributes::getReal(itsAttr[BEAM_RINIT]));
    ring->setLatticeRInit(Attributes::getReal(itsAttr[LATTICE_RINIT]));
    ring->setLatticePhiInit
                         (Attributes::getReal(itsAttr[LATTICE_PHIINIT])*degree);
    ring->setLatticeThetaInit
                       (Attributes::getReal(itsAttr[LATTICE_THETAINIT])*degree);
    ring->setSymmetry(Attributes::getReal(itsAttr[SYMMETRY]));
    ring->setHarmonicNumber(Attributes::getReal(itsAttr[HARMONIC_NUMBER]));
    ring->setRFFreq(Attributes::getReal(itsAttr[RFFREQ]));
    ring->setIsClosed(!(Attributes::getString(itsAttr[IS_CLOSED])=="FALSE"));

    setElement(ring->makeWrappers());
}

