/* 
 *  Copyright (c) 2017, Chris Rogers
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

#include "Attributes/Attributes.h"  // used?
#include "Utilities/OpalException.h"  // used?

#include "AbsBeamline/EndFieldModel/Tanh.h" // classic
#include "AbsBeamline/ScalingFFAGMagnet.h" // classic
#include "Elements/OpalScalingFFAGMagnet.h"

extern Inform *gmsg;

OpalScalingFFAGMagnet::OpalScalingFFAGMagnet() :
    OpalElement(SIZE, "SCALINGFFAGMAGNET",
             "The \"ScalingFFAGMagnet\" element defines a FFAG scaling magnet with zero or non-zero spiral angle.") {
    itsAttr[B0] = Attributes::makeReal
                              ("B0", "The nominal dipole field of the magnet [?].");
    itsAttr[R0] = Attributes::makeReal("R0", "Radial scale [mm].");
    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX",
      "The scaling magnet field index.");
    itsAttr[TAN_DELTA] = Attributes::makeReal("TAN_DELTA",
      "Tangent of the spiral angle; set to 0 to make a radial sector magnet.");
    itsAttr[MAX_Y_POWER] = Attributes::makeReal("MAX_Y_POWER",
      "The maximum power in y that will be considered in the field expansion.");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
                                          "The end length of the spiral FFAG [rad].");
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
                                       "The centre length of the spiral FFAG [rad].");
    itsAttr[RADIAL_NEG_EXTENT] = Attributes::makeReal("RADIAL_NEG_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0-RADIAL_NEG_EXTENT [mm].");
    itsAttr[RADIAL_POS_EXTENT] = Attributes::makeReal("RADIAL_POS_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0+RADIAL_POS_EXTENT [mm].");
    itsAttr[PHI_START] = Attributes::makeReal("PHI_START",
                                          "Offset the centre of the spiral FFAG magnet (default is centre_length/2 + end_length).");
    itsAttr[PHI_END] = Attributes::makeReal("PHI_END",
                                       "Offset to the end of the magnet, i.e. placement of the next element. Default is centre_length + 2*end_length.");
    itsAttr[AZIMUTHAL_EXTENT] = Attributes::makeReal("AZIMUTHAL_EXTENT",
                                       "The field will be assumed zero if particles are more than AZIMUTHAL_EXTENT from the magnet centre (psi=0). Default is CENTRE_LENGTH/2.+5.*END_LENGTH [rad].");
    registerRealAttribute("B0");
    registerRealAttribute("R0");
    registerRealAttribute("FIELD_INDEX");
    registerRealAttribute("TAN_DELTA");
    registerRealAttribute("MAX_Y_POWER");
    registerRealAttribute("END_LENGTH");
    registerRealAttribute("CENTRE_LENGTH");
    registerRealAttribute("RADIAL_NEG_EXTENT");
    registerRealAttribute("RADIAL_POS_EXTENT");
    registerRealAttribute("PHI_START");
    registerRealAttribute("PHI_END");
    registerRealAttribute("AZIMUTHAL_EXTENT");
    registerOwnership();

    ScalingFFAGMagnet* magnet = new ScalingFFAGMagnet("ScalingFFAGMagnet");
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet->makeAlignWrapper());
}


OpalScalingFFAGMagnet::OpalScalingFFAGMagnet(const std::string &name,
                                             OpalScalingFFAGMagnet *parent) :
    OpalElement(name, parent) {
    ScalingFFAGMagnet* magnet = new ScalingFFAGMagnet(name);
    magnet->setEndField(new endfieldmodel::Tanh(1., 1., 1));
    setElement(magnet->makeAlignWrapper());
}


OpalScalingFFAGMagnet::~OpalScalingFFAGMagnet() {
}


OpalScalingFFAGMagnet *OpalScalingFFAGMagnet::clone(const std::string &name) {
    return new OpalScalingFFAGMagnet(name, this);
}


void OpalScalingFFAGMagnet::
fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}


void OpalScalingFFAGMagnet::update() {
    ScalingFFAGMagnet *magnet = dynamic_cast<ScalingFFAGMagnet*>(getElement()->removeWrappers());
    magnet->setDipoleConstant(Attributes::getReal(itsAttr[B0]));
    magnet->setR0(Attributes::getReal(itsAttr[R0]));
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setTanDelta(Attributes::getReal(itsAttr[TAN_DELTA]));
    int maxOrder = floor(Attributes::getReal(itsAttr[MAX_Y_POWER]));
    magnet->setMaxOrder(maxOrder);
    endfieldmodel::Tanh* endField = dynamic_cast<endfieldmodel::Tanh*>(magnet->getEndField());
    endField->setLambda(Attributes::getReal(itsAttr[END_LENGTH]));
    endField->setX0(Attributes::getReal(itsAttr[CENTRE_LENGTH]));
    endField->setTanhDiffIndices(maxOrder+2);
    magnet->setRMin(Attributes::getReal(itsAttr[R0])-
                    Attributes::getReal(itsAttr[RADIAL_NEG_EXTENT]));
    magnet->setRMax(Attributes::getReal(itsAttr[R0])+
                    Attributes::getReal(itsAttr[RADIAL_POS_EXTENT]));
    Vector_t centre(-Attributes::getReal(itsAttr[R0]), 0, 0);
    magnet->setCentre(centre);
    // default length of the magnet element in radians
    double defaultLength = (endField->getLambda()*4.+endField->getX0());
    if (itsAttr[PHI_START]) {
        magnet->setPhiStart(Attributes::getReal(itsAttr[PHI_START]));
    } else {
        magnet->setPhiStart(defaultLength/2.);
    }
    if (itsAttr[PHI_END]) {
        magnet->setPhiEnd(Attributes::getReal(itsAttr[PHI_END]));
    } else {
        magnet->setPhiEnd(defaultLength);
    }
    if (itsAttr[AZIMUTHAL_EXTENT]) {
        magnet->setAzimuthalExtent(Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT]));
    } else {
        magnet->setAzimuthalExtent(Attributes::getReal(itsAttr[END_LENGTH])*5.+
                                   Attributes::getReal(itsAttr[CENTRE_LENGTH])*0.501);
    }
    magnet->initialise();
    setElement(magnet->makeAlignWrapper());

}
