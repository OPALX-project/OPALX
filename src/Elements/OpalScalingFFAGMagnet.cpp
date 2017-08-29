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
                              ("B0", "The nominal dipole field of the magnet [T].");
    itsAttr[R0] = Attributes::makeReal("R0", "Radial scale [m].");
    itsAttr[FIELD_INDEX] = Attributes::makeReal("FIELD_INDEX",
      "The scaling magnet field index.");
    itsAttr[TAN_DELTA] = Attributes::makeReal("TAN_DELTA",
      "Tangent of the spiral angle; set to 0 to make a radial sector magnet.");
    itsAttr[MAX_Y_POWER] = Attributes::makeReal("MAX_Y_POWER",
      "The maximum power in y that will be considered in the field expansion.");
    itsAttr[END_LENGTH] = Attributes::makeReal("END_LENGTH",
                                          "The end length of the spiral FFAG [m].");
    itsAttr[HEIGHT] = Attributes::makeReal("HEIGHT",
                                       "Full height of the magnet. Particles moving more than height/2. off the midplane (either above or below) are out of the aperture [m].");
    itsAttr[CENTRE_LENGTH] = Attributes::makeReal("CENTRE_LENGTH",
                                       "The centre length of the spiral FFAG [m].");
    itsAttr[RADIAL_NEG_EXTENT] = Attributes::makeReal("RADIAL_NEG_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0-RADIAL_NEG_EXTENT [m].");
    itsAttr[RADIAL_POS_EXTENT] = Attributes::makeReal("RADIAL_POS_EXTENT",
                                       "Particles are considered outside the tracking region if radius is greater than R0+RADIAL_POS_EXTENT [m].");
    itsAttr[MAGNET_START] = Attributes::makeReal("MAGNET_START",
                                          "Offset the centre of the spiral FFAG magnet (default is centre_length/2 + end_length). [m]");
    itsAttr[MAGNET_END] = Attributes::makeReal("MAGNET_END",
                                       "Offset to the end of the magnet, i.e. placement of the next element. Default is centre_length + 2*end_length.");
    itsAttr[AZIMUTHAL_EXTENT] = Attributes::makeReal("AZIMUTHAL_EXTENT",
                                       "The field will be assumed zero if particles are more than AZIMUTHAL_EXTENT from the magnet centre (psi=0). Default is CENTRE_LENGTH/2.+5.*END_LENGTH [m].");
    registerRealAttribute("B0");
    registerRealAttribute("R0");
    registerRealAttribute("FIELD_INDEX");
    registerRealAttribute("TAN_DELTA");
    registerRealAttribute("MAX_Y_POWER");
    registerRealAttribute("END_LENGTH");
    registerRealAttribute("CENTRE_LENGTH");
    registerRealAttribute("RADIAL_NEG_EXTENT");
    registerRealAttribute("RADIAL_POS_EXTENT");
    registerRealAttribute("HEIGHT");
    registerRealAttribute("MAGNET_START");
    registerRealAttribute("MAGNET_END");
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

    // use L = r0*theta; we define the magnet ito length for UI but ito angles
    // internally; and use m as external default unit and mm internally
    double metres = 1e3;
    // get r0 in m
    double r0 = Attributes::getReal(itsAttr[R0]);
    magnet->setR0(r0*metres);
    // get B0 in T
    magnet->setDipoleConstant(Attributes::getReal(itsAttr[B0]));

    // dimensionless quantities
    magnet->setFieldIndex(Attributes::getReal(itsAttr[FIELD_INDEX]));
    magnet->setTanDelta(Attributes::getReal(itsAttr[TAN_DELTA]));
    int maxOrder = floor(Attributes::getReal(itsAttr[MAX_Y_POWER]));
    magnet->setMaxOrder(maxOrder);

    // get centre length and end length in radians
    endfieldmodel::Tanh* endField = dynamic_cast<endfieldmodel::Tanh*>(magnet->getEndField());
    double end_length = Attributes::getReal(itsAttr[END_LENGTH])/r0;
    double centre_length = Attributes::getReal(itsAttr[CENTRE_LENGTH])/r0;
    endField->setLambda(end_length);
    endField->setX0(centre_length);
    endField->setTanhDiffIndices(maxOrder+2);

    // get rmin and rmax bounding box edge in mm
    double rmin = r0-Attributes::getReal(itsAttr[RADIAL_NEG_EXTENT]);
    double rmax = r0+Attributes::getReal(itsAttr[RADIAL_POS_EXTENT]);
    magnet->setRMin(rmin*metres);
    magnet->setRMax(rmax*metres);
    Vector_t centre(-r0*metres, 0, 0);
    magnet->setCentre(centre);

    // we store maximum vertical displacement (which is half the height)
    double height = Attributes::getReal(itsAttr[HEIGHT])*metres;
    magnet->setVerticalExtent(height/2.);
    std::cerr << "OpalScalingFFAGMagnet::update height " << height << std::endl;

    // get default length of the magnet element in radians
    double defaultLength = (endField->getLambda()*4.+endField->getX0());

    // get start of the magnet element in radians
    if (itsAttr[MAGNET_START]) {
        double phi_start = Attributes::getReal(itsAttr[MAGNET_START])/r0;
        magnet->setPhiStart(phi_start);
    } else {
        magnet->setPhiStart(defaultLength/2.);
    }

    // get end of the magnet element in radians
    if (itsAttr[MAGNET_END]) {
        double phi_end = Attributes::getReal(itsAttr[MAGNET_END])/r0;
        magnet->setPhiEnd(phi_end);
    } else {
        magnet->setPhiEnd(defaultLength);
    }
    // get azimuthal extent in radians
    if (itsAttr[AZIMUTHAL_EXTENT]) {
        magnet->setAzimuthalExtent(Attributes::getReal(itsAttr[AZIMUTHAL_EXTENT])/r0);
    } else {
        magnet->setAzimuthalExtent(end_length*5.+centre_length*0.501);
    }
    magnet->initialise();
    setElement(magnet->makeAlignWrapper());

}
