// ------------------------------------------------------------------------
// $RCSfile: Ip.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BeamBeam.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

extern Inform* gmsg;

// Class BeamBeam
// ------------------------------------------------------------------------

BeamBeam::BeamBeam() : BeamBeam("") {}

BeamBeam::BeamBeam(const BeamBeam& right) : ElementBase(right), nSlices_m(right.nSlices_m) {}

BeamBeam::BeamBeam(const std::string& name) : ElementBase(name), nSlices_m(1) {}

BeamBeam::~BeamBeam() {}

void BeamBeam::accept(BeamlineVisitor& visitor) const { visitor.visitBeamBeam(*this); }

void BeamBeam::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

// set the number of slices for map tracking
void BeamBeam::setNSlices(const std::size_t& nSlices) { nSlices_m = nSlices; }

// get the number of slices for map tracking
std::size_t BeamBeam::getNSlices() const { return nSlices_m; }

void BeamBeam::finalise() {}

void BeamBeam::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

ElementType BeamBeam::getType() const { return ElementType::BEAMBEAM; }
