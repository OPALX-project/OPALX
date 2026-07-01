// ------------------------------------------------------------------------
// $RCSfile: Drift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
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

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

extern Inform* gmsg;

// Class Drift
// ------------------------------------------------------------------------

Drift::Drift() : Drift("") {}

Drift::Drift(const Drift& right) : ElementBase(right), nSlices_m(right.nSlices_m) {}

Drift::Drift(const std::string& name) : ElementBase(name), nSlices_m(1) {}

Drift::~Drift() {}

void Drift::accept(BeamlineVisitor& visitor) const { visitor.visitDrift(*this); }

void Drift::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

// set the number of slices for map tracking
void Drift::setNSlices(const std::size_t& nSlices) { nSlices_m = nSlices; }

// get the number of slices for map tracking
std::size_t Drift::getNSlices() const { return nSlices_m; }

void Drift::finalise() {}

bool Drift::bends() const { return false; }

void Drift::getFieldExtent(double& zBegin, double& zEnd) const {
    // Local-chart field-support interval (a drift carries no field; report the body span).
    zBegin = 0.0;
    zEnd   = getElementLength();
}

ElementType Drift::getType() const { return ElementType::DRIFT; }
