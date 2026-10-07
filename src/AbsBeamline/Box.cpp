//
// Class Box
//   Interface for a rectangular block that absorbs every particle inside it.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/Box.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

Box::Box() : Box("") {}

Box::Box(const Box& right)
    : ElementBase(right), halfWidth_m(right.halfWidth_m), halfHeight_m(right.halfHeight_m) {}

Box::Box(const std::string& name) : ElementBase(name), halfWidth_m(0.0), halfHeight_m(0.0) {}

Box::~Box() {}

void Box::accept(BeamlineVisitor& visitor) const { visitor.visitBox(*this); }

void Box::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

void Box::finalise() {}

void Box::getFieldExtent(double& zBegin, double& zEnd) const {
    // Local-chart support interval (field-free; report the body span).
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

ElementType Box::getType() const { return ElementType::BOX; }

size_t Box::markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) {
    if (!pc || !getFlagDeleteOnTransverseExit()) {
        return 0;
    }
    const size_t nLocal = pc->getLocalNum();
    if (nLocal == 0) {
        return 0;
    }

    // Same body window as ElementBase::markOutsideAperture.
    const double zBegin = getGeometry().getStartZ();
    const double zEnd   = zBegin + getGeometry().getElementLength();

    // Members copied to locals; the device kernel must not capture `this`.
    const double xLimit = halfWidth_m;
    const double yLimit = halfHeight_m;

    auto Rview   = pc->R.getView();
    auto invalid = pc->InvalidMask.getView();

    size_t localMarked = 0;
    Kokkos::parallel_reduce(
            "Box::markOutsideAperture", nLocal,
            KOKKOS_LAMBDA(const size_t i, size_t& count) {
                const bool hit = Rview(i)[2] >= zBegin && Rview(i)[2] < zEnd
                                 && Kokkos::fabs(Rview(i)[0]) < xLimit
                                 && Kokkos::fabs(Rview(i)[1]) < yLimit;
                const bool newlyMarked = hit && !invalid(i);
                invalid(i)             = invalid(i) || hit;
                count += newlyMarked ? 1 : 0;
            },
            localMarked);
    Kokkos::fence();

    return localMarked;
}
