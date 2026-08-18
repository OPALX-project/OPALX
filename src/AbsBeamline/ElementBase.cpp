// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "AbsBeamline/ElementBase.h"

#include "Channels/Channel.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/LogicalError.h"
#include "Utility/Inform.h"

#include <filesystem>
#include <iostream>
#include <vector>

extern Inform* gmsg;

const std::vector<double> ElementBase::defaultAperture_m = std::vector<double>({1e6, 1e6});

const std::map<ElementType, std::string> ElementBase::elementTypeToString_s = {
        {ElementType::ANY, "Any"},
        {ElementType::BEAMLINE, "Beamline"},
        {ElementType::COLLIMATOR, "Collimator"},
        {ElementType::DRIFT, "Drift"},
        {ElementType::LASER, "Laser"},
        {ElementType::MARKER, "Marker"},
        {ElementType::MONITOR, "Monitor"},
        {ElementType::MULTIPOLE, "Multipole"},
        {ElementType::RFCAVITY, "RFCavity"},
        {ElementType::TRAVELINGWAVE, "TravelingWave"},
        {ElementType::SBEND, "SBEND"},
        {ElementType::RBEND, "RBEND"},
        {ElementType::RBEND3D, "RBEND3D"},
        {ElementType::RING, "Ring"},
        {ElementType::SOURCE, "SOURCE"},
        {ElementType::SOLENOID, "SOLENOID"},
        {ElementType::PROBE, "Probe"},
        {ElementType::VACUUM, "Vacuum"},
        {ElementType::CONSTANTEFIELDCAVITY, "ConstantEFieldCavity"}};

ElementBase::ElementBase() : ElementBase("") {}

ElementBase::ElementBase(const ElementBase& right)
    : std::enable_shared_from_this<ElementBase>(),
      shareFlag(true),
      csTrafoGlobal2Local_m(right.csTrafoGlobal2Local_m),
      misalignment_m(right.misalignment_m),
      rotationZAxis_m(right.rotationZAxis_m),
      aperture_m(right.aperture_m),
      RefPartBunch_m(nullptr),
      online_m(right.online_m),
      elementID(right.elementID),
      userAttribs(right.userAttribs),
      positionIsFixed(right.positionIsFixed),
      elementPosition_m(right.elementPosition_m),
      elemedgeSet_m(right.elemedgeSet_m),
      outputfn_m(right.outputfn_m),
      deleteOnTransverseExit_m(right.deleteOnTransverseExit_m) {}

ElementBase::ElementBase(const std::string& name)
    : shareFlag(true),
      csTrafoGlobal2Local_m(),
      misalignment_m(),
      rotationZAxis_m(0.0),
      RefPartBunch_m(nullptr),
      online_m(false),
      elementID(name),
      userAttribs(),
      positionIsFixed(false),
      elementPosition_m(0.0),
      elemedgeSet_m(false),
      outputfn_m("") {
    setAperture(ApertureType::ELLIPTICAL, defaultAperture_m);
}

ElementBase::~ElementBase() {}

const std::string& ElementBase::getName() const { return elementID; }

void ElementBase::setName(const std::string& name) { elementID = name; }

void ElementBase::setOutputFN(const std::string fn) { outputfn_m = fn; }

std::string ElementBase::getOutputFN() const {
    if (outputfn_m.empty()) {
        return getName();
    } else {
        std::filesystem::path filePath(outputfn_m);
        return filePath.replace_extension().native();
    }
}

double ElementBase::getAttribute(const std::string& aKey) const {
    const ConstChannel* aChannel = getConstChannel(aKey);

    if (aChannel != nullptr) {
        double val = *aChannel;
        delete aChannel;
        return val;
    } else {
        return 0.0;
    }
}

bool ElementBase::hasAttribute(const std::string& aKey) const {
    const ConstChannel* aChannel = getConstChannel(aKey);

    if (aChannel != nullptr) {
        delete aChannel;
        return true;
    } else {
        return false;
    }
}

void ElementBase::removeAttribute(const std::string& aKey) { userAttribs.removeAttribute(aKey); }

void ElementBase::setAttribute(const std::string& aKey, double val) {
    Channel* aChannel = getChannel(aKey, true);

    if (aChannel != nullptr && aChannel->isSettable()) {
        *aChannel = val;
        delete aChannel;
    } else
        std::cout << "Channel nullptr or not Settable" << std::endl;
}

Channel* ElementBase::getChannel(const std::string& aKey, bool create) {
    return userAttribs.getChannel(aKey, create);
}

const ConstChannel* ElementBase::getConstChannel(const std::string& aKey) const {
    // Use const_cast to allow calling the non-const method GetChannel().
    // The const return value of this method will nevertheless inhibit set().
    return const_cast<ElementBase*>(this)->getChannel(aKey);
}

std::string ElementBase::getTypeString(ElementType type) { return elementTypeToString_s.at(type); }

ElementBase* ElementBase::copyStructure() {
    if (isSharable()) {
        return this;
    } else {
        return clone();
    }
}

void ElementBase::makeSharable() { shareFlag = true; }

bool ElementBase::update(const AttributeSet& set) {
    for (AttributeSet::const_iterator i = set.begin(); i != set.end(); ++i) {
        setAttribute(i->first, i->second);
    }

    return true;
}

BoundingBox ElementBase::getBoundingBoxInLabCoords() const {
    CoordinateSystemTrafo toBegin = getGeometry().getEdgeToBegin() * csTrafoGlobal2Local_m;
    CoordinateSystemTrafo toEnd   = getGeometry().getEdgeToEnd() * csTrafoGlobal2Local_m;

    const double& x = aperture_m.second[0];
    const double& y = aperture_m.second[1];

    std::vector<Vector_t<double, 3>> corners(8);
    for (int i = -1; i < 2; i += 2) {
        for (int j = -1; j < 2; j += 2) {
            unsigned int idx = (i + 1) / 2 + (j + 1);
            corners[idx]     = toBegin.transformFrom(Vector_t<double, 3>({i * x, j * y, 0.0}));
            corners[idx + 4] = toEnd.transformFrom(Vector_t<double, 3>({i * x, j * y, 0.0}));
        }
    }

    return BoundingBox::getBoundingBox(corners);
}

/* ============================== Field / physics interface ================== */

ElementType ElementBase::getType() const { return ElementType::ANY; }

void ElementBase::goOnline(const double&) { online_m = true; }

void ElementBase::goOffline() { online_m = false; }

void ElementBase::apply(const std::shared_ptr<ParticleContainer_t>& /*pc*/) {}

void ElementBase::apply(
        const Vector_t<double, 3>& /*R*/, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    // A bare element carries no field. The bounds/aperture test the previous
    // bool overload performed here only produced the (unused) return value.
}

size_t ElementBase::markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) {
    if (!pc || !getFlagDeleteOnTransverseExit()) {
        return 0;
    }
    const size_t nLocal = pc->getLocalNum();
    if (nLocal == 0) {
        return 0;
    }

    // The aperture is a geometric property of the element body, so gate on the
    // geometric extent [0, L] (element-local frame), consistent with
    // applyToReferenceParticle() below. The field-support window
    // (getFieldExtent) can be narrower or offset from the body -- e.g. Solenoid
    // returns its field-map range and Monitor a plane-centered window -- which
    // would leave part of the body unchecked.
    const double zBegin = 0.0;
    const double zEnd   = getGeometry().getElementLength();

    // Members copied to locals; the device kernel must not capture `this`.
    const ApertureType type = aperture_m.first;
    const double xLimit     = aperture_m.second[0];
    const double yLimit     = aperture_m.second[1];

    auto Rview   = pc->R.getView();
    auto invalid = pc->InvalidMask.getView();

    size_t localMarked = 0;
    Kokkos::parallel_reduce(
            "ElementBase::markOutsideAperture", nLocal,
            KOKKOS_LAMBDA(const size_t i, size_t& count) {
                // z-window [zBegin, zEnd), matching isInside; z is measured in
                // the element-local frame.
                const bool inZ = Rview(i)[2] >= zBegin && Rview(i)[2] < zEnd;
                const bool hit = inZ
                                 && !ApertureHelper::isInsideAperture(
                                         Rview(i)[0], Rview(i)[1], type, xLimit, yLimit);
                const bool newlyMarked = hit && !invalid(i);
                invalid(i)             = invalid(i) || hit;
                count += newlyMarked ? 1 : 0;
            },
            localMarked);
    Kokkos::fence();

    return localMarked;
}

bool ElementBase::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    if (R(2) >= 0.0 && R(2) < getGeometry().getElementLength()) {
        if (!ApertureHelper::isInsideAperture(R, aperture_m)) {
            return true;
        }
    }
    return false;
}
/* ========================================================================== */
