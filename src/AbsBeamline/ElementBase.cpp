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
#include "Utilities/Options.h"
#include "Utility/Inform.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <vector>

extern Inform* gmsg;

const std::vector<double> ElementBase::defaultAperture_m = std::vector<double>({1e6, 1e6});

const std::map<ElementType, std::string> ElementBase::elementTypeToString_s = {
        {ElementType::ANY, "Any"},
        {ElementType::CYCLOTRONSECTOR, "CyclotronSector"},
        {ElementType::BEAMLINE, "Beamline"},
        {ElementType::BOX, "Box"},
        {ElementType::COLLIMATOR, "Collimator"},
        {ElementType::DRIFT, "Drift"},
        {ElementType::FIELDMAP, "Fieldmap"},
        {ElementType::LASER, "Laser"},
        {ElementType::MARKER, "Marker"},
        {ElementType::MONITOR, "Monitor"},
        {ElementType::MULTIPOLE, "Multipole"},
        {ElementType::MULTIPOLET, "MULTIPOLET"},
        {ElementType::RFCAVITY, "RFCavity"},
        {ElementType::TRAVELINGWAVE, "TravelingWave"},
        {ElementType::SBEND, "SBEND"},
        {ElementType::SCALINGFFAMAGNET, "ScalingFFAMagnet"},
        {ElementType::RBEND, "RBEND"},
        {ElementType::RBEND3D, "RBEND3D"},
        {ElementType::SOURCE, "SOURCE"},
        {ElementType::SOLENOID, "SOLENOID"},
        {ElementType::PROBE, "Probe"},
        {ElementType::VACUUM, "Vacuum"},
        {ElementType::CONSTANTEFIELDCAVITY, "ConstantEFieldCavity"},
        {ElementType::CONSTANTFOCUSING, "ConstantFocusing"}};

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
      beamlineMembership_m(right.beamlineMembership_m),
      userAttribs(right.userAttribs),
      positionIsFixed(right.positionIsFixed),
      elementPosition_m(right.elementPosition_m),
      elemedgeSet_m(right.elemedgeSet_m),
      outputfn_m(right.outputfn_m),
      deleteOnTransverseExit_m(right.deleteOnTransverseExit_m),
      reach_m(right.reach_m) {}

ElementBase::ElementBase(const std::string& name)
    : shareFlag(true),
      csTrafoGlobal2Local_m(),
      misalignment_m(),
      rotationZAxis_m(0.0),
      RefPartBunch_m(nullptr),
      online_m(false),
      elementID(name),
      beamlineMembership_m(),
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

const BeamlineMembership& ElementBase::getBeamlineMembership() const {
    return beamlineMembership_m;
}

BeamlineTopology ElementBase::getBeamlineTopology() const { return beamlineMembership_m.topology; }

const std::string& ElementBase::getBeamlineOwnerName() const {
    return beamlineMembership_m.ownerName;
}

void ElementBase::setBeamlineMembership(BeamlineTopology topology, std::string ownerName) {
    if (topology == BeamlineTopology::RING && ownerName.empty()) {
        throw GeneralOpalException(
                "ElementBase::setBeamlineMembership()",
                "RING membership requires a non-empty owner name");
    }
    if (topology == BeamlineTopology::LINE && !ownerName.empty()) {
        throw GeneralOpalException(
                "ElementBase::setBeamlineMembership()",
                "LINE membership cannot have an owner name");
    }

    beamlineMembership_m = {topology, std::move(ownerName)};
}

void ElementBase::clearBeamlineMembership() { beamlineMembership_m = {}; }

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

namespace {
    /// An APERTURE half-width of this size or more counts as no aperture (the default is 1e6 m).
    constexpr double unboundedAperture = 1e5;
}  // namespace

std::optional<ApertureShape> ElementBase::getFieldRange() const { return std::nullopt; }

std::optional<ApertureShape> ElementBase::getFieldMapRange() const { return std::nullopt; }

bool ElementBase::hasApertureSet() const {
    return std::abs(aperture_m.second[0]) < unboundedAperture
           || std::abs(aperture_m.second[1]) < unboundedAperture;
}

ElementRegions ElementBase::getRegions() const {
    const Geometry& geometry = getGeometry();
    const ElementType type   = getType();
    // Sector bends and curved MULTIPOLET measure their windows along the design arc.
    const bool arc = (type == ElementType::SBEND || type == ElementType::MULTIPOLET)
                     && geometry.kind() == GeometryKind::SBend && geometry.getCurvature() != 0.0;

    ElementRegions regions;
    getFieldExtent(regions.fieldBegin, regions.fieldEnd);
    regions.curvature  = arc ? geometry.getCurvature() : 0.0;
    regions.bodyLength = geometry.getElementLength();
    regions.bodyBegin  = arc ? 0.0 : geometry.getStartZ();
    regions.bodyEnd    = regions.bodyBegin + geometry.getElementLength();

    const double reach = reach_m.value_or(Options::elementMappingReach);
    regions.reachX     = reach;
    regions.reachY     = reach;

    const double apertureX = std::abs(aperture_m.second[0]);
    const double apertureY = std::abs(aperture_m.second[1]);
    if (hasApertureSet()) {
        // The APERTURE is the bore; the element is material out to its reach.
        regions.aperture.type    = aperture_m.first;
        regions.aperture.halfX   = apertureX < unboundedAperture ? apertureX : reach;
        regions.aperture.halfY   = apertureY < unboundedAperture ? apertureY : reach;
        regions.deletesParticles = true;
    } else if (const auto mapRange = getFieldMapRange()) {
        // A field map without APERTURE: its extent is the whole element and nothing is material.
        regions.aperture = *mapRange;
        regions.reachX   = std::abs(mapRange->centreX) + mapRange->halfX;
        regions.reachY   = std::abs(mapRange->centreY) + mapRange->halfY;
    } else if (const auto range = getFieldRange()) {
        regions.aperture         = *range;
        regions.deletesParticles = true;
    } else {
        // No aperture: the field acts over the whole reach and nothing is deleted.
        regions.aperture.type  = ApertureType::RECTANGULAR;
        regions.aperture.halfX = reach;
        regions.aperture.halfY = reach;
    }
    return regions;
}

void ElementBase::checkRegions() const {
    const ElementRegions regions = getRegions();
    if (!regions.deletesParticles || regions.aperture.fitsInside(regions.reachX, regions.reachY))
        return;

    const auto& a = regions.aperture;
    std::ostringstream message;
    message << "Element '" << getName() << "': "
            << (hasApertureSet() ? "its APERTURE"
                                 : "the transverse range of its field (HAPERT/VAPERT or "
                                   "WIDTH/HEIGHT), which is its aperture without an APERTURE,")
            << " reaches |x| = " << std::abs(a.centreX) + a.halfX
            << " m, |y| = " << std::abs(a.centreY) + a.halfY << " m, beyond its reach of "
            << regions.reachX << " m. Set REACH on the element or raise OPTION "
            << "ELEMENTMAPPINGREACH.";
    throw GeneralOpalException("ElementBase::checkRegions", message.str());
}

std::pair<double, double> ElementBase::getSupportHalfWidths(const double reach) const {
    const auto halfWidth = [reach](const double aperture) {
        return std::abs(aperture) >= unboundedAperture ? reach
                                                       : std::max(std::abs(aperture), reach);
    };
    return {halfWidth(aperture_m.second[0]), halfWidth(aperture_m.second[1])};
}

BoundingBox ElementBase::getLocalSupportBox(const double reach) const {
    const Geometry& geometry = getGeometry();
    double zBegin            = 0.0;
    double zEnd              = 0.0;
    getFieldExtent(zBegin, zEnd);

    // A field map without APERTURE is just its map: the map's extent across, its field window
    // along z.
    if (!hasApertureSet()) {
        if (const auto map = getFieldMapRange()) {
            return BoundingBox::getBoundingBox(
                    {Vector_t<double, 3>(
                             {map->centreX - map->halfX, map->centreY - map->halfY, zBegin}),
                     Vector_t<double, 3>(
                             {map->centreX + map->halfX, map->centreY + map->halfY, zEnd})});
        }
    }

    const auto [halfX, halfY] = getSupportHalfWidths(reach);
    zBegin                    = std::min(zBegin, geometry.getStartZ());
    zEnd                      = std::max(zEnd, geometry.getStartZ() + geometry.getElementLength());

    if (geometry.kind() == GeometryKind::SBend && geometry.getCurvature() != 0.0) {
        return getArcSupportBox(
                geometry.getCurvature(), geometry.getElementLength(), zBegin, zEnd, halfX, halfY);
    }
    return BoundingBox::getBoundingBox(
            {Vector_t<double, 3>({-halfX, -halfY, zBegin}),
             Vector_t<double, 3>({halfX, halfY, zEnd})});
}

std::vector<SupportBox> ElementBase::getLocalSupportBoxes(const double reach) const {
    return {toSupportBox(getLocalSupportBox(reach))};
}

SupportBox ElementBase::toSupportBox(const BoundingBox& box) {
    const auto [lower, upper] = box.getCorners();
    SupportBox result;
    result.centre = 0.5 * (lower + upper);
    result.half   = 0.5 * (upper - lower);
    for (unsigned i = 0; i < 3; ++i) {
        result.axes[i]    = Vector_t<double, 3>(0.0);
        result.axes[i](i) = 1.0;
    }
    return result;
}

std::vector<SupportBox> ElementBase::getArcSupportBoxes(
        const double curvature, const double bodyLength, const double sBegin, const double sEnd,
        const double halfX, const double halfY) {
    const Vector_t<double, 3> vertical({0.0, 1.0, 0.0});
    const auto normalAt = [curvature](const double s) {
        return Vector_t<double, 3>({std::cos(curvature * s), 0.0, std::sin(curvature * s)});
    };
    const auto tangentAt = [curvature](const double s) {
        return Vector_t<double, 3>({-std::sin(curvature * s), 0.0, std::cos(curvature * s)});
    };
    // The design orbit is p(s) = C + n(s) / h around the centre of curvature C = (-1/h, 0, 0).
    const Vector_t<double, 3> centreOfCurvature({-1.0 / curvature, 0.0, 0.0});

    std::vector<SupportBox> boxes;
    const auto addStraight = [&](const Vector_t<double, 3>& start,
                                 const Vector_t<double, 3>& normal,
                                 const Vector_t<double, 3>& tangent, const double length) {
        boxes.push_back(
                {start + 0.5 * length * tangent,
                 {normal, vertical, tangent},
                 Vector_t<double, 3>({halfX, halfY, 0.5 * length})});
    };

    if (sBegin < 0.0) {
        addStraight(
                Vector_t<double, 3>({0.0, 0.0, sBegin}), Vector_t<double, 3>({1.0, 0.0, 0.0}),
                Vector_t<double, 3>({0.0, 0.0, 1.0}), -sBegin);
    }

    const double first = std::clamp(sBegin, 0.0, bodyLength);
    const double last  = std::clamp(sEnd, 0.0, bodyLength);
    if (last > first) {
        // A point of the band is C + a (cos(delta) n_m + sin(delta) t_m) with a = 1/h + x,
        // |x| <= halfX, and |delta| <= half the segment angle around its middle direction n_m.
        constexpr double maximumSegmentAngle = 0.1;
        const int segments                   = std::clamp(
                static_cast<int>(
                        std::ceil(std::abs(curvature) * (last - first) / maximumSegmentAngle)),
                1, 256);
        const double aLow  = 1.0 / curvature - halfX;
        const double aHigh = 1.0 / curvature + halfX;
        for (int i = 0; i < segments; ++i) {
            const double s0       = first + (last - first) * i / segments;
            const double s1       = first + (last - first) * (i + 1) / segments;
            const double middle   = 0.5 * (s0 + s1);
            const double halfTurn = 0.5 * std::abs(curvature * (s1 - s0));
            const double c        = std::cos(halfTurn);
            const double along[4] = {aLow * c, aLow, aHigh * c, aHigh};
            const double low      = *std::min_element(along, along + 4);
            const double high     = *std::max_element(along, along + 4);
            const double across   = std::max(std::abs(aLow), std::abs(aHigh)) * std::sin(halfTurn);
            boxes.push_back(
                    {centreOfCurvature + 0.5 * (low + high) * normalAt(middle),
                     {normalAt(middle), vertical, tangentAt(middle)},
                     Vector_t<double, 3>({0.5 * (high - low), halfY, across})});
        }
    }

    if (sEnd > bodyLength) {
        addStraight(
                centreOfCurvature + normalAt(bodyLength) / curvature, normalAt(bodyLength),
                tangentAt(bodyLength), sEnd - bodyLength);
    }
    return boxes;
}

BoundingBox ElementBase::getArcSupportBox(
        const double curvature, const double bodyLength, const double sBegin, const double sEnd,
        const double halfX, const double halfY) {
    BoundingBox box;
    const auto addSection = [&](const Vector_t<double, 3>& centre,
                                const Vector_t<double, 3>& normal) {
        for (const double x : {-halfX, halfX}) {
            for (const double y : {-halfY, halfY}) {
                box.enlargeToContainPosition(centre + x * normal + Vector_t<double, 3>({0, y, 0}));
            }
        }
    };
    const auto centreAt = [curvature](const double s) {
        const double phi = curvature * s;
        return Vector_t<double, 3>(
                {(std::cos(phi) - 1.0) / curvature, 0.0, std::sin(phi) / curvature});
    };
    const auto normalAt = [curvature](const double s) {
        const double phi = curvature * s;
        return Vector_t<double, 3>({std::cos(phi), 0.0, std::sin(phi)});
    };

    // Upstream of the entrance the chart is the identity; the section at s = 0 closes it.
    if (sBegin < 0.0) {
        addSection(Vector_t<double, 3>({0.0, 0.0, sBegin}), Vector_t<double, 3>({1.0, 0.0, 0.0}));
    }

    // Sample the curved body finely enough that the band bulges by less than a micrometre per
    // metre of radius between samples, then add that bulge as a margin.
    const double first     = std::clamp(sBegin, 0.0, bodyLength);
    const double last      = std::clamp(sEnd, 0.0, bodyLength);
    const double turn      = std::abs(curvature * (last - first));
    const int samples      = std::max(2, static_cast<int>(std::ceil(turn / 1e-3)) + 1);
    const double stepAngle = turn / (samples - 1);
    for (int i = 0; i < samples; ++i) {
        const double s = first + (last - first) * i / (samples - 1);
        addSection(centreAt(s), normalAt(s));
    }

    // Downstream of the exit the chart continues along the exit tangent.
    if (sEnd > bodyLength) {
        const double phi = curvature * bodyLength;
        const Vector_t<double, 3> tangent({-std::sin(phi), 0.0, std::cos(phi)});
        addSection(centreAt(bodyLength) + (sEnd - bodyLength) * tangent, normalAt(bodyLength));
    }

    const double bulge  = (1.0 / std::abs(curvature) + halfX) * (1.0 - std::cos(0.5 * stepAngle));
    auto [lower, upper] = box.getCorners();
    lower -= Vector_t<double, 3>({bulge, 0.0, bulge});
    upper += Vector_t<double, 3>({bulge, 0.0, bulge});
    return BoundingBox::getBoundingBox({lower, upper});
}

/* ============================== Field / physics interface ================== */

ElementType ElementBase::getType() const { return ElementType::ANY; }

void ElementBase::goOffline() { online_m = false; }

void ElementBase::apply(
        const std::shared_ptr<ParticleContainer_t>& /*pc*/, const StepTime& /*time*/) {}

void ElementBase::takeDesignParticle(const PartBunch_t* bunch) {
    if (bunch == nullptr) return;
    const auto pc = bunch->getParticleContainer();
    if (!pc || pc->getReference() == nullptr) return;
    designMass_m = pc->getReference()->getM();
    if (std::abs(pc->getReference()->getQ()) > 1.0e-15) designCharge_m = pc->getReference()->getQ();
}

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

    // Inside the body window and the reach, outside the aperture: the particle is in material.
    const ElementRegions regions = getRegions();
    auto Rview                   = pc->R.getView();
    auto invalid                 = pc->InvalidMask.getView();

    size_t localMarked = 0;
    Kokkos::parallel_reduce(
            "ElementBase::markOutsideAperture", nLocal,
            KOKKOS_LAMBDA(const size_t i, size_t& count) {
                const bool hit         = regions.inLossRegion(Rview(i));
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
    return getRegions().inLossRegion(R);
}
/* ========================================================================== */
