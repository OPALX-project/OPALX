//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "AbsBeamline/MultipoleT.h"

#include "AbsBeamline/BeamlineVisitor.h"

#include "BeamlineGeometry/Geometry.h"

#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>

MultipoleT::MultipoleT() : MultipoleT("") {}

MultipoleT::MultipoleT(const std::string& name)
    : ElementBase(name), fringeLeft_m(0.0), fringeRight_m(0.0), maxFOrder_m(3) {
    rebuildTanhTable();
}

MultipoleT::MultipoleT(const MultipoleT& right)
    : ElementBase(right),
      profile_m(right.profile_m),
      tanhTable_m(right.tanhTable_m),
      tanhTableHost_m(right.tanhTableHost_m),
      fringeLeft_m(right.fringeLeft_m),
      fringeRight_m(right.fringeRight_m),
      maxFOrder_m(right.maxFOrder_m),
      scalingName_m(right.scalingName_m),
      scalingTD_m(right.scalingTD_m) {}

MultipoleT::~MultipoleT() = default;

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    initialiseTimeDependencies();
    visitor.visitMultipoleT(*this);
}

ElementType MultipoleT::getType() const {
    return ElementType::MULTIPOLET;
}

void MultipoleT::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    online_m       = true;
}

void MultipoleT::finalise() {
    online_m = false;
}

void MultipoleT::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    auto Rview          = pc->R.getView();
    auto Bview          = pc->B.getView();
    const size_t nLocal = pc->getLocalNum();

    // Field-support extent (the same source isInside() uses), read on the host: getFieldExtent
    // is not device callable.
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);

    // Profile, geometry and fringe lengths built once on the host and captured by value; the
    // kernel must not capture `this` or the Geometry.
    const double time                                = (RefPartBunch_m != nullptr)
                                                               ? RefPartBunch_m->getT()
                                                               : 0.0;
    const MultipoleTFieldModel::FieldInputs inputs   = makeFieldInputs(time);
    const Kokkos::View<double**> tanhTable           = tanhTable_m;

    Kokkos::parallel_for(
            "MultipoleT::apply", nLocal, KOKKOS_LAMBDA(const size_t i) {
                // Convert (x,y,z) -> (radial x, y, arc s) in the element's frame.
                const Vector_t<double, 3> arc = GeometryHelper::toBendArcCoords(
                        Rview(i), inputs.curvature, inputs.bodyLength);

                if (arc(2) < zBegin || arc(2) > zEnd) {
                    return;  // return this particle's lambda
                }

                // The model works in the local basis at that arc position, so rotate its field
                // into the entrance frame the tracker hands us.
                const Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
                        MultipoleTFieldModel::field(arc, inputs, tanhTable), arc(2),
                        inputs.curvature, inputs.bodyLength);

                Bview(i)(0) += Bf(0);
                Bview(i)(1) += Bf(1);
                Bview(i)(2) += Bf(2);
            });
}

void MultipoleT::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const Vector_t<double, 3> arc = bendCoords(R);
    if (!isInsideArc(arc)) {
        return;
    }
    if (!ApertureHelper::isInsideAperture(arc, aperture_m)) {
        return;
    }

    computeFieldHost(R, t, B);
    (void)E;
}

bool MultipoleT::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>&, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const Vector_t<double, 3> arc = bendCoords(R);
    if (!isInsideArc(arc)) {
        return false;
    }
    if (!ApertureHelper::isInsideAperture(arc, aperture_m)) {
        return true;
    }

    computeFieldHost(R, t, B);
    (void)E;
    return false;
}

size_t MultipoleT::markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) {
    if (!pc || !getFlagDeleteOnTransverseExit()) {
        return 0;
    }
    const size_t nLocal = pc->getLocalNum();
    if (nLocal == 0) {
        return 0;
    }

    // Members copied to locals; the device kernel must not capture `this`.
    const ApertureType type = aperture_m.first;
    const double xLimit     = aperture_m.second[0];
    const double yLimit     = aperture_m.second[1];
    const double curvature  = getGeometry().getCurvature();
    const double bodyLength = getGeometry().getElementLength();

    // Geometric body window [0, body length), not the field extent: the aperture is hardware.
    const double zBegin = 0.0;
    const double zEnd   = bodyLength;

    auto Rview   = pc->R.getView();
    auto invalid = pc->InvalidMask.getView();

    size_t localMarked = 0;
    Kokkos::parallel_reduce(
            "MultipoleT::markOutsideAperture", nLocal,
            KOKKOS_LAMBDA(const size_t i, size_t& count) {
                // Convert (x,y,z) -> (radial x, y, arc s): the z-window and the aperture are
                // measured relative to the design orbit, matching isInside.
                const Vector_t<double, 3> arc =
                        GeometryHelper::toBendArcCoords(Rview(i), curvature, bodyLength);

                const bool inZ = arc(2) >= zBegin && arc(2) < zEnd;
                const bool hit =
                        inZ
                        && !ApertureHelper::isInsideAperture(arc(0), arc(1), type, xLimit, yLimit);
                const bool newlyMarked = hit && !invalid(i);
                invalid(i)             = invalid(i) || hit;
                count += newlyMarked ? 1 : 0;
            },
            localMarked);
    Kokkos::fence();

    return localMarked;
}

void MultipoleT::getFieldExtent(double& zBegin, double& zEnd) const {
    // The tanh profile has fallen to about 6e-6 of the flat top at FringeReach fringe lengths.
    // A zero fringe length is a hard edge, so this is then the plain body extent [0, L].
    zBegin = -MultipoleTFieldModel::FringeReach * fringeLeft_m;
    zEnd   = getGeometry().getElementLength()
           + MultipoleTFieldModel::FringeReach * fringeRight_m;
}

bool MultipoleT::isInside(const Vector_t<double, 3>& r) const {
    // Selection uses the arc length s and the radial offset (not the straight-frame z/x), so a
    // curved body stays selected as the orbit turns through it and the aperture is measured
    // relative to the design orbit.
    const Vector_t<double, 3> arc = bendCoords(r);
    return isInsideArc(arc) && ApertureHelper::isInsideAperture(arc, aperture_m);
}

void MultipoleT::setTransverseProfile(const std::vector<double>& profile) {
    if (profile.size() > MultipoleTFieldModel::NumPoles) {
        throw OpalException(
                "MultipoleT::setTransverseProfile",
                "The transverse profile takes at most "
                        + std::to_string(MultipoleTFieldModel::NumPoles)
                        + " coefficients (dipole to decapole).");
    }
    for (unsigned int i = 0; i < MultipoleTFieldModel::NumPoles; ++i) {
        profile_m[i] = (i < profile.size()) ? profile[i] : 0.0;
    }
}

void MultipoleT::setFringeLengths(const double left, const double right) {
    if (left < 0.0 || right < 0.0) {
        throw OpalException(
                "MultipoleT::setFringeLengths",
                "The fringe lengths must not be negative (zero is a hard edge).");
    }
    fringeLeft_m  = left;
    fringeRight_m = right;
}

void MultipoleT::setMaxFOrder(const unsigned int order) {
    if (order < 1 || order > MultipoleTFieldModel::MaxFOrder) {
        throw OpalException(
                "MultipoleT::setMaxFOrder",
                "The number of terms in the field expansion must be between 1 and "
                        + std::to_string(MultipoleTFieldModel::MaxFOrder) + ".");
    }
    if (order == maxFOrder_m && tanhTable_m.extent(0) > 0) {
        return;
    }
    maxFOrder_m = order;
    rebuildTanhTable();
}

void MultipoleT::setScalingName(const std::string& name) {
    // Element names are stored in upper case.
    scalingName_m = name;
    std::ranges::transform(scalingName_m, scalingName_m.begin(), [](const unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
}

void MultipoleT::computeFieldHost(
        const Vector_t<double, 3>& R, const double t, Vector_t<double, 3>& B) const {
    const MultipoleTFieldModel::FieldInputs inputs = makeFieldInputs(t);
    const Vector_t<double, 3> arc =
            GeometryHelper::toBendArcCoords(R, inputs.curvature, inputs.bodyLength);
    const Vector_t<double, 3> Bf = GeometryHelper::rotateArcFieldToEntry(
            MultipoleTFieldModel::field(arc, inputs, tanhTableHost_m), arc(2), inputs.curvature,
            inputs.bodyLength);
    for (unsigned int d = 0; d < 3; ++d) {
        B(d) += Bf(d);
    }
}

MultipoleTFieldModel::FieldInputs MultipoleT::makeFieldInputs(const double t) const {
    MultipoleTFieldModel::FieldInputs in{};
    // TP holds the physical mid-plane field, so the only scaling is the time dependence.
    const double scaling = getScaling(t);
    for (unsigned int i = 0; i < MultipoleTFieldModel::NumPoles; ++i) {
        in.profile[i] = profile_m[i] * scaling;
    }
    in.maxFOrder   = maxFOrder_m;
    in.bodyLength  = getGeometry().getElementLength();
    in.curvature   = getGeometry().getCurvature();
    in.fringeLeft  = fringeLeft_m;
    in.fringeRight = fringeRight_m;
    return in;
}

Vector_t<double, 3> MultipoleT::bendCoords(const Vector_t<double, 3>& r) const {
    // The stored frame is the design-orbit entrance tangent, so the arc coordinate is measured
    // directly (no pole-face de-tilt).
    return GeometryHelper::toBendArcCoords(
            r, getGeometry().getCurvature(), getGeometry().getElementLength());
}

bool MultipoleT::isInsideArc(const Vector_t<double, 3>& arc) const {
    double zBegin = 0.0;
    double zEnd   = 0.0;
    getFieldExtent(zBegin, zEnd);
    return arc(2) >= zBegin && arc(2) < zEnd;
}

double MultipoleT::getScaling(const double t) const {
    return scalingTD_m ? scalingTD_m->getValue(t) : 1.0;
}

void MultipoleT::initialiseTimeDependencies() const {
    scalingTD_m.reset();
    if (!scalingName_m.empty()) {
        scalingTD_m = AbstractTimeDependence::getTimeDependence(scalingName_m);
    }
}

void MultipoleT::rebuildTanhTable() {
    // The recursions ask for derivatives of the fringe up to order 2 * maxFOrder + 1.
    const unsigned int numDerivatives = 2 * maxFOrder_m + 1;
    tanhTable_m = Kokkos::View<double**>("MultipoleT::tanh", numDerivatives + 1, numDerivatives + 2);
    tanhTableHost_m = Kokkos::create_mirror_view(tanhTable_m);
    MultipoleTFieldModel::tanhDerivativeTable(numDerivatives, tanhTableHost_m);
    Kokkos::deep_copy(tanhTable_m, tanhTableHost_m);
}
