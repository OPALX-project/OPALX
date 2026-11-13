//
// Cubic Spline Interpolation to replace GSL spline
//
// Copyright (c) 2023, Paul Scherrer Institute, Villigen PSI, Switzerland
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

#include "MultipoleT.h"
#include "BeamlineVisitor.h"
#include "MultipoleTCurvedConstRadius.h"
#include "MultipoleTStraight.h"

#include "BeamlineGeometry/Geometry.h"
#include "PartBunch/PartBunch.h"

#include <algorithm>
#include <cmath>

MultipoleT::MultipoleT(const std::string& name) : ElementBase(name) { chooseImplementation(); }

MultipoleT::MultipoleT(const MultipoleT& right)
    : ElementBase(right),
      config_m(right.config_m),
      scalingName_m(right.scalingName_m),
      scalingTD_m(right.scalingTD_m) {
    RefPartBunch_m = right.RefPartBunch_m;
    chooseImplementation();
}

ElementBase* MultipoleT::clone() const { return new MultipoleT(*this); }

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    initialiseTimeDependencies();
    visitor.visitMultipoleT(*this);
}

ElementType MultipoleT::getType() const { return ElementType::MULTIPOLET; }

double MultipoleT::getScaling(const double t) const {
    double scaling = 1.0;
    if (scalingTD_m) {
        scaling = scalingTD_m->getValue(t);
    }
    return scaling;
}

void MultipoleT::apply(const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& time) {
    validateConfiguration();
    implementation_->getField(
            pc->R.getView(), pc->E.getView(), pc->B.getView(), getScaling(time.t),
            pc->getLocalNum(), getRegions());
}

void MultipoleT::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    validateConfiguration();
    if (!getRegions().inFieldRegion(R)) return;
    implementation_->getField(R, E, B, getScaling(t));
}

void MultipoleT::setFringeField(
        const double& s0, const double& lambda_left, const double& lambda_right) {
    config_m.fringeS0_m          = s0;
    config_m.fringeLambdaLeft_m  = lambda_left;
    config_m.fringeLambdaRight_m = lambda_right;
    implementation_->initialise();
}

std::tuple<double, double, double> MultipoleT::getFringeField() const {
    return {config_m.fringeS0_m, config_m.fringeLambdaLeft_m, config_m.fringeLambdaRight_m};
}

void MultipoleT::finalise() { RefPartBunch_m = nullptr; }

void MultipoleT::setElementLength(const double length) {
    // Only the config is authoritative; initialise() pushes it into the geometry.
    config_m.length_m = length;
    implementation_->initialise();
}

void MultipoleT::setBendAngle(const double angle, const bool variableRadius) {
    // Record information
    config_m.bendAngle_m      = angle;
    config_m.variableRadius_m = variableRadius;
    chooseImplementation();
}

void MultipoleT::chooseImplementation() {
    if (config_m.bendAngle_m == 0.0) {
        implementation_ = std::make_unique<MultipoleTStraight>(this);
        // This is where the variable radius code is to be patched in.
        //} else if (config_m.variableRadius_m) {
        //    implementation_ = std::make_unique<MultipoleTCurvedVarRadius>(this);
    } else {
        implementation_ = std::make_unique<MultipoleTCurvedConstRadius>(this);
    }
    implementation_->initialise();
}

void MultipoleT::setAperture(const double& vertAp, const double& horizAp) {
    config_m.verticalAperture_m   = vertAp;
    config_m.horizontalAperture_m = horizAp;
}

void MultipoleT::setBoundingBoxLength(const double boundingBoxLength) {
    config_m.boundingBoxLength_m = boundingBoxLength;
}

void MultipoleT::setTransProfile(const std::vector<double>& profile) {
    config_m.transverseProfileMaxOrder_m = 1;
    for (unsigned int i = 0; i < MultipoleTConfig::NumPoles; ++i) {
        if (i < profile.size() && profile[i] != 0.0) {
            config_m.transverseProfile_m[i] = profile[i];
            config_m.transverseProfileMaxOrder_m =
                    std::max(config_m.transverseProfileMaxOrder_m, i);
        } else {
            config_m.transverseProfile_m[i] = 0.0;
        }
    }
}

void MultipoleT::setMaxOrder(const size_t orderZ, const size_t orderX) {
    config_m.maxFOrder_m = orderZ;
    config_m.maxXOrder_m = orderX;
    implementation_->initialise();
}

void MultipoleT::setRotation(const double rot) { config_m.rotation_m = rot; }

void MultipoleT::setEntranceAngle(const double entranceAngle) {
    config_m.entranceAngle_m = entranceAngle;
}

void MultipoleT::setEntryOffset(const double offset) { config_m.entryOffset_m = offset; }

void MultipoleT::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    implementation_->initialise();
}

void MultipoleT::setScalingName(const std::string& name) {
    // Element names are stored in upper case
    scalingName_m = name;
    std::ranges::transform(scalingName_m, scalingName_m.begin(), [](const unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
}

void MultipoleT::initialiseTimeDependencies() const {
    scalingTD_m.reset();
    if (!scalingName_m.empty()) {
        scalingTD_m = AbstractTimeDependence::getTimeDependence(scalingName_m);
    }
}

Geometry& MultipoleT::getGeometry() { return *implementation_->getGeometry(); }

const Geometry& MultipoleT::getGeometry() const { return *implementation_->getGeometry(); }

Vector_t<double, 3> MultipoleT::bendCoords(const Vector_t<double, 3>& r) const {
    // The stored frame is the design-orbit entrance tangent, so the arc coordinate is
    // measured directly. A straight MultipoleT has zero curvature and r comes back unchanged.
    return GeometryHelper::toBendArcCoords(
            r, getGeometry().getCurvature(), getGeometry().getElementLength());
}

BoundingBox MultipoleT::getLocalSupportBox(const double reach) const {
    BoundingBox box = ElementBase::getLocalSupportBox(reach);

    // Kernel gate in the magnet frame of MultipoleTBase::rotateFrame(). Without BBLENGTH the
    // kernel has no z cut, so the selection extent bounds z as it did for the reference path.
    const double length = getGeometry().getElementLength();
    double zBegin       = 0.0;
    double zEnd         = 0.0;
    getFieldExtent(zBegin, zEnd);
    if (config_m.boundingBoxLength_m != 0.0) {
        zBegin = std::min(zBegin, 0.5 * (length - config_m.boundingBoxLength_m));
        zEnd   = std::max(zEnd, 0.5 * (length + config_m.boundingBoxLength_m));
    }
    const double halfX     = 0.5 * config_m.horizontalAperture_m;
    const double halfY     = 0.5 * config_m.verticalAperture_m;
    const double curvature = getGeometry().getCurvature();
    const BoundingBox kernel =
            curvature != 0.0 ? getArcSupportBox(curvature, length, zBegin, zEnd, halfX, halfY)
                             : BoundingBox::getBoundingBox(
                                       {Vector_t<double, 3>({-halfX, -halfY, zBegin}),
                                        Vector_t<double, 3>({halfX, halfY, zEnd})});

    // Undo the EANGLE turn about y, then the ROTATION about z.
    const double cosE         = std::cos(config_m.entranceAngle_m);
    const double sinE         = std::sin(config_m.entranceAngle_m);
    const double cosR         = std::cos(config_m.rotation_m);
    const double sinR         = std::sin(config_m.rotation_m);
    const auto [lower, upper] = kernel.getCorners();
    for (const double x : {lower(0), upper(0)}) {
        for (const double y : {lower(1), upper(1)}) {
            for (const double z : {lower(2), upper(2)}) {
                const double x1 = cosE * x - sinE * z;
                const double z1 = sinE * x + cosE * z;
                box.enlargeToContainPosition(
                        Vector_t<double, 3>({cosR * x1 - sinR * y, sinR * x1 + cosR * y, z1}));
            }
        }
    }
    return box;
}

std::vector<SupportBox> MultipoleT::getLocalSupportBoxes(const double reach) const {
    const double curvature = getGeometry().getCurvature();
    if (curvature == 0.0 || config_m.rotation_m != 0.0 || config_m.entranceAngle_m != 0.0) {
        return ElementBase::getLocalSupportBoxes(reach);
    }
    const double length = getGeometry().getElementLength();
    double zBegin       = 0.0;
    double zEnd         = 0.0;
    getFieldExtent(zBegin, zEnd);
    zBegin = std::min(zBegin, 0.0);
    zEnd   = std::max(zEnd, length);
    if (config_m.boundingBoxLength_m != 0.0) {
        zBegin = std::min(zBegin, 0.5 * (length - config_m.boundingBoxLength_m));
        zEnd   = std::max(zEnd, 0.5 * (length + config_m.boundingBoxLength_m));
    }
    const auto [halfX, halfY] = getSupportHalfWidths(reach);
    return getArcSupportBoxes(
            curvature, length, zBegin, zEnd, std::max(halfX, 0.5 * config_m.horizontalAperture_m),
            std::max(halfY, 0.5 * config_m.verticalAperture_m));
}

std::optional<ApertureShape> MultipoleT::getFieldRange() const {
    // A turned kernel gate is not an axis-aligned shape in the element frame; the kernel keeps
    // its own gate then.
    if (config_m.rotation_m != 0.0 || config_m.entranceAngle_m != 0.0) return std::nullopt;
    ApertureShape range;
    range.halfX = 0.5 * config_m.horizontalAperture_m;
    range.halfY = 0.5 * config_m.verticalAperture_m;
    return range;
}

bool MultipoleT::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const ElementRegions regions = getRegions();
    if (!regions.inFieldRegion(R)) return regions.inLossRegion(R);

    // The reference particle has to see the same field as every other particle. Without this
    // the reference walks a straight line through a curved magnet and never leaves it.
    apply(R, P, t, E, B);
    return false;
}

void MultipoleT::validateConfiguration() const {
    if (2 * config_m.maxFOrder_m + 1 > MultipoleTBase::MaxDerivatives) {
        throw OpalException(
                "MultipoleT::validateConfiguration",
                "Max F order too large for this implementation");
    }
}
