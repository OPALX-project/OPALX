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

using namespace endfieldmodel;

MultipoleT::MultipoleT(const std::string& name) : Component(name) { chooseImplementation(); }

MultipoleT::MultipoleT(const MultipoleT& right)
    : Component(right),
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

double MultipoleT::getScaling(const double t) const {
    double scaling = 1.0;
    if (scalingTD_m) {
        scaling = scalingTD_m->getValue(t);
    }
    return scaling;
}

bool MultipoleT::apply() {
    validateConfiguration();
    const auto pc = RefPartBunch_m->getParticleContainer();
    implementation_->getField(
            pc->R.getView(), pc->E.getView(), pc->B.getView(), getScaling(RefPartBunch_m->getT()),
            pc->getLocalNum());
    return false;
}

bool MultipoleT::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
        Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    validateConfiguration();
    return implementation_->getField(R, E, B, getScaling(t));
}

bool MultipoleT::apply(
        const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    validateConfiguration();
    const auto pc = RefPartBunch_m->getParticleContainer();
    Vector_t<double, 3> R{};
    Kokkos::deep_copy(
            Kokkos::View<Vector_t<double, 3>, Kokkos::HostSpace>(&R),
            Kokkos::subview(pc->R.getView(), i));
    return implementation_->getField(R, E, B, getScaling(t));
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
    // Base class first
    Component::setElementLength(length);
    // Then me
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

bool MultipoleT::bends() const {
    return config_m.transverseProfile_m[0] != 0 || config_m.bendAngle_m != 0.0;
}

void MultipoleT::initialise(PartBunch_t* bunch, double& /*startField*/, double& /*endField*/) {
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

BGeometryBase& MultipoleT::getGeometry() { return *implementation_->getGeometry(); }

const BGeometryBase& MultipoleT::getGeometry() const { return *implementation_->getGeometry(); }

void MultipoleT::validateConfiguration() const {
    if (2 * config_m.maxFOrder_m + 1 > MultipoleTBase::MaxDerivatives) {
        throw OpalException(
            "MultipoleT::validateConfiguration",
            "Max F order too large for this implementation");
    }
}
