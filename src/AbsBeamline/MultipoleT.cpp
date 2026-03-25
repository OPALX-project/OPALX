/*
 *  Copyright (c) 2017, Titus Dascalu
 *  Copyright (c) 2018, Martin Duy Tat
 *  Copyright (c) 2019-2023, Chris Rogers
 *  Copyright (c) 2025, Jon Thompson
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

Vector_t<double, 3> MultipoleT::rotateFrame(const Vector_t<double, 3>& R) const {
    Vector_t<double, 3> R_prime(3), R_pprime(3);
    /** Apply two 2D rotation matrices to coordinate vector
     * Rotate around central axis => skew fields
     * Rotate azimuthally => entrance angle
     */
    // 1st rotation
    R_prime[0] = R[0] * std::cos(config_m.rotation_m) + R[1] * std::sin(config_m.rotation_m);
    R_prime[1] = -R[0] * std::sin(config_m.rotation_m) + R[1] * std::cos(config_m.rotation_m);
    R_prime[2] = R[2];
    // 2nd rotation
    R_pprime[0] = R_prime[2] * std::sin(config_m.entranceAngle_m)
                  + R_prime[0] * std::cos(config_m.entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * std::cos(config_m.entranceAngle_m)
                  - R_prime[0] * std::sin(config_m.entranceAngle_m);
    return R_pprime;
}

double MultipoleT::getScaling(const double t) const {
    double scaling = 1.0;
    if (scalingTD_m) {
        scaling = scalingTD_m->getValue(t);
    }
    return scaling;
}

bool MultipoleT::apply() {
    const auto pc = RefPartBunch_m->getParticleContainer();
    apply(
        pc->R.getView(), pc->E.getView(), pc->B.getView(), RefPartBunch_m->getT(),
        pc->getLocalNum());
    return false;
}

void MultipoleT::apply(
    const Kokkos::View<Vector_t<double, 3>*>& R, Kokkos::View<Vector_t<double, 3>*>& E,
    Kokkos::View<Vector_t<double, 3>*>& B, const double t, const size_t count) const {
    implementation_->getField(R, E, B, getScaling(t), count);
}

bool MultipoleT::apply(
    const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& t,
    Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    return implementation_->getField(R, E, B, getScaling(t));
}

bool MultipoleT::apply(
    const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    const auto pc = RefPartBunch_m->getParticleContainer();
    return implementation_->getField(pc->R.getView()(i), E, B, getScaling(t));
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
    implementation_->setMaxOrder(config_m.maxFOrder_m, config_m.maxXOrder_m);
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

Vector_t<double, 3> MultipoleT::localCartesianToOpalCartesian(const Vector_t<double, 3>& r) const {
    return implementation_->localCartesianToOpalCartesian(r);
}

double MultipoleT::localCartesianRotation() const {
    return implementation_->localCartesianRotation();
}
