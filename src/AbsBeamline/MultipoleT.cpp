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
#include "AbsBeamline/MultipoleTFunctions/tanhDeriv.h"
#include "MultipoleTStraight.h"
#include "MultipoleTCurvedConstRadius.h"
#include "MultipoleTCurvedVarRadius.h"

using namespace endfieldmodel;

MultipoleT::MultipoleT(const std::string& name)
    : Component(name) {
    chooseImplementation();
}

MultipoleT::MultipoleT(const MultipoleT& right)
    : Component(right),
    config_m(right.config_m),
    scalingName_m(right.scalingName_m),
    scalingTD_m(right.scalingTD_m) {
    RefPartBunch_m = right.RefPartBunch_m;
    chooseImplementation();
}

ElementBase* MultipoleT::clone() const {
    return new MultipoleT(*this);
}

void MultipoleT::accept(BeamlineVisitor& visitor) const {
    initialiseTimeDepencencies();
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
    R_pprime[0] = R_prime[2] * std::sin(config_m.entranceAngle_m) +
            R_prime[0] * std::cos(config_m.entranceAngle_m);
    R_pprime[1] = R_prime[1];
    R_pprime[2] = R_prime[2] * std::cos(config_m.entranceAngle_m) -
            R_prime[0] * std::sin(config_m.entranceAngle_m);
    return R_pprime;
}

bool MultipoleT::insideAperture(const Vector_t<double, 3>& R) const {
    return std::abs(R[1]) <= config_m.verticalAperture_m / 2.0 &&
            std::abs(R[0]) <= config_m.horizontalAperture_m / 2.0;
}

bool MultipoleT::insideBoundingBox(const Vector_t<double, 3>& R) const {
    return config_m.boundingBoxLength_m == 0.0 || fabs(R[2]) <= config_m.boundingBoxLength_m / 2.0;
}

Vector_t<double, 3> MultipoleT::toMagnetCoords(const Vector_t<double, 3>& R) {
    /** Rotate coordinates around the central axis of the magnet */
    /* TODO:  I'm not sure this is the correct thing to do */
    Vector_t<double, 3> result = rotateFrame(R);
    /** Go to local Frenet-Serret coordinates */
    result[2] *= -1; // OPAL uses a different sign convention...
    implementation_->transformCoords(result);
    return result;
}

Vector_t<double, 3> MultipoleT::getField(const Vector_t<double, 3>& magnetCoords) {
    Vector_t<double, 3> result;
    /** Calculate B-field in the local Frenet-Serret frame */
    result[0] = implementation_->getBx(magnetCoords);
    result[1] = implementation_->getBz(magnetCoords);
    result[2] = implementation_->getBs(magnetCoords);
    /** Transform B-field from local to lab coordinates */
    implementation_->transformBField(result, magnetCoords);
    result[2] *= -1; // OPAL uses a different sign convention...
    return result;
}

double MultipoleT::getScaling(const double t) const {
    double scaling = 1.0;
    if(scalingTD_m) {
        scaling = scalingTD_m->getValue(t);
    }
    return scaling;
}

bool MultipoleT::apply() {
    const auto pc = RefPartBunch_m->getParticleContainer();
    apply(pc->R.getView(), pc->E.getView(), pc->B.getView(), RefPartBunch_m->getT());
    return false;
}

void MultipoleT::apply(const Kokkos::View<Vector_t<double, 3>*>& R,
        Kokkos::View<Vector_t<double, 3>*>& E, Kokkos::View<Vector_t<double, 3>*>& B,
        const double t) const {
    implementation_->getField(R, E, B, getScaling(t));
}

bool MultipoleT::apply(const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/,
        const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B) {
    implementation_->getField(R, E, B, getScaling(t));
    return false;
}

bool MultipoleT::apply(const size_t& i, const double& t, Vector_t<double, 3>& E,
        Vector_t<double, 3>& B) {
    const auto pc = RefPartBunch_m->getParticleContainer();
    implementation_->getField(pc->R.getView()(i), E, B, getScaling(t));
    return false;
}

void MultipoleT::setFringeField(const double& s0, const double& lambda_l, const double& lambda_r) {
    config_m.fringeS0_m = s0;
    config_m.fringeLambdaLeft_m = lambda_l;
    config_m.fringeLambdaRight_m = lambda_r;
    implementation_->initialise();
}

std::tuple<double, double, double> MultipoleT::getFringeField() const {
    return {config_m.fringeS0_m, config_m.fringeLambdaLeft_m, config_m.fringeLambdaRight_m};
}

double MultipoleT::getFringeDeriv(const std::size_t& n, const double& s) {
    double result{};
#if 0
    if(n <= 10) {
        result = (fringeField_l.getTanh(s, static_cast<int>(n)) -
            fringeField_r.getNegTanh(s, static_cast<int>(n))) / 2;
    } else {
        result = tanhderiv::integrate(
                s, fringeField_l.getX0(), fringeField_l.getLambda(),
                fringeField_r.getLambda(), static_cast<int>(n));
    }
#endif
    return result;
}

double MultipoleT::getTransDeriv(const std::size_t& n, const double& x) const {
    double func = 0.0;
    std::size_t transMaxOrder = getTransMaxOrder();
    if(n > transMaxOrder) {
        return func;
    }
    auto temp = getTransProfile();
    for(std::size_t i = 1; i <= n; i++) {
        for(std::size_t j = 0; j <= transMaxOrder; j++) {
            if(j <= config_m.transverseProfileMaxOrder_m - i) {
                temp[j] = temp[j + 1] * static_cast<double>(j + 1);
            } else {
                temp[j] = 0.0;
            }
        }
    }
    std::size_t k = transMaxOrder - n + 1;
    while(k != 0) {
        k--;
        func = func * x + temp[k];
    }
    return func;
}

double MultipoleT::getFnDerivX(const std::size_t& n, const double& x, const double& s) {
    if(n == 0) {
        return getTransDeriv(1, x) * getFringeDeriv(0, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * implementation_->getFn(n, x - 2. * stepSize, s);
    deriv += -8. * implementation_->getFn(n, x - stepSize, s);
    deriv += 8. * implementation_->getFn(n, x + stepSize, s);
    deriv += -1. * implementation_->getFn(n, x + 2. * stepSize, s);
    deriv /= 12 * stepSize;
    return deriv;
}

double MultipoleT::getFnDerivS(const std::size_t& n, const double& x, const double& s) {
    if(n == 0) {
        return getTransDeriv(0, x) * getFringeDeriv(1, s);
    }
    double deriv = 0.0;
    double stepSize = 1e-3;
    deriv += 1. * implementation_->getFn(n, x, s - 2. * stepSize);
    deriv += -8. * implementation_->getFn(n, x, s - stepSize);
    deriv += 8. * implementation_->getFn(n, x, s + stepSize);
    deriv += -1. * implementation_->getFn(n, x, s + 2. * stepSize);
    deriv /= 12 * stepSize;
    return deriv;
}

void MultipoleT::finalise() {
    RefPartBunch_m = nullptr;
}

void MultipoleT::setElementLength(double length) {
    // Base class first
    Component::setElementLength(length);
    // Then me
    config_m.length_m = length;
    implementation_->initialise();
}

void MultipoleT::setBendAngle(double angle, bool variableRadius) {
    // Record information
    config_m.bendAngle_m = angle;
    config_m.variableRadius_m = variableRadius;
    chooseImplementation();
}

void MultipoleT::chooseImplementation() {
    if(config_m.bendAngle_m == 0.0) {
        implementation_ = std::make_unique<MultipoleTStraight>(this);
    } else if(config_m.variableRadius_m) {
        implementation_ = std::make_unique<MultipoleTCurvedVarRadius>(this);
    } else {
        implementation_ = std::make_unique<MultipoleTCurvedConstRadius>(this);
    }
    implementation_->initialise();
}

void MultipoleT::setAperture(const double& vertAp, const double& horizAp) {
    config_m.verticalAperture_m = vertAp;
    config_m.horizontalAperture_m = horizAp;
}

void MultipoleT::setBoundingBoxLength(double boundingBoxLength) {
    config_m.boundingBoxLength_m = boundingBoxLength;
}

void MultipoleT::setTransProfile(const std::vector<double>& profile) {
    config_m.transverseProfileMaxOrder_m = 1;
    for(unsigned int i = 0; i < MultipoleTConfig::NumPoles; ++i) {
        if(i < profile.size() && profile[i] != 0.0) {
            config_m.transverseProfile_m[i] = profile[i];
            config_m.transverseProfileMaxOrder_m =
                    std::max(config_m.transverseProfileMaxOrder_m, i);
        } else {
            config_m.transverseProfile_m[i] = 0.0;
        }
    }
}

void MultipoleT::setMaxOrder(size_t orderZ, size_t orderX) {
    config_m.maxFOrder_m = orderZ;
    config_m.maxXOrder_m = orderX;
    implementation_->setMaxOrder(config_m.maxFOrder_m, config_m.maxXOrder_m);
    implementation_->initialise();
}

void MultipoleT::setRotation(double rot) {
    config_m.rotation_m = rot;
}

void MultipoleT::setEntranceAngle(double entranceAngle) {
    config_m.entranceAngle_m = entranceAngle;
}

void MultipoleT::setEntryOffset(double offset) {
    config_m.entryOffset_m = offset;
}

bool MultipoleT::bends() const {
    return config_m.transverseProfile_m[0] != 0 || config_m.bendAngle_m != 0.0;
}

void MultipoleT::initialise(PartBunch_t* bunch,
        double& /*startField*/, double& /*endField*/) {
    RefPartBunch_m = bunch;
    implementation_->initialise();
}

void MultipoleT::setScalingName(const std::string& name) {
    // Element names are stored in upper case
    scalingName_m = name;
    std::ranges::transform(scalingName_m, scalingName_m.begin(),
            [](const unsigned char c) {
                return static_cast<char>(std::toupper(c));
            });
}

void MultipoleT::initialiseTimeDepencencies() const {
    scalingTD_m.reset();
    if(!scalingName_m.empty()) {
        scalingTD_m = AbstractTimeDependence::getTimeDependence(scalingName_m);
    }
}

BGeometryBase& MultipoleT::getGeometry() {
    return *implementation_->getGeometry();
}

const BGeometryBase& MultipoleT::getGeometry() const {
    return *implementation_->getGeometry();
}

Vector_t<double, 3> MultipoleT::localCartesianToOpalCartesian(const Vector_t<double, 3>& r) {
    return implementation_->localCartesianToOpalCartesian(r);
}

double MultipoleT::localCartesianRotation() {
    return implementation_->localCartesianRotation();
}
