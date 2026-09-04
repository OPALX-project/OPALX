/*
 *  Copyright (c) 2017, Chris Rogers
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

#include <cmath>

#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"
#include "AbsBeamline/ScalingFFAMagnet.h"

ScalingFFAMagnet::ScalingFFAMagnet(const std::string& name)
    : ElementBase(name), planarArcGeometry_m(Geometry::makeSBend(1., 1.)) {}

ScalingFFAMagnet::ScalingFFAMagnet(const ScalingFFAMagnet& right)
    : ElementBase(right),
      planarArcGeometry_m(right.planarArcGeometry_m),
      config_m(right.config_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

ScalingFFAMagnet::~ScalingFFAMagnet() {}

ScalingFFAMagnet* ScalingFFAMagnet::clone() const {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

void ScalingFFAMagnet::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    // Kernel launch over all particles
    const ScalingFFAMagnetConfig config = getConfig();
    getFieldValue<>(config, *tanh, pc);
}

bool ScalingFFAMagnet::getFieldValue(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    return getFieldValue<>(config_m, *tanh, R, B);
}

void ScalingFFAMagnet::initialise() { calculateDfCoefficients(); }

void ScalingFFAMagnet::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    initialise();
}

void ScalingFFAMagnet::finalise() { RefPartBunch_m = nullptr; }

Geometry& ScalingFFAMagnet::getGeometry() { return planarArcGeometry_m; }

const Geometry& ScalingFFAMagnet::getGeometry() const { return planarArcGeometry_m; }

void ScalingFFAMagnet::accept(BeamlineVisitor& visitor) const {
    visitor.visitScalingFFAMagnet(*this);
    setupEndField();
}

void ScalingFFAMagnet::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) {
    getFieldValue<>(config_m, *tanh, R, B);
}

void ScalingFFAMagnet::calculateDfCoefficients() {
    config_m.dfCoefficients_m    = std::vector<std::vector<double> >(config_m.maxOrder_m + 1);
    config_m.dfCoefficients_m[0] = std::vector<double>(1, 1.);  // f_0 = 1.*0th derivative
    for (size_t n = 0; n < config_m.maxOrder_m; n += 2) {       // n indexes the power in z
        config_m.dfCoefficients_m[n + 1] = std::vector<double>(config_m.dfCoefficients_m[n].size() + 1, 0);
        for (size_t i = 0; i < config_m.dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            config_m.dfCoefficients_m[n + 1][i + 1] = config_m.dfCoefficients_m[n][i] / (n + 1);
        }
        if (n + 1 == config_m.maxOrder_m) {
            break;
        }
        config_m.dfCoefficients_m[n + 2] = std::vector<double>(config_m.dfCoefficients_m[n].size() + 2, 0);
        for (size_t i = 0; i < config_m.dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            config_m.dfCoefficients_m[n + 2][i] =
                    -(config_m.k_m - n) * (config_m.k_m - n) / (n + 1) * config_m.dfCoefficients_m[n][i] / (n + 2);
        }
        for (size_t i = 0; i < config_m.dfCoefficients_m[n + 1].size(); ++i) {  // i indexes the derivative
            config_m.dfCoefficients_m[n + 2][i] +=
                    2 * (config_m.k_m - n) * config_m.tanDelta_m * config_m.dfCoefficients_m[n + 1][i] / (n + 2);
            config_m.dfCoefficients_m[n + 2][i + 1] -=
                    (1 + config_m.tanDelta_m * config_m.tanDelta_m) * config_m.dfCoefficients_m[n + 1][i] / (n + 2);
        }
    }
}

extern Inform* gmsg;

// Note this is tested in OpalScalingFFAMagnetTest.*
void ScalingFFAMagnet::setupEndField() const {
    if (config_m.endFieldName_m == "") {  // no end field is defined
        return;
    }


    tanh->rescale(1.0 / getR0());
    double defaultExtent = (tanh->getEndLength() * 4. + tanh->getCentreLength());
    if (config_m.phiStart_m < 0.0) {
        config_m.phiStart_m  = defaultExtent / 2.0;
    } else {
        config_m.phiStart_m  = getPhiStart() + tanh->getCentreLength() * 0.5;
    }
    if (config_m.phiEnd_m < 0.0) {
        config_m.phiEnd_m = defaultExtent;
    }
    if (config_m.azimuthalExtent_m < 0.0) {
        config_m.azimuthalExtent_m  = tanh->getEndLength() * 5. + tanh->getCentreLength() / 2.0;
    }
    planarArcGeometry_m.setElementLength(config_m.r0_m * config_m.phiEnd_m);  // length = phi r
    planarArcGeometry_m.setCurvature(1. / config_m.r0_m);
    std::cerr << "ScalingFFAMagnet::setupEndField " << std::endl;
    std::cerr << "    name: " << getName() << std::endl;
    std::cerr << "    R0: " << config_m.r0_m << " phiend " << config_m.phiEnd_m << std::endl;
    std::cerr << "    geom: " << planarArcGeometry_m.getChordLength() << " " << planarArcGeometry_m.getBendAngle() << std::endl;
}

template <>
std::shared_ptr<endfieldmodel::Tanh> ScalingFFAMagnet::getEndField<endfieldmodel::Tanh>() const {
    return tanh;
}

template <>
void ScalingFFAMagnet::setEndField(std::shared_ptr<endfieldmodel::Tanh> endField) {
    tanh = endField;
}


