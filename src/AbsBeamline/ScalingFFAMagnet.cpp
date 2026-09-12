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
#include "AbsBeamline/EndFieldModel/EndFieldModelManager.h"
#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Utilities/GeneralOpalException.h"

extern Inform* gmsg;

ScalingFFAMagnet::ScalingFFAMagnet(const std::string& name)
    : ElementBase(name), planarArcGeometry_m(Geometry::makeSBend(1., 1.)) {}

ScalingFFAMagnet::ScalingFFAMagnet(const ScalingFFAMagnet& right)
    : ElementBase(right),
      planarArcGeometry_m(right.planarArcGeometry_m),
      config_m(right.config_m),
      endFieldName_m(right.endFieldName_m),
      dfCoefficients_m(right.dfCoefficients_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

ScalingFFAMagnet* ScalingFFAMagnet::clone() const {
    auto* magnet = new ScalingFFAMagnet(*this);
    magnet->efm_m = efm_m;
    magnet->initialise();
    return magnet;
}

void ScalingFFAMagnet::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    // Kernel launch over all particles
    getFieldValue(config_m, efm_m, pc);
}

void ScalingFFAMagnet::getCylindricalCoordinates(const Vector_t<double, 3>& R,
                                                 Vector_t<double, 5>& Rcyl) {
    getCylindricalCoordinates(config_m, R, Rcyl);
}

void ScalingFFAMagnet::getFieldValue(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    Vector_t<double, 5> Rffa;
    Vector_t<double, 3> Bcyl;
    getCylindricalCoordinates(config_m, R, Rffa);
    Vector_t<double, 3> Rcyl = {Rffa[0], Rffa[1], Rffa[2]};
    getFieldValueCylindrical(Rcyl, Bcyl);
    rotateBfield(config_m, Rffa, Bcyl, B);
    // std::cerr << "ScalingFFAManget::getFieldValue Rcyl " << Rcyl << std::endl;
}

void ScalingFFAMagnet::getFieldValueCylindrical(const Vector_t<double, 3>& Rcyl, Vector_t<double, 3>& Bcyl) const {
    Kokkos::View<double*, Kokkos::HostSpace> derivatives(
            "single_derivatives", config_m.maxOrder_m + 1);
    Vector_t<double, 5> Rffa;
    Rffa[0] = Rcyl[0];
    Rffa[1] = Rcyl[1];
    Rffa[2] = Rcyl[2];
    Rffa[3] = std::abs(Rcyl[0]/config_m.r0_m); // rnorm
    Rffa[4] = Rcyl[2]-config_m.tanDelta_m * std::log(Rffa[3])-config_m.phiStart_m; // phispiral
    for (size_t i = 0; i <= config_m.maxOrder_m; ++i)
        derivatives(i) = efm_m->function(Rffa[4], i);
    getFieldValueCylindricalImpl(config_m, derivatives, Rffa, Bcyl);
}

void ScalingFFAMagnet::initialise() {
    calculateDfCoefficients();
    if (efm_m) {
        efm_m->setMaximumDerivative(config_m.maxOrder_m);
    }
}

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
    getFieldValue(R, B);
}

void ScalingFFAMagnet::calculateDfCoefficients() {
    dfCoefficients_m    = std::vector<std::vector<double> >(config_m.maxOrder_m + 1);
    dfCoefficients_m[0] = std::vector<double>(1, 1.);  // f_0 = 1.*0th derivative
    for (size_t n = 0; n < config_m.maxOrder_m; n += 2) {       // n indexes the power in z
        dfCoefficients_m[n + 1] = std::vector<double>(dfCoefficients_m[n].size() + 1, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 1][i + 1] = dfCoefficients_m[n][i] / (n + 1);
        }
        if (n + 1 == config_m.maxOrder_m) {
            break;
        }
        dfCoefficients_m[n + 2] = std::vector<double>(dfCoefficients_m[n].size() + 2, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 2][i] =
                    -(config_m.k_m - n) * (config_m.k_m - n) / (n + 1) * dfCoefficients_m[n][i] / (n + 2);
        }
        for (size_t i = 0; i < dfCoefficients_m[n + 1].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 2][i] +=
                    2 * (config_m.k_m - n) * config_m.tanDelta_m * dfCoefficients_m[n + 1][i] / (n + 2);
            dfCoefficients_m[n + 2][i + 1] -=
                    (1 + config_m.tanDelta_m * config_m.tanDelta_m) * dfCoefficients_m[n + 1][i] / (n + 2);
        }
    }
    for (size_t i = 0; i < ScalingFFAMagnetConfig::CoefficientCount; ++i) {
        config_m.dfCoefficients_m[i] = 0.;
    }
    for (size_t n = 0; n < dfCoefficients_m.size(); ++n) {
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            config_m.dfCoefficients_m[n * (ScalingFFAMagnetConfig::MaxOrder + 1) + i] =
                    dfCoefficients_m[n][i];
        }
    }
}

void ScalingFFAMagnet::setMaxOrder(size_t maxOrder) {
    if (maxOrder > ScalingFFAMagnetConfig::MaxOrder) {
        throw GeneralOpalException(
                "ScalingFFAMagnet::setMaxOrder",
                "GPU-compatible field expansions are limited to order 20");
    }
    config_m.maxOrder_m = maxOrder;
}

// Note this is tested in OpalScalingFFAMagnetTest.*
void ScalingFFAMagnet::setupEndField() const {
    if (endFieldName_m == "") {  // no end field is defined
        return;
    }
    auto efmMan = endfieldmodel::EndFieldModelManager::getEFMManager();
    std::shared_ptr<endfieldmodel::EndFieldModel> efm =
                             efmMan->getEndFieldModel(endFieldName_m);
    efm->rescale(1.0 / getR0());
    double defaultExtent = efm->getEndLength()*4. + efm->getCentreLength();
    if (config_m.phiStart_m < 0.0) {
        config_m.phiStart_m  = defaultExtent / 2.0;
    } else {
        config_m.phiStart_m  = getPhiStart() + efm->getCentreLength() * 0.5;
    }
    if (config_m.phiEnd_m < 0.0) {
        config_m.phiEnd_m = defaultExtent;
    }
    if (config_m.azimuthalExtent_m < 0.0) {
        config_m.azimuthalExtent_m  = efm->getEndLength() * 5. + efm->getCentreLength() / 2.0;
    }
    planarArcGeometry_m.setElementLength(config_m.r0_m * config_m.phiEnd_m);  // length = phi r
    planarArcGeometry_m.setCurvature(1. / config_m.r0_m);
    efm_m = efm;
}
