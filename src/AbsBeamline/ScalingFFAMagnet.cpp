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
#include "Physics/Units.h"
#include "AbsBeamline/ScalingFFAMagnet.h"

template <class EFM>
ScalingFFAMagnet<EFM>::ScalingFFAMagnet(const std::string& name)
    : ElementBase(name), planarArcGeometry_m(Geometry::makeSBend(1., 1.)) {}

template <class EFM>
ScalingFFAMagnet<EFM>::ScalingFFAMagnet(const ScalingFFAMagnet<EFM>& right)
    : ElementBase(right),
      planarArcGeometry_m(right.planarArcGeometry_m),
      maxOrder_m(right.maxOrder_m),
      tanDelta_m(right.tanDelta_m),
      k_m(right.k_m),
      Bz_m(right.Bz_m),
      r0_m(right.r0_m),
      rMin_m(right.rMin_m),
      rMax_m(right.rMax_m),
      phiStart_m(right.phiStart_m),
      phiEnd_m(right.phiEnd_m),
      azimuthalExtent_m(right.azimuthalExtent_m),
      verticalExtent_m(right.verticalExtent_m),
      centre_m(right.centre_m),
      endField_m(right.endField_m),
      endFieldName_m(right.endFieldName_m),
      dfCoefficients_m(right.dfCoefficients_m) {
    RefPartBunch_m = right.RefPartBunch_m;
    Bz_m           = right.Bz_m;
    r0_m           = right.r0_m;
}

template <class EFM>
ScalingFFAMagnet<EFM>::~ScalingFFAMagnet() {}

template <class EFM>
ScalingFFAMagnet<EFM>* ScalingFFAMagnet<EFM>::clone() const {
    ScalingFFAMagnet* magnet = new ScalingFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}


template <class EFM>
void ScalingFFAMagnet<EFM>::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    // Kernel launch over all particles
    const Kokkos::View<Vector_t<double, 3>*> R = pc->R.getView();
    const Kokkos::View<Vector_t<double, 3>*> B = pc->B.getView();
    const size_t count = pc->getLocalNum();
    const ScalingFFAMagnetConfig<EFM> config = getConfig();
        Kokkos::parallel_for(
            "ScalingFFAMagnet<>::getFieldValue()", count, KOKKOS_LAMBDA(const size_t i) {
                getFieldValue(config, R(i), B(i));
            });
}

template <class EFM>
void ScalingFFAMagnet<EFM>::initialise() { calculateDfCoefficients(); }

template <class EFM>
void ScalingFFAMagnet<EFM>::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    initialise();
}

template <class EFM>
void ScalingFFAMagnet<EFM>::finalise() { RefPartBunch_m = nullptr; }

template <class EFM>
Geometry& ScalingFFAMagnet<EFM>::getGeometry() { return planarArcGeometry_m; }

template <class EFM>
const Geometry& ScalingFFAMagnet<EFM>::getGeometry() const { return planarArcGeometry_m; }

template <class EFM>
void ScalingFFAMagnet<EFM>::accept(BeamlineVisitor& visitor) const {
    visitor.visitScalingFFAMagnet(*this);
}

template <class EFM>
void ScalingFFAMagnet<EFM>::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& B) {
    getFieldValue(R, B);
}

template <class EFM>
void ScalingFFAMagnet<EFM>::calculateDfCoefficients() {
    dfCoefficients_m    = std::vector<std::vector<double> >(maxOrder_m + 1);
    dfCoefficients_m[0] = std::vector<double>(1, 1.);  // f_0 = 1.*0th derivative
    for (size_t n = 0; n < maxOrder_m; n += 2) {       // n indexes the power in z
        dfCoefficients_m[n + 1] = std::vector<double>(dfCoefficients_m[n].size() + 1, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 1][i + 1] = dfCoefficients_m[n][i] / (n + 1);
        }
        if (n + 1 == maxOrder_m) {
            break;
        }
        dfCoefficients_m[n + 2] = std::vector<double>(dfCoefficients_m[n].size() + 2, 0);
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 2][i] =
                    -(k_m - n) * (k_m - n) / (n + 1) * dfCoefficients_m[n][i] / (n + 2);
        }
        for (size_t i = 0; i < dfCoefficients_m[n + 1].size(); ++i) {  // i indexes the derivative
            dfCoefficients_m[n + 2][i] +=
                    2 * (k_m - n) * tanDelta_m * dfCoefficients_m[n + 1][i] / (n + 2);
            dfCoefficients_m[n + 2][i + 1] -=
                    (1 + tanDelta_m * tanDelta_m) * dfCoefficients_m[n + 1][i] / (n + 2);
        }
    }
}

template <class EFM>
void ScalingFFAMagnet<EFM>::setEndField(EFM endField) {
    endField_m = endField;
}

extern Inform* gmsg;

// Note this is tested in OpalScalingFFAMagnetTest.*
template <class EFM>
void ScalingFFAMagnet<EFM>::setupEndField() {
    if (endFieldName_m == "") {  // no end field is defined
        return;
    }

    //std::shared_ptr<endfieldmodel::EndFieldModel> efm =
    //        endfieldmodel::EndFieldModel::getEndFieldModel(endFieldName_m);
    EFM newEFM = endField_m;
    newEFM.rescale(Units::m2mm * 1.0 / getR0());
    setEndField(newEFM);

    double defaultExtent = (newEFM.getEndLength() * 4. + newEFM.getCentreLength());
    if (phiStart_m < 0.0) {
        setPhiStart(defaultExtent / 2.0);
    } else {
        setPhiStart(getPhiStart() + newEFM.getCentreLength() * 0.5);
    }
    if (phiEnd_m < 0.0) {
        setPhiEnd(defaultExtent);
    }
    if (azimuthalExtent_m < 0.0) {
        setAzimuthalExtent(newEFM.getEndLength() * 5. + newEFM.getCentreLength() / 2.0);
    }
    planarArcGeometry_m.setElementLength(r0_m * phiEnd_m);  // length = phi r
    planarArcGeometry_m.setCurvature(1. / r0_m);
}
