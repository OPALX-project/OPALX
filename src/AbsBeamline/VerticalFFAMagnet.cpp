//
// Source file for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "AbsBeamline/BeamlineVisitor.h"
#include "AbsBeamline/VerticalFFAMagnet.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>

template <class EFM>
VerticalFFAMagnet<EFM>::VerticalFFAMagnet(const std::string& name)
    : ElementBase(name), straightGeometry_m(Geometry::makeStraight(1.)) {}

template <class EFM>
VerticalFFAMagnet<EFM>::VerticalFFAMagnet(const VerticalFFAMagnet& right)
    : ElementBase(right),
      config_m(right.config_m),
      endField_m(right.endField_m),
      dfCoefficients_m(right.dfCoefficients_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

template <class EFM>
VerticalFFAMagnet<EFM>::~VerticalFFAMagnet() {}

template <class EFM>
ElementBase* VerticalFFAMagnet<EFM>::clone() const {
    VerticalFFAMagnet* magnet = new VerticalFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

template <class EFM>
void VerticalFFAMagnet<EFM>::initialise() {
    calculateDfCoefficients();
    endField_m.setMaximumDerivative(config_m.maxOrder_m + 1);
    config_m.endField_m = endField_m.getDeviceData();
    straightGeometry_m.setElementLength(config_m.bbLength_m);  // length = phi r
}

template <class EFM>
void VerticalFFAMagnet<EFM>::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    initialise();
}

template <class EFM>
void VerticalFFAMagnet<EFM>::finalise() { RefPartBunch_m = nullptr; }

template <class EFM>
Geometry& VerticalFFAMagnet<EFM>::getGeometry() { return straightGeometry_m; }

template <class EFM>
const Geometry& VerticalFFAMagnet<EFM>::getGeometry() const { return straightGeometry_m; }

template <class EFM>
void VerticalFFAMagnet<EFM>::accept(BeamlineVisitor& visitor) const {
    visitor.visitVerticalFFAMagnet(*this);
}

template <class EFM>
bool VerticalFFAMagnet<EFM>::getFieldValue(const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    return getFieldValue(config_m, R, B);
}

template <class EFM>
void VerticalFFAMagnet<EFM>::calculateDfCoefficients() {
    dfCoefficients_m    = std::vector<std::vector<double> >(config_m.maxOrder_m + 1);
    dfCoefficients_m[0] = std::vector<double>(1, 1.);
    if (config_m.maxOrder_m > 0) {
        dfCoefficients_m[1] = std::vector<double>();
    }
    // n indexes like the polynomial order of the midplane expansion
    // e.g. Bz = exp(mz) f_n y^n
    // where y is distance from the midplane and z is height
    for (size_t n = 2; n < dfCoefficients_m.size(); n += 2) {
        const std::vector<double>& oldCoefficients = dfCoefficients_m[n - 2];
        std::vector<double> coefficients(oldCoefficients.size() + 2, 0);
        // j indexes the derivative of f_0
        for (size_t j = 0; j < oldCoefficients.size(); ++j) {
            coefficients[j] += -1. / (n) / (n - 1) * config_m.k_m * config_m.k_m * oldCoefficients[j];
            coefficients[j + 2] += -1. / (n) / (n - 1) * oldCoefficients[j];
        }
        dfCoefficients_m[n] = coefficients;
    }
    for (size_t i = 0; i < VerticalFFAMagnetConfig<EFM>::CoefficientCount; ++i) {
        config_m.dfCoefficients_m[i] = 0.;
    }
    for (size_t n = 0; n < dfCoefficients_m.size(); ++n) {
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            config_m.dfCoefficients_m[
                    n * (VerticalFFAMagnetConfig<EFM>::MaxOrder + 1) + i] =
                    dfCoefficients_m[n][i];
        }
    }
}

template <class EFM>
void VerticalFFAMagnet<EFM>::setEndField(EFM endField) {
    endField_m = endField;
    endField_m.setMaximumDerivative(config_m.maxOrder_m + 1);
    config_m.endField_m = endField_m.getDeviceData();
}

template <class EFM>
void VerticalFFAMagnet<EFM>::setMaxOrder(size_t maxOrder) {
    if (maxOrder > VerticalFFAMagnetConfig<EFM>::MaxOrder) {
        throw GeneralOpalException(
                "VerticalFFAMagnet::setMaxOrder",
                "GPU-compatible field expansions are limited to order 20");
    }
    endField_m.setMaximumDerivative(maxOrder + 1);
    config_m.endField_m = endField_m.getDeviceData();
    config_m.maxOrder_m = maxOrder;
}

template class VerticalFFAMagnet<endfieldmodel::Tanh>;
