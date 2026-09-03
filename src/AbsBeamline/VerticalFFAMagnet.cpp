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

#include <cmath>

template <class EFM>
VerticalFFAMagnet<EFM>::VerticalFFAMagnet(const std::string& name)
    : ElementBase(name), straightGeometry_m(Geometry::makeStraight(1.)) {}

template <class EFM>
VerticalFFAMagnet<EFM>::VerticalFFAMagnet(const VerticalFFAMagnet& right)
    : ElementBase(right),
      config_m(right.config_m) {
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
    config_m.dfCoefficients_m    = std::vector<std::vector<double> >(config_m.maxOrder_m + 1);
    config_m.dfCoefficients_m[0] = std::vector<double>(1, 1.);
    if (config_m.maxOrder_m > 0) {
        config_m.dfCoefficients_m[1] = std::vector<double>();
    }
    // n indexes like the polynomial order of the midplane expansion
    // e.g. Bz = exp(mz) f_n y^n
    // where y is distance from the midplane and z is height
    for (size_t n = 2; n < config_m.dfCoefficients_m.size(); n += 2) {
        const std::vector<double>& oldCoefficients = config_m.dfCoefficients_m[n - 2];
        std::vector<double> coefficients(oldCoefficients.size() + 2, 0);
        // j indexes the derivative of f_0
        for (size_t j = 0; j < oldCoefficients.size(); ++j) {
            coefficients[j] += -1. / (n) / (n - 1) * config_m.k_m * config_m.k_m * oldCoefficients[j];
            coefficients[j + 2] += -1. / (n) / (n - 1) * oldCoefficients[j];
        }
        config_m.dfCoefficients_m[n] = coefficients;
    }
}

template <class EFM>
void VerticalFFAMagnet<EFM>::setEndField(EFM endField) {
    config_m.endField_m  = endField;
    config_m.endField_m.setMaximumDerivative(config_m.maxOrder_m);
}

template <class EFM>
void VerticalFFAMagnet<EFM>::setMaxOrder(size_t maxOrder) {
    config_m.endField_m.setMaximumDerivative(maxOrder);
    config_m.maxOrder_m = maxOrder;
}
