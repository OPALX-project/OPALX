//
// Source file for VerticalFFAMagnet Component
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#include "AbsBeamline/VerticalFFAMagnet.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/GeneralOpalException.h"

#include <cmath>


VerticalFFAMagnet::VerticalFFAMagnet(const std::string& name)
    : ElementBase(name), straightGeometry_m(Geometry::makeStraight(1.)) {}

VerticalFFAMagnet::VerticalFFAMagnet(const VerticalFFAMagnet& right)
    : ElementBase(right),
      config_m(right.config_m),
      endField_m(right.endField_m),
      dfCoefficients_m(right.dfCoefficients_m) {
    RefPartBunch_m = right.RefPartBunch_m;
}

VerticalFFAMagnet::~VerticalFFAMagnet() {}

ElementBase* VerticalFFAMagnet::clone() const {
    VerticalFFAMagnet* magnet = new VerticalFFAMagnet(*this);
    magnet->initialise();
    return magnet;
}

void VerticalFFAMagnet::initialise() {
    calculateDfCoefficients();
    endField_m->setMaximumDerivative(config_m.maxOrder_m + 1);
    // config_m.endField_m = endField_m.getDeviceData();
    straightGeometry_m.setElementLength(config_m.bbLength_m);  // length = phi r
}

void VerticalFFAMagnet::initialise(PartBunch_t* bunch) {
    RefPartBunch_m = bunch;
    initialise();
}

void VerticalFFAMagnet::finalise() {
    RefPartBunch_m = nullptr;
}

Geometry& VerticalFFAMagnet::getGeometry() {
    return straightGeometry_m;
}


const Geometry& VerticalFFAMagnet::getGeometry() const {
    return straightGeometry_m;
}

void VerticalFFAMagnet::accept(BeamlineVisitor& visitor) const {
    visitor.visitVerticalFFAMagnet(*this);
}

bool VerticalFFAMagnet::getFieldValue(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& B) const {
    return getFieldValue(config_m, R, B);
}

void VerticalFFAMagnet::calculateDfCoefficients() {
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
            coefficients[j] +=
                    -1. / (n) / (n - 1) * config_m.k_m * config_m.k_m * oldCoefficients[j];
            coefficients[j + 2] += -1. / (n) / (n - 1) * oldCoefficients[j];
        }
        dfCoefficients_m[n] = coefficients;
    }
    for (size_t i = 0; i < VerticalFFAMagnetConfig::CoefficientCount; ++i) {
        config_m.dfCoefficients_m[i] = 0.;
    }
    for (size_t n = 0; n < dfCoefficients_m.size(); ++n) {
        for (size_t i = 0; i < dfCoefficients_m[n].size(); ++i) {
            config_m.dfCoefficients_m[n * (VerticalFFAMagnetConfig::MaxOrder + 1) + i] =
                    dfCoefficients_m[n][i];
        }
    }
}

void VerticalFFAMagnet::setEndField(std::shared_ptr<endfieldmodel::EndFieldModel> endField) {
    endField_m = endField;
    endField_m->setMaximumDerivative(config_m.maxOrder_m + 1);
}

void VerticalFFAMagnet::setMaxOrder(size_t maxOrder) {
    if (maxOrder > VerticalFFAMagnetConfig::MaxOrder) {
        throw GeneralOpalException(
                "VerticalFFAMagnet::setMaxOrder",
                "GPU-compatible field expansions are limited to order 20");
    }
    endField_m->setMaximumDerivative(maxOrder + 1);
    // config_m.endField_m = endField_m.getDeviceData();
    config_m.maxOrder_m = maxOrder;
}

