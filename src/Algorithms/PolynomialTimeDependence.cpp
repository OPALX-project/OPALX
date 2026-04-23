//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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

#include "Algorithms/PolynomialTimeDependence.h"
#include <iostream>
#include "Utility/Inform.h"

PolynomialTimeDependence* PolynomialTimeDependence::clone() {
    const auto temp(coeffs);
    auto* d = new PolynomialTimeDependence(temp);
    return d;
}

double PolynomialTimeDependence::getValue(const double time) {
    double x       = 0.;
    double t_power = 1.;
    for (const double coeff : coeffs) {
        x += coeff * t_power;
        t_power *= time;
    }
    return x;
}

double PolynomialTimeDependence::getIntegral(const double time) {
    double result{};
    double t_power{time};
    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * t_power / static_cast<double>(i + 1);
        t_power *= time;
    }
    return result;
}

Inform& PolynomialTimeDependence::print(Inform& os) const {
    const Inform::FmtFlags_t ff = os.flags();
    os << std::scientific;
    for (unsigned int i = 0; i < this->coeffs.size(); i++) {
        os << "c" << i << "= " << this->coeffs[i] << " ";
    }
    os << endl;
    os.flags(ff);
    return os;
}
