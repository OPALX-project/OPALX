//
// Function phaseSpaceRate
//   The equations of motion of a charged particle in given electric and magnetic fields.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPALX.
//
// OPALX is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPALX. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_PHASE_SPACE_RATE_H
#define OPALX_PHASE_SPACE_RATE_H

#include "Physics/Physics.h"
#include "VectorMath.h"

/// Time derivative of position R (m/s) and momentum P = beta gamma (1/s).
struct PhaseSpaceRate {
    Vector_t<double, 3> dR;
    Vector_t<double, 3> dP;
};

/**
 * @brief @f$\dot R = c P / \gamma@f$ and @f$\dot P = (q c / M)\,(E + \dot R \times B)@f$ in
 *        the fields E (V/m) and B (T), for charge q (e) and rest energy M (eV), the units of
 *        BorisPusher::kick().
 */
KOKKOS_INLINE_FUNCTION PhaseSpaceRate phaseSpaceRate(
        const Vector_t<double, 3>& P, const Vector_t<double, 3>& E, const Vector_t<double, 3>& B,
        double charge, double mass) {
    const double gamma              = Kokkos::sqrt(1.0 + dot(P, P));
    const Vector_t<double, 3> dR    = (Physics::c / gamma) * P;
    const Vector_t<double, 3> force = E + cross(dR, B);
    return {dR, (charge * Physics::c / mass) * force};
}

#endif  // OPALX_PHASE_SPACE_RATE_H
