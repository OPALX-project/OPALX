//
// Classes BorisMidpoint and BorisPrestep
//   The Boris schemes: one kick between two half drifts, with the self field solved in the middle
//   of the step or before it.
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
#ifndef OPALX_BORIS_INTEGRATORS_H
#define OPALX_BORIS_INTEGRATORS_H

#include "Steppers/NumericalIntegrator.h"

#include <utility>

/**
 * @brief Drift-kick-drift with the self field solved at @f$R_{n+1/2}@f$; second order. The
 *        default scheme, "BORIS".
 */
class BorisMidpoint final : public NumericalIntegrator {
public:
    explicit BorisMidpoint(Context context) : NumericalIntegrator(std::move(context)) {}

    void advance() override {
        halfDrift();
        solveSelfFields();
        emit();
        addExternalFields();
        kick();
        halfDrift();
    }
};

/**
 * @brief The self field is solved at @f$R_n@f$ and carried by the particles through the first
 *        half drift; first order in the self-field term. Scheme "BORIS-PRESTEP".
 *
 * This is the timing of the OPAL-cycl leapfrog (field solved before the step and held), not of
 * OPAL-T, which solved at the midpoint.
 */
class BorisPrestep final : public NumericalIntegrator {
public:
    explicit BorisPrestep(Context context) : NumericalIntegrator(std::move(context)) {}

    void advance() override {
        solveSelfFields();
        halfDrift();
        emit();
        addExternalFields();
        kick();
        halfDrift();
    }
};

#endif  // OPALX_BORIS_INTEGRATORS_H
