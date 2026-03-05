//
// Class BorisPusher
//   Boris-Buneman time integrator.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_PartPusher_H
#define CLASSIC_PartPusher_H

#include "Algorithms/PartData.h"

#include "Physics/Physics.h"

#include "Expression/IpplExpressions.h"

/*


*/

extern Inform* gmsg;

class BorisPusher {
public:

    KOKKOS_INLINE_FUNCTION BorisPusher(const PartData& ref);
    KOKKOS_INLINE_FUNCTION BorisPusher();

    KOKKOS_INLINE_FUNCTION void initialise(const PartData* ref);

    KOKKOS_INLINE_FUNCTION void kick(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
        const Vector_t<double, 3>& Bf, const double& dt) const;

    KOKKOS_INLINE_FUNCTION void kick(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
        const Vector_t<double, 3>& Bf, const double& dt, const double& mass,
        const double& charge) const;

    KOKKOS_INLINE_FUNCTION void push(Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& dt) const;

private:
    const PartData* itsReference;
};

KOKKOS_INLINE_FUNCTION BorisPusher::BorisPusher(const PartData& ref) : itsReference(&ref) {
}

KOKKOS_INLINE_FUNCTION BorisPusher::BorisPusher() : itsReference(nullptr) {
}

KOKKOS_INLINE_FUNCTION void BorisPusher::initialise(const PartData* ref) {
    itsReference = ref;
}

KOKKOS_INLINE_FUNCTION void BorisPusher::kick(
    const Vector_t<double, 3>& R, Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
    const Vector_t<double, 3>& Bf, const double& dt) const {
    kick(R, P, Ef, Bf, dt, itsReference->getM(), itsReference->getQ());
}

KOKKOS_INLINE_FUNCTION void BorisPusher::kick(
    const Vector_t<double, 3>& /*R*/, Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
    const Vector_t<double, 3>& Bf, const double& dt, const double& mass,
    const double& charge) const {
    // Implementation follows chapter 4-4, p. 61 - 63 from
    // Birdsall, C. K. and Langdon, A. B. (1985). Plasma physics
    // via computer simulation.
    //
    // Up to finite precision effects, the new implementation is equivalent to the
    // old one, but uses less floating point operations.
    //
    // Relativistic variant implemented below is described in
    // chapter 15-4, p. 356 - 357.
    // However, since other units are used here, small
    // modifications are required. The relativistic variant can be derived
    // from the nonrelativistic one by replacing
    //     mass
    // by
    //     gamma * rest mass
    // and transforming the units.
    //
    // Parameters:
    //     R = x / (c * dt): Scaled position x, not used in here
    //     P = v / c * gamma: Scaled velocity v
    //     Ef: Electric field
    //     Bf: Magnetic field
    //     dt: Timestep
    //     mass = rest energy = rest mass * c * c
    //     charge

    // Half step E
    P += 0.5 * dt * charge * Physics::c / mass * Ef;

    // Full step B

    /*
        LF
    double const gamma = sqrt(1.0 + dot(P, P).apply());
    Vector_t<double, 3> const t = dt * charge * c * c / (gamma * mass) * Bf;
    P +=  cross(P, t);
    */

    double gamma                = Kokkos::sqrt(1.0 + dot(P, P));
    Vector_t<double, 3> const t = 0.5 * dt * charge * Physics::c * Physics::c / (gamma * mass) * Bf;
    Vector_t<double, 3> const w = P + cross(P, t);
    Vector_t<double, 3> const s = 2.0 / (1.0 + dot(t, t)) * t;
    P += cross(w, s);

    /* a poor Leap-Frog
        P += 1.0 * dt * charge * c * c / (gamma * mass) * cross(P,Bf);
    */

    // Half step E
    P += 0.5 * dt * charge * Physics::c / mass * Ef;
}

KOKKOS_INLINE_FUNCTION void BorisPusher::push(
    Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& /* dt */) const {
    /** \f[ \vec{x}_{n+1/2} = \vec{x}_{n} + \frac{1}{2}\vec{v}_{n-1/2}\quad (= \vec{x}_{n} +
     * \frac{\Delta t}{2} \frac{\vec{\beta}_{n-1/2}\gamma_{n-1/2}}{\gamma_{n-1/2}}) \f]
     *
     * \code
     * R[i] += 0.5 * P[i] * recpgamma;
     * \endcode
     */
    R += 0.5 * P / Kokkos::sqrt(1.0 + dot(P));
}

#endif
