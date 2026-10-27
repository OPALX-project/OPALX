//
// Class SpinTBMTPusher
//   Thomas-BMT spin precession integrator using analytic Rodrigues rotation.
//   Preserves |Pol| exactly per step.
//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPALX_SpinTBMTPusher_H
#define OPALX_SpinTBMTPusher_H

#include "Physics/Physics.h"

#include "Expression/IpplExpressions.h"
#include "VectorMath.h"

class SpinTBMTPusher {
public:
    KOKKOS_INLINE_FUNCTION SpinTBMTPusher() = default;

    /// Per-particle Thomas-BMT update.
    ///
    /// Inputs (all along the same Cartesian axes as Pol, in OPALX conventions; the
    /// tracker passes bunch-frame P, E, B and Pol):
    ///   P       : momentum in beta*gamma (dimensionless)
    ///   Ef, Bf  : electric and magnetic field at the particle position [V/m, T]
    ///   dt      : timestep [s]
    ///   mass    : rest energy [eV]  (m * c^2)
    ///   charge  : particle charge [proton charges]
    ///   anom    : magnetic moment anomaly G = (g - 2) / 2 [dimensionless]
    /// In/out:
    ///   Pol     : polarization 3-vector (rest-frame Pauli expectations along those axes,
    ///             |Pol| in [0, 1]). Updated in place.
    template <typename SpinVec>
    KOKKOS_INLINE_FUNCTION void evolve(
            SpinVec& Pol, const Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
            const Vector_t<double, 3>& Bf, const double& dt, const double& mass,
            const double& charge, const double& anom) const {
        // Thomas-BMT angular velocity vector Omega from lab-frame E, B and lab time, written
        // along the axes of the inputs (Omega is built from dot and cross products only, so
        // any fixed rotation of P, E, B and Pol together gives the same result):
        //   Omega = -(q/m) [ (G + 1/gamma) B
        //                    - (G gamma / (gamma+1)) (beta . B) beta
        //                    - (G + 1/(gamma+1)) (beta x E) / c ]
        //
        // OPALX convention: mass is rest energy in eV (i.e. m*c^2 numerically),
        // so q/m becomes charge*c^2/mass — same idiom used in BorisPusher.
        const double bg2               = dot(P, P);
        const double gamma             = Kokkos::sqrt(1.0 + bg2);
        const Vector_t<double, 3> beta = P / gamma;

        const double betaDotB             = dot(beta, Bf);
        const Vector_t<double, 3> bcrossE = cross(beta, Ef);

        const double k1 = anom + 1.0 / gamma;
        const double k2 = (anom * gamma) / (gamma + 1.0);
        const double k3 = anom + 1.0 / (gamma + 1.0);

        const double prefactor = -charge * Physics::c * Physics::c / mass;
        const Vector_t<double, 3> Omega =
                prefactor * (k1 * Bf - k2 * betaDotB * beta - (k3 / Physics::c) * bcrossE);

        const double omegaMag = Kokkos::sqrt(dot(Omega, Omega));
        const double phi      = omegaMag * dt;

        if (phi < 1.0e-14) {
            return;
        }

        // One rotation by phi about Omega (Rodrigues, exact at any angle for this Omega):
        //   S' = cos(phi) S + sin(phi) (axis x S) + (1 - cos(phi)) (axis . S) axis
        const double cosPhi            = Kokkos::cos(phi);
        const double sinPhi            = Kokkos::sin(phi);
        const Vector_t<double, 3> axis = Omega / omegaMag;

        Vector_t<double, 3> S{
                static_cast<double>(Pol[0]), static_cast<double>(Pol[1]),
                static_cast<double>(Pol[2])};
        const Vector_t<double, 3> axisCrossS = cross(axis, S);
        const double axisDotS                = dot(axis, S);
        S = cosPhi * S + sinPhi * axisCrossS + (1.0 - cosPhi) * axisDotS * axis;

        Pol[0] = S[0];
        Pol[1] = S[1];
        Pol[2] = S[2];
    }
};

#endif
