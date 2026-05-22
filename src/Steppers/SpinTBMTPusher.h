//
// Class SpinTBMTPusher
//   Thomas-BMT spin precession integrator using analytic Rodrigues rotation.
//   Preserves |Pol| exactly per step; sub-steps when the per-step rotation
//   angle exceeds a small-angle threshold.
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

    /// Sub-step threshold: maximum rotation angle (rad) per Rodrigues call.
    /// Above this we split dt into uniform sub-steps to keep the rotation accurate.
    static constexpr double maxAngle_m = 0.1;

    /// Per-particle Thomas-BMT update.
    ///
    /// Inputs (all lab-frame, in OPALX conventions):
    ///   P       : momentum in beta*gamma (dimensionless)
    ///   Ef, Bf  : electric and magnetic field at the particle position [V/m, T]
    ///   dt      : timestep [s]
    ///   mass    : rest energy [eV]  (m * c^2)
    ///   charge  : particle charge [C]
    ///   anom    : magnetic moment anomaly G = (g - 2) / 2 [dimensionless]
    /// In/out:
    ///   Pol     : polarization 3-vector (rest-frame Pauli expectations along lab axes,
    ///             |Pol| in [0, 1]). Updated in place.
    template <typename SpinVec>
    KOKKOS_INLINE_FUNCTION void evolve(
            SpinVec& Pol, const Vector_t<double, 3>& P, const Vector_t<double, 3>& Ef,
            const Vector_t<double, 3>& Bf, const double& dt, const double& mass,
            const double& charge, const double& anom) const {
        // Thomas-BMT angular velocity vector Omega in the lab frame, lab time:
        //   Omega = -(q/m) [ (G + 1/gamma) B
        //                    - (G gamma / (gamma+1)) (beta . B) beta
        //                    - (G + 1/(gamma+1)) (beta x E) / c ]
        //
        // OPALX convention: mass is rest energy in eV (i.e. m*c^2 numerically),
        // so q/m becomes charge*c^2/mass — same idiom used in BorisPusher.
        const double bg2   = dot(P, P);
        const double gamma = Kokkos::sqrt(1.0 + bg2);
        const Vector_t<double, 3> beta = P / gamma;

        const double betaDotB = dot(beta, Bf);
        const Vector_t<double, 3> bcrossE = cross(beta, Ef);

        const double k1 = anom + 1.0 / gamma;
        const double k2 = (anom * gamma) / (gamma + 1.0);
        const double k3 = anom + 1.0 / (gamma + 1.0);

        const double prefactor = -charge * Physics::c * Physics::c / mass;
        const Vector_t<double, 3> Omega =
                prefactor
                * (k1 * Bf - k2 * betaDotB * beta - (k3 / Physics::c) * bcrossE);

        const double omegaMag = Kokkos::sqrt(dot(Omega, Omega));
        const double phi      = omegaMag * dt;

        if (phi < 1.0e-14) {
            return;
        }

        // Sub-stepping for large rotations: split phi into N steps of size <= maxAngle_m.
        // Rodrigues' formula is exact at any angle, but float precision and the
        // assumption of constant Omega across the step degrade for large phi.
        const int nSub = static_cast<int>(Kokkos::ceil(phi / maxAngle_m));
        const double dphi = phi / static_cast<double>(nSub);
        const double cosDphi = Kokkos::cos(dphi);
        const double sinDphi = Kokkos::sin(dphi);

        const Vector_t<double, 3> axis = Omega / omegaMag;

        Vector_t<double, 3> S{static_cast<double>(Pol[0]), static_cast<double>(Pol[1]),
                              static_cast<double>(Pol[2])};
        for (int i = 0; i < nSub; ++i) {
            // Rodrigues: S' = cos(dphi) S + sin(dphi) (axis x S) + (1-cos(dphi))(axis.S) axis
            const Vector_t<double, 3> axisCrossS = cross(axis, S);
            const double axisDotS                = dot(axis, S);
            S = cosDphi * S + sinDphi * axisCrossS + (1.0 - cosDphi) * axisDotS * axis;
        }

        Pol[0] = static_cast<float>(S[0]);
        Pol[1] = static_cast<float>(S[1]);
        Pol[2] = static_cast<float>(S[2]);
    }
};

#endif
