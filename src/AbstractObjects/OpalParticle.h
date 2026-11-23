//
// Class OpalParticle
//   This class represents the canonical coordinates of a particle.
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
#ifndef OPALX_OpalParticle_HH
#define OPALX_OpalParticle_HH

#include "Ippl.h"
#include "OPALTypes.h"

//
#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>

class OpalParticle {
public:
    // Particle coordinate numbers.
    enum { X, Y, L, INVALID };

    /// Constructor.
    //  Construct particle with the given coordinates.
    OpalParticle(
            int64_t id, double x, double px, double y, double py, double z, double pz, double time,
            double q, double m);

    OpalParticle(
            int64_t id, Vector_t<double, 3> const& R, Vector_t<double, 3> const& P, double time,
            double q, double m);

    /// Constructor including the per-particle polarization vector.
    /// Sets the hasPol flag so the loss sink emits the polx/poly/polz columns.
    OpalParticle(
            int64_t id, Vector_t<double, 3> const& R, Vector_t<double, 3> const& P, double time,
            double q, double m, Vector_t<double, 3> const& Pol);

    OpalParticle();

    /// Get the id of the particle
    int64_t getId() const;

    /// Get coordinate.
    //  Access coordinate by index for constant particle.
    double operator[](unsigned int) const;

    /// Get horizontal position in m.
    double getX() const;

    /// Get horizontal momentum (no dimension).
    double getPx() const;

    /// Get vertical displacement in m.
    double getY() const;

    /// Get vertical momentum (no dimension).
    double getPy() const;

    /// Get longitudinal displacement c*t in m.
    double getZ() const;

    /// Get relative momentum error (no dimension).
    double getPz() const;

    /// Get momentum
    const Vector_t<double, 3>& getP() const;

    /// Get time
    double getTime() const;

    /// Get charge in Coulomb
    double getCharge() const;

    /// Get mass in GeV/c^2
    double getMass() const;

    /// Whether a polarization vector was attached to this particle.
    bool hasPol() const;

    /// Get the rest-frame polarization vector (rest-frame components along lab axes).
    /// Only meaningful when hasPol() is true.
    const Vector_t<double, 3>& getPol() const;

private:
    int64_t id_m;
    Vector_t<double, 3> R_m;
    Vector_t<double, 3> P_m;
    double time_m;
    double charge_m;
    double mass_m;
    Vector_t<double, 3> Pol_m{0.0, 0.0, 0.0};
    bool hasPol_m = false;
};

inline int64_t OpalParticle::getId() const { return id_m; }

inline double OpalParticle::operator[](unsigned int i) const {
    PAssert_LT(i, 6u);
    return i % 2 == 0 ? R_m[i / 2] : P_m[i / 2];
}

inline double OpalParticle::getX() const { return R_m[X]; }

inline double OpalParticle::getY() const { return R_m[Y]; }

inline double OpalParticle::getZ() const { return R_m[L]; }

inline double OpalParticle::getPx() const { return P_m[X]; }

inline double OpalParticle::getPy() const { return P_m[Y]; }

inline double OpalParticle::getPz() const { return P_m[L]; }

inline const Vector_t<double, 3>& OpalParticle::getP() const { return P_m; }

inline double OpalParticle::getTime() const { return time_m; }

inline double OpalParticle::getCharge() const { return charge_m; }

inline double OpalParticle::getMass() const { return mass_m; }

inline bool OpalParticle::hasPol() const { return hasPol_m; }

inline const Vector_t<double, 3>& OpalParticle::getPol() const { return Pol_m; }

#endif  // OPALX_OpalParticle_HH