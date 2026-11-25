//
// Make the per-beam values of a PartBunch from a BEAM object, for tests without a RUN.
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
// All rights reserved
//
#ifndef OPALX_UNIT_TESTS_BUNCH_BEAM_H
#define OPALX_UNIT_TESTS_BUNCH_BEAM_H

#include <cstddef>

#include "PartBunch/PartBunch.h"
#include "Structure/Beam.h"

namespace opalx::test {

    /**
     * @brief The values PartBunch needs for one container, for unit tests that make a bunch
     *        without a RUN: reference particle, PARTICLE name and POLARIZATION from @p beam.
     * @param numParticles      Particles to reserve room for.
     * @param chargePerParticle Macro charge per particle (C).
     * @param massPerParticle   Macro mass per particle (GeV/c^2).
     */
    inline BunchBeam bunchBeam(
            const Beam& beam, std::size_t numParticles, double chargePerParticle = 1.0,
            double massPerParticle = 1.0) {
        BunchBeam b;
        b.chargePerParticle = chargePerParticle;
        b.massPerParticle   = massPerParticle;
        b.reference         = beam.getReference();
        b.particleName      = beam.getParticleName();
        b.hasSpin           = beam.hasPolarization();
        b.numParticles      = numParticles;
        return b;
    }

}  // namespace opalx::test

#endif  // OPALX_UNIT_TESTS_BUNCH_BEAM_H
