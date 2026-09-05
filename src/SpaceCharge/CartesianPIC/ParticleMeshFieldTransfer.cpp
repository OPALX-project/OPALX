/**
 * @file ParticleMeshFieldTransfer.cpp
 * @brief Implements Cartesian PIC CIC deposition and gathering.
 */

#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"
#include "Interpolation/CIC.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <string>

namespace opalx::spacecharge {

    void ParticleMeshFieldTransfer::depositCharge(
            ParticleContainer& particles, FieldStorage& fieldStorage, DepositKind depositKind,
            const Selection& selection, const ChargeNormalization& normalization,
            const ImagePolicy& imagePolicy) const {
        validateSelection(particles, selection);

        if (normalization.gamma <= 0.0) {
            throw OpalException(
                    "ParticleMeshFieldTransfer::depositCharge",
                    "Charge-density normalization requires a positive gamma.");
        }
        if (depositKind == DepositKind::Image && !imagePolicy.enabled) {
            throw OpalException(
                    "ParticleMeshFieldTransfer::depositCharge",
                    "An image-only deposit requires an enabled image policy.");
        }

        ScalarField& rho             = fieldStorage.chargeDensity();
        PositionAttribute& positions = particles.R;
        rho                          = 0.0;

        if (depositKind == DepositKind::Primary || depositKind == DepositKind::PrimaryAndImage) {
            scatterScaledTimeStep(particles, positions, rho, selection);
        }

        if (depositKind == DepositKind::Image
            || (depositKind == DepositKind::PrimaryAndImage && imagePolicy.enabled)) {
            depositImage(particles, positions, rho, selection, imagePolicy);
        }

        normalizeChargeDensity(fieldStorage, normalization);
    }

    void ParticleMeshFieldTransfer::gatherVector(
            VectorAttribute& destination, VectorField& source, const PositionAttribute& positions,
            GatherMode mode) const {
        const bool addToDestination = mode == GatherMode::Add;
        gather(destination, source, positions, addToDestination);
    }

    void ParticleMeshFieldTransfer::scatterScaledTimeStep(
            ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection) const {
        // IPPL scatter accepts a particle attribute as its weight. Preserve the established
        // convention by temporarily storing dt*Q in dt, then restore dt after the scatter.
        particles.scaleDtByCharge();
        scatterCurrentWeights(particles, positions, rho, selection);
        particles.unscaleDtByCharge();
    }

    void ParticleMeshFieldTransfer::applyImageTransform(
            ParticleContainer& particles, PositionAttribute& positions, const Selection& selection,
            double planeZ) const {
        // Reflecting z and flipping Q are both self-inverse. Keeping them paired allows the same
        // operation to restore the particle state after successful deposition.
        reflectPositions(positions, selection, planeZ);
        flipChargeSign(particles, selection);
    }

    void ParticleMeshFieldTransfer::restoreImageTransform(
            ParticleContainer& particles, PositionAttribute& positions, const Selection& selection,
            double planeZ) const {
        applyImageTransform(particles, positions, selection, planeZ);
    }

    void ParticleMeshFieldTransfer::reflectPositions(
            PositionAttribute& positions, const Selection& selection, double planeZ) const {
        const auto positionView  = positions.getView();
        const RangePolicy policy = selection.policy();
        const Hash hash          = selection.hash();
        const bool useHash       = selection.kind() == Selection::Kind::Indexed;

        Kokkos::parallel_for(
                "ParticleMeshFieldTransfer::reflectImagePositions", policy,
                KOKKOS_LAMBDA(const size_type selectionIndex) {
                    const size_type particleIndex =
                            useHash ? static_cast<size_type>(hash(selectionIndex)) : selectionIndex;
                    positionView(particleIndex)[2] = 2.0 * planeZ - positionView(particleIndex)[2];
                });
    }

    void ParticleMeshFieldTransfer::flipChargeSign(
            ParticleContainer& particles, const Selection& selection) const {
        const auto chargeView = particles.getQView();
        if (particles.getQMStorageMode() == ParticleContainer::QMStorageMode::Attributes) {
            const RangePolicy policy = selection.policy();
            const Hash hash          = selection.hash();
            const bool useHash       = selection.kind() == Selection::Kind::Indexed;
            Kokkos::parallel_for(
                    "ParticleMeshFieldTransfer::flipSelectedChargeSigns", policy,
                    KOKKOS_LAMBDA(const size_type selectionIndex) {
                        const size_type particleIndex =
                                useHash ? static_cast<size_type>(hash(selectionIndex))
                                        : selectionIndex;
                        chargeView(particleIndex) = -chargeView(particleIndex);
                    });
            return;
        }

        Kokkos::parallel_for(
                "ParticleMeshFieldTransfer::flipSharedChargeSign", size_type(1),
                KOKKOS_LAMBDA(const size_type) { chargeView(0) = -chargeView(0); });
    }

    void ParticleMeshFieldTransfer::validateSelection(
            const ParticleContainer& particles, const Selection& selection) const {
        const size_type begin = static_cast<size_type>(selection.policy().begin());
        const size_type end   = static_cast<size_type>(selection.policy().end());
        if (begin > end) {
            throw OpalException(
                    "ParticleMeshFieldTransfer::validateSelection",
                    "The particle selection begins after it ends.");
        }

        if (selection.kind() == Selection::Kind::Direct) {
            if (end > static_cast<size_type>(particles.getLocalNum())) {
                throw OpalException(
                        "ParticleMeshFieldTransfer::validateSelection",
                        "A direct particle selection exceeds the local particle count.");
            }
            return;
        }

        if (end > static_cast<size_type>(selection.hash().extent(0))) {
            throw OpalException(
                    "ParticleMeshFieldTransfer::validateSelection",
                    "An indexed particle selection exceeds its hash extent.");
        }
    }

    void ParticleMeshFieldTransfer::scatterCurrentWeights(
            const ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection) const {
        const size_type begin = static_cast<size_type>(selection.policy().begin());
        const size_type end   = static_cast<size_type>(selection.policy().end());
        const bool allLocal   = selection.kind() == Selection::Kind::Direct && begin == 0
                              && end == static_cast<size_type>(particles.getLocalNum());

        if (allLocal) {
            // Avoid hash or custom-policy indirection for the common one-bin early-emission case.
            scatter(particles.dt, rho, positions);
            return;
        }
        if (selection.kind() == Selection::Kind::Indexed) {
            scatter(particles.dt, rho, positions, selection.policy(), selection.hash());
            return;
        }
        scatter(particles.dt, rho, positions, selection.policy());
    }

    void ParticleMeshFieldTransfer::depositImage(
            ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection, const ImagePolicy& imagePolicy) const {
        applyImageTransform(particles, positions, selection, imagePolicy.planeZ);
        scatterScaledTimeStep(particles, positions, rho, selection);
        restoreImageTransform(particles, positions, selection, imagePolicy.planeZ);
    }

    void ParticleMeshFieldTransfer::normalizeChargeDensity(
            FieldStorage& fieldStorage, const ChargeNormalization& normalization) const {
        // The deposited weight is dt*Q. Divide by the global time step, cell volume when required,
        // and the bin gamma to recover the charge-density convention expected by the Poisson
        // solver. P3M retains volume normalization but deliberately skips the periodic
        // neutralizing background.
        if (!std::isfinite(normalization.timeStep) || normalization.timeStep == 0.0
            || !std::isfinite(normalization.gamma) || normalization.gamma < 1.0) {
            throw OpalException(
                    "ParticleMeshFieldTransfer::normalizeChargeDensity",
                    "Deposition requires a finite nonzero time step and gamma >= 1.");
        }
        double normalizer = normalization.timeStep;
        if (normalization.normalizeByCellVolume) {
            const auto& spacing = fieldStorage.spacing();
            const double cellVolume =
                    std::reduce(spacing.begin(), spacing.end(), 1.0, std::multiplies<double>());
            normalizer *= cellVolume;
        }

        double shift = 0.0;
        if (normalization.subtractNeutralizingBackground) {
            // Periodic FFT requires zero net charge. Subtract only the charge represented by this
            // solve unit so binned and whole-bunch paths use the same physical background.
            const auto& lower = fieldStorage.lower();
            const auto& upper = fieldStorage.upper();
            double volume     = 1.0;
            for (size_type d = 0; d < 3; ++d) {
                volume *= upper[d] - lower[d];
            }
            shift = -(normalization.selectedCharge / volume) * normalization.couplingConstant
                    / normalization.gamma;
        }

        normalizer *= normalization.gamma;
        ScalarField& rho = fieldStorage.chargeDensity();
        rho              = rho * (normalization.couplingConstant / normalizer) + shift;
    }

}  // namespace opalx::spacecharge
