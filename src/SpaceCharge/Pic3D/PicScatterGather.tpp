/**
 * @file PicScatterGather.tpp
 * @brief Implements Cartesian PIC CIC deposition and gathering.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_TPP
#define OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_TPP

#include "Interpolation/CIC.h"
#include "Utilities/OpalException.h"

#include <exception>
#include <functional>
#include <numeric>
#include <string>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::depositCharge(
            ParticleContainer& particles, Workspace& workspace, DepositKind depositKind,
            const Selection& selection, const ChargeNormalization& normalization,
            const ImagePolicy& imagePolicy) const {
        validateSelection(particles, selection);

        if (normalization.gamma <= 0.0) {
            throw OpalException(
                    "PicScatterGather::depositCharge",
                    "Charge-density normalization requires a positive gamma.");
        }
        if (depositKind == DepositKind::Image && !imagePolicy.enabled) {
            throw OpalException(
                    "PicScatterGather::depositCharge",
                    "An image-only deposit requires an enabled image policy.");
        }

        ScalarField& rho             = workspace.chargeDensity();
        PositionAttribute& positions = particles.R;
        rho                          = 0.0;

        if (depositKind == DepositKind::Primary || depositKind == DepositKind::PrimaryAndImage) {
            scatterScaledTimeStep(particles, positions, rho, selection);
        }

        if (depositKind == DepositKind::Image
            || (depositKind == DepositKind::PrimaryAndImage && imagePolicy.enabled)) {
            depositImage(particles, positions, rho, selection, imagePolicy);
        }

        normalizeChargeDensity(workspace, normalization);
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::gatherVector(
            VectorAttribute& destination, VectorField& source, const PositionAttribute& positions,
            GatherMode mode) const {
        const bool addToDestination = mode == GatherMode::Add;
        gather(destination, source, positions, addToDestination);
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::scatterScaledTimeStep(
            ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection) const {
        particles.scaleDtByCharge();
        try {
            scatterCurrentWeights(particles, positions, rho, selection);
        } catch (...) {
            const std::exception_ptr original = std::current_exception();
            try {
                particles.unscaleDtByCharge();
            } catch (...) {
            }
            std::rethrow_exception(original);
        }
        particles.unscaleDtByCharge();
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::applyImageTransform(
            ParticleContainer& particles, PositionAttribute& positions, const Selection& selection,
            double planeZ) const {
        reflectPositions(positions, selection, planeZ);
        try {
            flipChargeSign(particles, selection);
        } catch (...) {
            const std::exception_ptr original = std::current_exception();
            try {
                reflectPositions(positions, selection, planeZ);
            } catch (...) {
            }
            std::rethrow_exception(original);
        }
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::restoreImageTransform(
            ParticleContainer& particles, PositionAttribute& positions, const Selection& selection,
            double planeZ) const {
        applyImageTransform(particles, positions, selection, planeZ);
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::reflectPositions(
            PositionAttribute& positions, const Selection& selection, double planeZ) const {
        const auto positionView  = positions.getView();
        const RangePolicy policy = selection.policy();
        const Hash hash          = selection.hash();
        const bool useHash       = selection.kind() == Selection::Kind::Indexed;

        Kokkos::parallel_for(
                "PicScatterGather::reflectImagePositions", policy,
                KOKKOS_LAMBDA(const size_type selectionIndex) {
                    const size_type particleIndex =
                            useHash ? static_cast<size_type>(hash(selectionIndex)) : selectionIndex;
                    positionView(particleIndex)[2] = 2.0 * planeZ - positionView(particleIndex)[2];
                });
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::flipChargeSign(
            ParticleContainer& particles, const Selection& selection) const {
        const auto chargeView = particles.getQView();
        if (particles.getQMStorageMode() == ParticleContainer::QMStorageMode::Attributes) {
            const RangePolicy policy = selection.policy();
            const Hash hash          = selection.hash();
            const bool useHash       = selection.kind() == Selection::Kind::Indexed;
            Kokkos::parallel_for(
                    "PicScatterGather::flipSelectedChargeSigns", policy,
                    KOKKOS_LAMBDA(const size_type selectionIndex) {
                        const size_type particleIndex =
                                useHash ? static_cast<size_type>(hash(selectionIndex))
                                        : selectionIndex;
                        chargeView(particleIndex) = -chargeView(particleIndex);
                    });
            return;
        }

        Kokkos::parallel_for(
                "PicScatterGather::flipSharedChargeSign", size_type(1),
                KOKKOS_LAMBDA(const size_type) { chargeView(0) = -chargeView(0); });
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::validateSelection(
            const ParticleContainer& particles, const Selection& selection) const {
        const size_type begin = static_cast<size_type>(selection.policy().begin());
        const size_type end   = static_cast<size_type>(selection.policy().end());
        if (begin > end) {
            throw OpalException(
                    "PicScatterGather::validateSelection",
                    "The particle selection begins after it ends.");
        }

        if (selection.kind() == Selection::Kind::Direct) {
            if (end > static_cast<size_type>(particles.getLocalNum())) {
                throw OpalException(
                        "PicScatterGather::validateSelection",
                        "A direct particle selection exceeds the local particle count.");
            }
            return;
        }

        if (end > static_cast<size_type>(selection.hash().extent(0))) {
            throw OpalException(
                    "PicScatterGather::validateSelection",
                    "An indexed particle selection exceeds its hash extent.");
        }
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::scatterCurrentWeights(
            const ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection) const {
        const size_type begin = static_cast<size_type>(selection.policy().begin());
        const size_type end   = static_cast<size_type>(selection.policy().end());
        const bool allLocal   = selection.kind() == Selection::Kind::Direct && begin == 0
                              && end == static_cast<size_type>(particles.getLocalNum());

        if (allLocal) {
            scatter(particles.dt, rho, positions);
            return;
        }
        if (selection.kind() == Selection::Kind::Indexed) {
            scatter(particles.dt, rho, positions, selection.policy(), selection.hash());
            return;
        }
        scatter(particles.dt, rho, positions, selection.policy());
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::depositImage(
            ParticleContainer& particles, PositionAttribute& positions, ScalarField& rho,
            const Selection& selection, const ImagePolicy& imagePolicy) const {
        applyImageTransform(particles, positions, selection, imagePolicy.planeZ);
        try {
            scatterScaledTimeStep(particles, positions, rho, selection);
        } catch (...) {
            const std::exception_ptr original = std::current_exception();
            try {
                restoreImageTransform(particles, positions, selection, imagePolicy.planeZ);
            } catch (...) {
            }
            std::rethrow_exception(original);
        }
        restoreImageTransform(particles, positions, selection, imagePolicy.planeZ);
    }

    template <typename T, unsigned Dim>
    void PicScatterGather<T, Dim>::normalizeChargeDensity(
            Workspace& workspace, const ChargeNormalization& normalization) const {
        double normalizer = normalization.timeStep;
        if (normalization.normalizeByCellVolume) {
            const auto& spacing = workspace.spacing();
            const double cellVolume =
                    std::reduce(spacing.begin(), spacing.end(), 1.0, std::multiplies<double>());
            normalizer *= cellVolume;
        }

        double shift = 0.0;
        if (normalization.subtractNeutralizingBackground) {
            const auto& lower = workspace.lower();
            const auto& upper = workspace.upper();
            double volume     = 1.0;
            for (size_type d = 0; d < Dim; ++d) {
                volume *= upper[d] - lower[d];
            }
            shift = -(normalization.selectedCharge / volume) * normalization.couplingConstant
                    / normalization.gamma;
        }

        normalizer *= normalization.gamma;
        ScalarField& rho = workspace.chargeDensity();
        rho              = rho * (normalization.couplingConstant / normalizer) + shift;
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_SCATTER_GATHER_TPP
