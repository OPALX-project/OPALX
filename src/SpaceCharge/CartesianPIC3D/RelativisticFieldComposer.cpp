/**
 * @file RelativisticFieldComposer.cpp
 * @brief Implements CartesianPIC3D field conversion, accumulation, and final gather operations.
 */

#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/FieldMirror.hpp"
#include "Utilities/OpalException.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace opalx::spacecharge {

    FieldCompositionPolicy FieldCompositionPolicy::centreOfMomentum(
            const std::array<double, 3>& momentum, double meanGamma) {
        const double norm      = std::hypot(momentum[0], momentum[1], momentum[2]);
        const double coldGamma = std::hypot(1.0, norm);
        // A small negative energy gap can result from the existing parallel averages. Only
        // admit their binary64 roundoff; do not clip the measured moments or invariant.
        const double roundoff =
                64.0 * std::numeric_limits<double>::epsilon() * std::max(meanGamma, coldGamma);
        const double invariant = 1.0 + (meanGamma - coldGamma) * (meanGamma + coldGamma);
        if (!std::isfinite(meanGamma) || !std::isfinite(norm) || meanGamma < 1.0
            || meanGamma + roundoff < coldGamma || !std::isfinite(invariant) || invariant <= 0.0) {
            throw OpalException(
                    "FieldCompositionPolicy::centreOfMomentum",
                    "The mean particle energy and momentum do not define a physical CM frame.");
        }
        FieldCompositionPolicy result;
        result.meanMomentum = momentum;
        // Equivalent to G/sqrt(D), using G^2 = D + |mean(p)|^2. This keeps gamma >= 1
        // for stationary and nearly stationary warm bunches without clipping moments.
        result.gamma             = std::hypot(1.0, norm / std::sqrt(invariant));
        result.meanParticleGamma = meanGamma;
        return result;
    }

    LorentzBoost LorentzBoost::from(const FieldCompositionPolicy& policy) {
        LorentzBoost boost;
        const double velocityGamma =
                policy.meanParticleGamma > 0.0 ? policy.meanParticleGamma : policy.gamma;
        const double invGamma = (velocityGamma > 0.0) ? (1.0 / velocityGamma) : 0.0;
        Vector meanMomentum(0.0);
        for (unsigned dimension = 0; dimension < 3; ++dimension) {
            meanMomentum[dimension] = policy.meanMomentum[dimension];
        }

        boost.velocity            = Physics::c * meanMomentum * invGamma;
        const double velocityNorm = Kokkos::sqrt(boost.velocity.dot(boost.velocity));
        boost.velocityUnit =
                (velocityNorm > 0.0) ? Vector(boost.velocity / velocityNorm) : Vector(0.0);
        boost.gamma          = policy.gamma;
        boost.gammaMinusOne  = policy.gamma - 1.0;
        boost.magneticFactor = policy.magneticSign * (policy.gamma / (Physics::c * Physics::c));
        return boost;
    }

    void RelativisticFieldComposer::clearAccumulation(FieldStorage& fieldStorage) const {
        fieldStorage.accumulatedElectricField() = 0.0;
        fieldStorage.accumulatedMagneticField() = 0.0;
    }

    void RelativisticFieldComposer::accumulate(
            FieldStorage& fieldStorage, const Policy& policy) const {
        const LorentzBoost boost = LorentzBoost::from(policy);

        const VectorField& sourceField     = fieldStorage.electricField();
        VectorField& electricTotal         = fieldStorage.accumulatedElectricField();
        VectorField& magneticTotal         = fieldStorage.accumulatedMagneticField();
        using ReadOnlyFieldView            = VectorField::view_type::const_type;
        const ReadOnlyFieldView sourceView = sourceField.getView();
        auto electricTotalView             = electricTotal.getView();
        auto magneticTotalView             = magneticTotal.getView();

        if (policy.sourceRule == FieldSourceRule::Direct) {
            ippl::parallel_for(
                    "RelativisticFieldComposer::accumulateDirect",
                    sourceField.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(const ippl::RangePolicy<3>::index_array_type& idx) {
                        const Vector electricPrime = apply(sourceView, idx);
                        Vector electricLab;
                        Vector magneticLab;
                        boost.toLab(electricPrime, electricLab, magneticLab);
                        Vector accumulatedElectric = apply(electricTotalView, idx);
                        Vector accumulatedMagnetic = apply(magneticTotalView, idx);
                        accumulatedElectric += electricLab;
                        accumulatedMagnetic += magneticLab;
                        apply(electricTotalView, idx) = accumulatedElectric;
                        apply(magneticTotalView, idx) = accumulatedMagnetic;
                    });
            return;
        }

        if (policy.sourceRule != FieldSourceRule::ShiftedGreenImageZ) {
            throw OpalException(
                    "RelativisticFieldComposer::accumulate", "Unknown backend-field source rule.");
        }

        VectorField& mirroredField             = fieldStorage.mirrorScratchFor(sourceField);
        IpplTimings::TimerRef mirrorFieldTimer = IpplTimings::getTimer("mirrorField");
        IpplTimings::startTimer(mirrorFieldTimer);
        opalx::detail::mirrorField(sourceField, mirroredField, 2);
        IpplTimings::stopTimer(mirrorFieldTimer);

        const ReadOnlyFieldView mirroredView = mirroredField.getView();
        const int flipAxis                   = 2;
        const int flipAxisValue              = flipAxis;
        ippl::parallel_for(
                "RelativisticFieldComposer::accumulateShiftedGreenImageZ",
                sourceField.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<3>::index_array_type& idx) {
                    Vector electricPrime = mirroredView(idx[0], idx[1], idx[2]);
                    for (unsigned dimension = 0; dimension < 3; ++dimension) {
                        if (static_cast<int>(dimension) != flipAxisValue) {
                            electricPrime[dimension] = -electricPrime[dimension];
                        }
                    }

                    Vector electricLab;
                    Vector magneticLab;
                    boost.toLab(electricPrime, electricLab, magneticLab);
                    Vector accumulatedElectric = apply(electricTotalView, idx);
                    Vector accumulatedMagnetic = apply(magneticTotalView, idx);
                    accumulatedElectric += electricLab;
                    accumulatedMagnetic += magneticLab;
                    apply(electricTotalView, idx) = accumulatedElectric;
                    apply(magneticTotalView, idx) = accumulatedMagnetic;
                });
    }

    void RelativisticFieldComposer::addParticleContribution(
            const VectorView& electricRest, VectorAttribute& electricDestination,
            VectorAttribute& magneticDestination, std::size_t localCount,
            const Policy& policy) const {
        if (policy.sourceRule != FieldSourceRule::Direct) {
            throw OpalException(
                    "RelativisticFieldComposer::addParticleContribution",
                    "Particle contributions support only the direct source rule.");
        }
        if (electricRest.extent(0) < localCount) {
            throw OpalException(
                    "RelativisticFieldComposer::addParticleContribution",
                    "The rest-frame particle field is shorter than the local particle count.");
        }
        const LorentzBoost boost            = LorentzBoost::from(policy);
        const VectorView::const_type source = electricRest;
        auto electricView                   = electricDestination.getView();
        auto magneticView                   = magneticDestination.getView();
        Kokkos::parallel_for(
                "RelativisticFieldComposer::addParticleContribution",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, localCount),
                KOKKOS_LAMBDA(const std::size_t index) {
                    Vector electricLab;
                    Vector magneticLab;
                    boost.toLab(source(index), electricLab, magneticLab);
                    electricView(index) += electricLab;
                    magneticView(index) += magneticLab;
                });
    }

    void RelativisticFieldComposer::gatherElectrostatic(
            ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& destination,
            const PositionAttribute& positions, FieldStorage& fieldStorage) const {
        particleMeshTransfer.gatherVector(
                destination, fieldStorage.electricField(), positions,
                ParticleMeshTransfer::GatherMode::Replace);
    }

    void RelativisticFieldComposer::gatherAccumulated(
            ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& electricDestination,
            VectorAttribute& magneticDestination, const PositionAttribute& positions,
            FieldStorage& fieldStorage) const {
        particleMeshTransfer.gatherVector(
                electricDestination, fieldStorage.accumulatedElectricField(), positions,
                ParticleMeshTransfer::GatherMode::Replace);
        particleMeshTransfer.gatherVector(
                magneticDestination, fieldStorage.accumulatedMagneticField(), positions,
                ParticleMeshTransfer::GatherMode::Replace);
    }

}  // namespace opalx::spacecharge
