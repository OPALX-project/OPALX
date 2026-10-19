/**
 * @file RelativisticFieldComposer.cpp
 * @brief Implements CartesianPIC3D field conversion, accumulation, and final gather operations.
 */

#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/FieldMirror.hpp"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    LorentzBoost LorentzBoost::from(const FieldCompositionPolicy& policy) {
        LorentzBoost boost;
        const double invGamma = (policy.gamma > 0.0) ? (1.0 / policy.gamma) : 0.0;
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

    void RelativisticFieldComposer::gatherConverted(
            ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& electricDestination,
            VectorAttribute& magneticDestination, const PositionAttribute& positions,
            FieldStorage& fieldStorage, std::size_t localCount, const Policy& policy) const {
        if (policy.sourceRule != FieldSourceRule::Direct) {
            throw OpalException(
                    "RelativisticFieldComposer::gatherConverted",
                    "Per-particle conversion supports only the direct source rule.");
        }
        particleMeshTransfer.gatherVector(
                electricDestination, fieldStorage.electricField(), positions,
                ParticleMeshTransfer::GatherMode::Replace);

        const LorentzBoost boost = LorentzBoost::from(policy);
        auto electricView        = electricDestination.getView();
        auto magneticView        = magneticDestination.getView();
        Kokkos::parallel_for(
                "RelativisticFieldComposer::convertGathered",
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, localCount),
                KOKKOS_LAMBDA(const std::size_t index) {
                    const Vector electricRest = electricView(index);
                    Vector electricLab;
                    Vector magneticLab;
                    boost.toLab(electricRest, electricLab, magneticLab);
                    electricView(index) = electricLab;
                    magneticView(index) = magneticLab;
                });
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
