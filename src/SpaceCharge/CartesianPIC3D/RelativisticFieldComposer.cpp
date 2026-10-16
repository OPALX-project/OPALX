/**
 * @file RelativisticFieldComposer.cpp
 * @brief Implements CartesianPIC3D field conversion, accumulation, and final gather operations.
 */

#include "SpaceCharge/CartesianPIC3D/RelativisticFieldComposer.h"
#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC3D/FieldMirror.hpp"
#include "SpaceCharge/RestFrame.h"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    void RelativisticFieldComposer::clearAccumulation(FieldStorage& fieldStorage) const {
        fieldStorage.accumulatedElectricField() = 0.0;
        fieldStorage.accumulatedMagneticField() = 0.0;
    }

    void RelativisticFieldComposer::accumulate(
            FieldStorage& fieldStorage, const Policy& policy) const {
        const double bFieldSign    = policy.magneticSign;
        const RestFrameBoost boost = RestFrameBoost::fromMoments(policy.meanMomentum, policy.gamma);

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
                        const Vector electricLab   = boost.electric(electricPrime);
                        const Vector magneticLab   = boost.magnetic(electricPrime, bFieldSign);
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

                    const Vector electricLab = boost.electric(electricPrime);
                    const Vector magneticLab = boost.magnetic(electricPrime, bFieldSign);
                    Vector accumulatedElectric = apply(electricTotalView, idx);
                    Vector accumulatedMagnetic = apply(magneticTotalView, idx);
                    accumulatedElectric += electricLab;
                    accumulatedMagnetic += magneticLab;
                    apply(electricTotalView, idx) = accumulatedElectric;
                    apply(magneticTotalView, idx) = accumulatedMagnetic;
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
