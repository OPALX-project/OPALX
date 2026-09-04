/**
 * @file RelativisticFieldComposer.tpp
 * @brief Implements Cartesian PIC field conversion, accumulation, and final gather operations.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_TPP
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_TPP

#include "Physics/Physics.h"
#include "SpaceCharge/CartesianPIC/FieldMirror.hpp"
#include "Utilities/OpalException.h"

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    void RelativisticFieldComposer<T, Dim>::clearAccumulation(FieldStorage& fieldStorage) const {
        fieldStorage.accumulatedElectricField() = 0.0;
        fieldStorage.accumulatedMagneticField() = 0.0;
    }

    template <typename T, unsigned Dim>
    void RelativisticFieldComposer<T, Dim>::accumulate(
            FieldStorage& fieldStorage, const Policy& policy) const {
        const double gammaBin   = policy.gamma;
        const double bFieldSign = policy.magneticSign;
        const double invGamma   = (gammaBin > 0.0) ? (1.0 / gammaBin) : 0.0;
        Vector meanMomentum(0.0);
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            meanMomentum[dimension] = policy.meanMomentum[dimension];
        }

        const Vector velocity      = Physics::c * meanMomentum * invGamma;
        const double velocityNorm  = Kokkos::sqrt(velocity.dot(velocity));
        const Vector velocityUnit  = (velocityNorm > 0.0) ? (velocity / velocityNorm) : Vector(0.0);
        const double gammaMinusOne = gammaBin - 1.0;
        const double gammaOverCSq  = gammaBin / (Physics::c * Physics::c);

        const VectorField& sourceField = fieldStorage.electricField();
        VectorField& electricTotal     = fieldStorage.accumulatedElectricField();
        VectorField& magneticTotal     = fieldStorage.accumulatedMagneticField();
        auto sourceView                = sourceField.getView();
        auto electricTotalView         = electricTotal.getView();
        auto magneticTotalView         = magneticTotal.getView();

        if (policy.sourceRule == FieldSourceRule::Direct) {
            ippl::parallel_for(
                    "RelativisticFieldComposer::accumulateDirect",
                    sourceField.getFieldRangePolicy(),
                    KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                        Vector electricPrime            = apply(sourceView, idx);
                        const T electricDotVelocityUnit = electricPrime.dot(velocityUnit);
                        Vector electricLab =
                                gammaBin * electricPrime
                                - gammaMinusOne * electricDotVelocityUnit * velocityUnit;
                        Vector magneticLab =
                                bFieldSign * gammaOverCSq * cross(velocity, electricPrime);
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
        try {
            opalx::detail::mirrorField(sourceField, mirroredField, Dim - 1);
        } catch (...) {
            IpplTimings::stopTimer(mirrorFieldTimer);
            throw;
        }
        IpplTimings::stopTimer(mirrorFieldTimer);

        auto mirroredView       = mirroredField.getView();
        const int flipAxis      = static_cast<int>(Dim) - 1;
        const int flipAxisValue = flipAxis;
        ippl::parallel_for(
                "RelativisticFieldComposer::accumulateShiftedGreenImageZ",
                sourceField.getFieldRangePolicy(),
                KOKKOS_LAMBDA(const ippl::RangePolicy<Dim>::index_array_type& idx) {
                    Vector electricPrime = mirroredView(idx[0], idx[1], idx[2]);
                    for (unsigned dimension = 0; dimension < Dim; ++dimension) {
                        if (static_cast<int>(dimension) != flipAxisValue) {
                            electricPrime[dimension] = -electricPrime[dimension];
                        }
                    }

                    const T electricDotVelocityUnit = electricPrime.dot(velocityUnit);
                    Vector electricLab              = gammaBin * electricPrime
                                         - gammaMinusOne * electricDotVelocityUnit * velocityUnit;
                    Vector magneticLab = bFieldSign * gammaOverCSq * cross(velocity, electricPrime);
                    Vector accumulatedElectric = apply(electricTotalView, idx);
                    Vector accumulatedMagnetic = apply(magneticTotalView, idx);
                    accumulatedElectric += electricLab;
                    accumulatedMagnetic += magneticLab;
                    apply(electricTotalView, idx) = accumulatedElectric;
                    apply(magneticTotalView, idx) = accumulatedMagnetic;
                });
    }

    template <typename T, unsigned Dim>
    void RelativisticFieldComposer<T, Dim>::gatherElectrostatic(
            ParticleMeshTransfer& particleMeshTransfer, VectorAttribute& destination,
            const PositionAttribute& positions, FieldStorage& fieldStorage) const {
        particleMeshTransfer.gatherVector(
                destination, fieldStorage.electricField(), positions,
                ParticleMeshTransfer::GatherMode::Replace);
    }

    template <typename T, unsigned Dim>
    void RelativisticFieldComposer<T, Dim>::gatherAccumulated(
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

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_RELATIVISTIC_FIELD_COMPOSER_TPP
