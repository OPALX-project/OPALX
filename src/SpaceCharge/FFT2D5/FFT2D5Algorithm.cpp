/**
 * @file FFT2D5Algorithm.cpp
 * @brief Implements the independent reference-path FFT2D5 space-charge algorithm.
 */

#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"

#include "PartBunch/BunchStateHandler.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include <Kokkos_NumericTraits.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

namespace opalx::spacecharge {

    FFT2D5Algorithm::FFT2D5Algorithm(
            FFT2D5Config config, std::span<ParticleContainer* const> particles,
            std::shared_ptr<const BunchStateHandler> bunchState)
        : config_m(std::move(config)), bunchState_m(std::move(bunchState)) {
        validateSpaceChargeConfig(SpaceChargeConfig(config_m));
        particles_m.reserve(particles.size());
        for (ParticleContainer* particleContainer : particles) {
            if (particleContainer == nullptr) {
                throw OpalException(
                        "FFT2D5Algorithm::FFT2D5Algorithm",
                        "FFT2D5 cannot borrow a null particle container.");
            }
            particles_m.push_back(particleContainer);
        }
        if (particles_m.empty()) {
            throw OpalException(
                    "FFT2D5Algorithm::FFT2D5Algorithm",
                    "FFT2D5 requires at least one particle container.");
        }
        if (bunchState_m == nullptr) {
            throw OpalException(
                    "FFT2D5Algorithm::FFT2D5Algorithm", "The bunch state handler is null.");
        }
    }

    SpaceChargeSolveResult FFT2D5Algorithm::solve(const SpaceChargeSolveContext& context) {
        if (context.trackingActive().size() != particles_m.size()) {
            throw OpalException(
                    "FFT2D5Algorithm::solve", "The container activity set has the wrong size.");
        }
        if (context.stepState().mpiSize != 1 || ippl::Comm->size() != 1) {
            throw OpalException(
                    "FFT2D5Algorithm::solve",
                    "FFT2D5 currently supports only one MPI rank. Distributed fields and ORB "
                    "load balancing are not implemented for this solver.");
        }
        if (bunchState_m->fixedCartesianDomain().has_value()) {
            throw OpalException(
                    "FFT2D5Algorithm::solve", "FFT2D5 does not support a fixed Cartesian domain.");
        }
        ensureInitialized();

        // Design-path points are fixed in lab coordinates; particle arrays live in each
        // container's moving reference frame. Use the primary container to define the
        // common aligned solve frame, then map both the path and every container into it.
        const auto labToSolve = context.stepState().frames.trackerToSolve
                                * particles_m.front()->getToLabTrafo().inverted();
        Kokkos::deep_copy(solveReferencePath_m, referencePath_m->deviceView());
        labToSolve.transformBunchTo(solveReferencePath_m, solveReferencePath_m.extent(0));

        for (std::size_t index = 0; index < particles_m.size(); ++index) {
            if (context.trackingActive()[index] != 0) {
                auto& particles           = *particles_m[index];
                const auto trackerToSolve = labToSolve * particles.getToLabTrafo();
                clearSelfFields(particles);
                enterSolveFrame({trackerToSolve, trackerToSolve.inverted()}, particles);
            }
        }
        const auto restoreFrames = [&] {
            for (std::size_t index = 0; index < particles_m.size(); ++index) {
                if (context.trackingActive()[index] == 0) {
                    continue;
                }
                auto& particles           = *particles_m[index];
                const auto trackerToSolve = labToSolve * particles.getToLabTrafo();
                leaveSolveFrame({trackerToSolve, trackerToSolve.inverted()}, particles);
                particles.markMomentsDirty();
                particles.updateMoments();
            }
        };
        try {
            doRunSolver(context);
        } catch (...) {
            restoreFrames();
            throw;
        }
        restoreFrames();
        SpaceChargeSolveResult result;
        result.backendSolves = fieldStorage_m->backendSolves();
        return result;
    }

    void FFT2D5Algorithm::preparePathGeometryAndProjections(
            const SpaceChargeSolveContext& context) {
        const auto ref = solveReferencePath_m;
        if (ref.extent(0) < 2) return;
        const auto segments = ref.extent(0) - 1;
        if (segmentGeometry_m.extent(0) != segments) {
            segmentGeometry_m = SegmentGeometryView_t("FFT2D5SegmentGeometry", segments);
        }
        const auto geometry = segmentGeometry_m;
        // This short serial prefix preserves increasing segment summation order.
        Kokkos::parallel_for(
                "FFT2D5::preparePathGeometry", 1, KOKKOS_LAMBDA(const size_t) {
                    T prefix = 0.0;
                    for (size_t i = 0; i < segments; ++i) {
                        const Vector3D_t direction = ref(i + 1) - ref(i);
                        const T length             = euclidean_norm(direction);
                        for (unsigned d = 0; d < 3; ++d)
                            geometry(i, d) = direction[d];
                        geometry(i, 3) = direction.dot(direction);
                        geometry(i, 4) = length;
                        geometry(i, 5) = prefix;
                        prefix += length;
                    }
                });
        Kokkos::fence();
        projectedPositions_m.resize(particles_m.size());
        projectedSegments_m.resize(particles_m.size());
        for (size_t container = 0; container < particles_m.size(); ++container) {
            if (context.trackingActive()[container] == 0) continue;
            const auto count = particles_m[container]->getLocalNum();
            if (projectedPositions_m[container].extent(0) != count) {
                projectedPositions_m[container] = VectorView_t("FFT2D5ProjectedPositions", count);
                projectedSegments_m[container] =
                        SegmentIndexView_t("FFT2D5ProjectedSegments", count);
            }
            const auto positions = particles_m[container]->R.getView();
            const auto invalid   = particles_m[container]->InvalidMask.getView();
            const auto projected = projectedPositions_m[container];
            const auto indices   = projectedSegments_m[container];
            Kokkos::parallel_for(
                    "FFT2D5::prepareParticleProjections", count, KOKKOS_LAMBDA(const size_t n) {
                        if (invalid(n)) {
                            projected(n) = Vector3D_t(0.0);
                            indices(n)   = segments;
                            return;
                        }
                        const Vector3D_t position = positions(n);
                        T bestDistance2           = Kokkos::Experimental::finite_max_v<T>;
                        size_t best               = 0;
                        T bestFraction            = 0.0;
                        Vector3D_t closest(0.0);
                        for (size_t i = 0; i < segments; ++i) {
                            const T length2 = geometry(i, 3);
                            if (length2 <= 0.0) continue;
                            const Vector3D_t direction(
                                    geometry(i, 0), geometry(i, 1), geometry(i, 2));
                            const Vector3D_t start = ref(i);
                            const T fraction       = Kokkos::clamp(
                                    direction.dot(position - start) / length2, 0.0, 1.0);
                            const Vector3D_t point      = start + fraction * direction;
                            const Vector3D_t difference = position - point;
                            const T distance2           = difference.dot(difference);
                            if (distance2 < bestDistance2) {
                                best          = i;
                                bestFraction  = fraction;
                                closest       = point;
                                bestDistance2 = distance2;
                            }
                        }
                        const Vector3D_t tangent =
                                Vector3D_t(geometry(best, 0), geometry(best, 1), geometry(best, 2))
                                / geometry(best, 4);
                        const Vector3D_t normal(0.0, 1.0, 0.0);
                        const Vector3D_t binormal = ippl::cross(normal, tangent);
                        const Vector3D_t offset   = position - closest;
                        projected(n)              = Vector3D_t(
                                offset.dot(binormal), offset.dot(normal),
                                geometry(best, 5) + bestFraction * geometry(best, 4));
                        indices(n) = best;
                    });
        }
        Kokkos::fence();
    }

    std::size_t FFT2D5Algorithm::sliceCount() const {
        return fieldStorage_m == nullptr ? 0 : fieldStorage_m->meshSize()[2];
    }

    const ReferencePath::View& FFT2D5Algorithm::referencePathView() const {
        if (referencePath_m == nullptr) {
            throw OpalException(
                    "FFT2D5Algorithm::referencePathView", "FFT2D5 has not been initialized.");
        }
        return referencePath_m->deviceView();
    }

    void FFT2D5Algorithm::ensureInitialized() {
        if (fieldStorage_m != nullptr) {
            return;
        }

        // Publish dependent storage together only after path loading and allocation succeed.
        referencePath_m =
                std::make_unique<ReferencePath>(ReferencePath::load(config_m.referencePathFile));
        solveReferencePath_m = ReferenceView_t(
                "FFT2D5SolveReferencePath", referencePath_m->deviceView().extent(0));
        Kokkos::deep_copy(solveReferencePath_m, referencePath_m->deviceView());
        fieldStorage_m = std::make_unique<FFT2D5FieldStorage>(config_m, referencePath_m->length());
        lineDensity_m  = LineDensityView_t(
                "FFT2D5LineDensity", fieldStorage_m->meshSize()[2] + LineDensityGhostCells);
        lineDensityGradient_m =
                LineDensityView_t("FFT2D5LineDensityGradient", fieldStorage_m->meshSize()[2]);
    }

}  // namespace opalx::spacecharge
