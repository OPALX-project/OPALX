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
#include <limits>
#include <numeric>
#include <utility>

namespace opalx::spacecharge {

    FFT2D5Algorithm::FFT2D5Algorithm(
            FFT2D5Config config, std::span<ParticleContainer* const> particles,
            std::shared_ptr<const BunchStateHandler> bunchState)
        : config_m(std::move(config)), bunchState_m(std::move(bunchState)) {
        validateSpaceChargeConfig(SpaceChargeConfig(config_m));
        fft2d5_detail::registerTimers();
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
        if (ippl::Comm->size() > 1 && !config_m.grid.decomposition[2]) {
            throw OpalException("FFT2D5Algorithm::solve", "Multiple ranks require PARFFTZ=TRUE.");
        }
        if (bunchState_m->fixedCartesianDomain().has_value()) {
            throw OpalException(
                    "FFT2D5Algorithm::solve", "FFT2D5 does not support a fixed Cartesian domain.");
        }
        if (size_t(ippl::Comm->size()) > config_m.grid.meshSize[2])
            throw OpalException(
                    "FFT2D5Algorithm::solve", "Need at least one z slice per MPI rank.");
        ensureInitialized();

        SpaceChargeSolveResult result;
        fft2d5_detail::timed(firstPhysicalSolve_m ? "2D5 First Solve" : "2D5 Solve", [&] {
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
            result.backendSolves = fieldStorage_m->backendSolves();
            if (ippl::Comm->size() > 1) {
                const auto local = result.backendSolves;
                ippl::Comm->allreduce(local, result.backendSolves, 1, std::plus<size_t>());
            }
        });
        firstPhysicalSolve_m = false;
        return result;
    }

    void FFT2D5Algorithm::migrateToSlabs(const SpaceChargeSolveContext& context) {
        opalx::spacecharge::fft2d5_detail::timed("2D5 Migration", [&] {
            const auto origin = fieldStorage_m->origin();
            const auto size   = fieldStorage_m->size();
            const auto length = referencePath_m->length();
            const auto closed = config_m.closedRing;
            for (size_t c = 0; c < particles_m.size(); ++c) {
                if (!context.trackingActive()[c]) continue;
                auto& pc             = *particles_m[c];
                const auto projected = projectedPositions_m[c];
                const auto invalid   = pc.InvalidMask.getView();
                size_t outside       = 0;
                Kokkos::parallel_reduce(
                        "FFT2D5::checkMigrationDomain", pc.getLocalNum(),
                        KOKKOS_LAMBDA(size_t i, size_t& count) {
                            if (!invalid(i)
                                && (projected(i)[0] < origin[0]
                                    || projected(i)[0] >= origin[0] + size[0]
                                    || projected(i)[1] < origin[1]
                                    || projected(i)[1] >= origin[1] + size[1]))
                                ++count;
                        },
                        outside);
                ippl::Comm->allreduce(outside, 1, std::plus<size_t>());
                if (outside)
                    throw OpalException(
                            "FFT2D5Algorithm::migrateToSlabs",
                            "Particle outside transverse FFT2D5 domain; enlarge "
                            "PIPESIZEX/PIPESIZEY.");
                auto saved = pc.migrationPositions().getView();
                auto r     = pc.R.getView();
                Kokkos::parallel_for(
                        "FFT2D5::saveCartesianPositions", pc.getLocalNum(),
                        KOKKOS_LAMBDA(size_t i) {
                            saved(i) = r(i);
                            auto q   = invalid(i) ? Vector3D_t(0.) : projected(i);
                            q[2] = q[2] >= length ? (closed ? 0. : length * (1. - 1.e-15)) : q[2];
                            r(i) = q;
                        });
                Kokkos::fence();
                pc.updateLayout(fieldStorage_m->layout(), fieldStorage_m->mesh());
                const auto restore = [&] {
                    const auto positions = pc.R.getView();
                    const auto cartesian = pc.migrationPositions().getView();
                    Kokkos::parallel_for(
                            "FFT2D5::restoreCartesianPositions", pc.getLocalNum(),
                            KOKKOS_LAMBDA(size_t i) { positions(i) = cartesian(i); });
                    Kokkos::fence();
                    pc.markMomentsDirty();
                };
                // Logical coordinates are already bounded. In particular, IPPL's
                // periodic particle BC maps the exact lower endpoint to the upper
                // endpoint; do not let it change the selected slab ownership.
                const auto particleBC = pc.getPL().getParticleBC();
                pc.setParticleBC(ippl::BC::NO);
                try {
                    pc.update();
                } catch (...) {
                    pc.setParticleBC(particleBC);
                    restore();
                    throw;
                }
                pc.setParticleBC(particleBC);
                restore();
            }
            preparePathGeometryAndProjections(context, false);
        });
    }

    void FFT2D5Algorithm::accumulateDensityHalos() {
        opalx::spacecharge::fft2d5_detail::timed("2D5 Rho Halo", [&] {
            auto& rho    = fieldStorage_m->chargeDensity();
            auto view    = rho.getView();
            const int nx = fieldStorage_m->meshSize()[0], ny = fieldStorage_m->meshSize()[1];
            // OPEN transverse boundaries: do not wrap transverse CIC spill when the
            // longitudinal layout is periodic. Only z ghosts enter halo accumulation.
            Kokkos::parallel_for(
                    "FFT2D5::clearTransverseSourceHalo",
                    Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                            {0, 0, 0}, {view.extent(0), view.extent(1), view.extent(2)}),
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                        if (i == 0 || j == 0 || i > size_t(nx) || j > size_t(ny))
                            view(i, j, k) = 0.;
                    });
            Kokkos::fence();
            rho.accumulateHalo();
        });
    }

    void FFT2D5Algorithm::fillFieldHalos() {
        opalx::spacecharge::fft2d5_detail::timed("2D5 E Halo", [&] {
            auto& field = fieldStorage_m->electricField();
            field.fillHalo();
            auto view    = field.getView();
            const int nx = fieldStorage_m->meshSize()[0], ny = fieldStorage_m->meshSize()[1];
            Kokkos::parallel_for(
                    "FFT2D5::clearTransverseFieldHalo",
                    Kokkos::MDRangePolicy<Kokkos::Rank<3>>(
                            {0, 0, 0}, {view.extent(0), view.extent(1), view.extent(2)}),
                    KOKKOS_LAMBDA(size_t i, size_t j, size_t k) {
                        if (i == 0 || j == 0 || i > size_t(nx) || j > size_t(ny))
                            view(i, j, k) = Vector3D_t(0.);
                    });
            Kokkos::fence();
        });
    }

    void FFT2D5Algorithm::preparePathGeometryAndProjections(
            const SpaceChargeSolveContext& context, bool rebuildGeometry) {
        opalx::spacecharge::fft2d5_detail::timed("2D5 Geometry", [&] {
            const auto ref = solveReferencePath_m;
            if (ref.extent(0) < 2) return;
            const auto segments = ref.extent(0) - 1;
            const T periodicLength =
                    !rebuildGeometry && config_m.closedRing ? referencePath_m->length() : 0.;
            if (segmentGeometry_m.extent(0) != segments) {
                segmentGeometry_m = SegmentGeometryView_t("FFT2D5SegmentGeometry", segments);
            }
            size_t leaves = 1;
            while (leaves < segments) {
                if (leaves > std::numeric_limits<size_t>::max() / 4)
                    throw OpalException(
                            "FFT2D5Algorithm::preparePathGeometryAndProjections",
                            "Reference path is too large for the bounds hierarchy.");
                leaves *= 2;
            }
            if (segmentBounds_m.extent(0) != 2 * leaves)
                segmentBounds_m = SegmentGeometryView_t("FFT2D5SegmentBounds", 2 * leaves);
            const auto bounds   = segmentBounds_m;
            const auto geometry = segmentGeometry_m;
            // Migration changes particle ownership, not the path geometry.
            if (rebuildGeometry) {
                // Independent segment arithmetic and endpoint boxes span the device.
                Kokkos::parallel_for(
                        "FFT2D5::preparePathSegments", leaves, KOKKOS_LAMBDA(const size_t i) {
                            if (i < segments) {
                                const Vector3D_t direction = ref(i + 1) - ref(i);
                                const T length             = euclidean_norm(direction);
                                for (unsigned d = 0; d < 3; ++d)
                                    geometry(i, d) = direction[d];
                                geometry(i, 3) = direction.dot(direction);
                                geometry(i, 4) = length;
                            }
                            for (unsigned d = 0; d < 3; ++d) {
                                bounds(leaves + i, d) =
                                        i < segments ? Kokkos::min(ref(i)[d], ref(i + 1)[d])
                                                     : Kokkos::Experimental::finite_max_v<T>;
                                bounds(leaves + i, d + 3) =
                                        i < segments ? Kokkos::max(ref(i)[d], ref(i + 1)[d])
                                                     : -Kokkos::Experimental::finite_max_v<T>;
                            }
                        });
                // The preceding kernel completes before this team reads its output.
                // Keep the prefix ordered; a parallel scan would reassociate additions.
                Kokkos::parallel_for(
                        "FFT2D5::preparePathHierarchy", Kokkos::TeamPolicy<>(1, Kokkos::AUTO),
                        KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                            Kokkos::single(Kokkos::PerTeam(team), [&] {
                                T prefix = 0.0;
                                for (size_t i = 0; i < segments; ++i) {
                                    geometry(i, 5) = prefix;
                                    prefix += geometry(i, 4);
                                }
                            });
                            for (size_t first = leaves / 2; first > 0; first /= 2) {
                                Kokkos::parallel_for(
                                        Kokkos::TeamThreadRange(team, first), [&](const size_t j) {
                                            const size_t i = first + j;
                                            for (unsigned d = 0; d < 3; ++d) {
                                                bounds(i, d) = Kokkos::min(
                                                        bounds(2 * i, d), bounds(2 * i + 1, d));
                                                bounds(i, d + 3) = Kokkos::max(
                                                        bounds(2 * i, d + 3),
                                                        bounds(2 * i + 1, d + 3));
                                            }
                                        });
                                team.team_barrier();
                            }
                        });
                Kokkos::fence();
            }
            projectedPositions_m.resize(particles_m.size());
            projectedSegments_m.resize(particles_m.size());
            for (size_t container = 0; container < particles_m.size(); ++container) {
                if (context.trackingActive()[container] == 0) continue;
                const auto count = particles_m[container]->getLocalNum();
                if (projectedPositions_m[container].extent(0) != count) {
                    projectedPositions_m[container] =
                            VectorView_t("FFT2D5ProjectedPositions", count);
                    projectedSegments_m[container] =
                            SegmentIndexView_t("FFT2D5ProjectedSegments", count);
                    Kokkos::deep_copy(projectedSegments_m[container], segments);
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
                            // The seed is only an upper bound, even after reordering/path changes.
                            const auto evaluate = [&](size_t i) {
                                const T length2 = geometry(i, 3);
                                if (length2 <= 0.0) return;
                                const Vector3D_t direction(
                                        geometry(i, 0), geometry(i, 1), geometry(i, 2));
                                const Vector3D_t start = ref(i);
                                const T fraction       = Kokkos::clamp(
                                        direction.dot(position - start) / length2, 0.0, 1.0);
                                const Vector3D_t point      = start + fraction * direction;
                                const Vector3D_t difference = position - point;
                                const T distance2           = difference.dot(difference);
                                if (distance2 < bestDistance2
                                    || (distance2 == bestDistance2 && i < best)) {
                                    best          = i;
                                    bestFraction  = fraction;
                                    closest       = point;
                                    bestDistance2 = distance2;
                                }
                            };
                            evaluate(indices(n) < segments ? indices(n) : 0);
                            // This scale bounds rounding in both the original segment arithmetic
                            // and the box-distance calculation. Borderline/overflow cases descend.
                            T scale = 0.0;
                            for (unsigned d = 0; d < 3; ++d) {
                                scale = Kokkos::max(scale, Kokkos::abs(position[d]));
                                scale = Kokkos::max(scale, Kokkos::abs(bounds(1, d)));
                                scale = Kokkos::max(scale, Kokkos::abs(bounds(1, d + 3)));
                            }
                            constexpr T epsilon = std::numeric_limits<T>::epsilon();
                            const T roundoff    = 256 * epsilon * scale * scale
                                               + 256 * std::numeric_limits<T>::min();
                            const auto lowerBound = [&](size_t node) {
                                T distance2 = 0.0;
                                for (unsigned d = 0; d < 3; ++d) {
                                    const T gap = Kokkos::max(
                                            T(0), Kokkos::max(
                                                          bounds(node, d) - position[d],
                                                          position[d] - bounds(node, d + 3)));
                                    distance2 += gap * gap;
                                }
                                return distance2;
                            };
                            size_t stack[64];
                            unsigned top  = 0;
                            stack[top++]  = 1;
                            bool fallback = false;
                            while (top != 0) {
                                const size_t node = stack[--top];
                                if (lowerBound(node) > bestDistance2 + roundoff) continue;
                                if (node >= leaves) {
                                    const size_t i = node - leaves;
                                    if (i < segments) evaluate(i);
                                } else {
                                    if (top + 2 > 64) {
                                        fallback = true;
                                        break;
                                    }
                                    const size_t left = 2 * node, right = left + 1;
                                    // Visit nearer child first to tighten the upper bound quickly.
                                    const bool leftFirst = lowerBound(left) <= lowerBound(right);
                                    stack[top++]         = leftFirst ? right : left;
                                    stack[top++]         = leftFirst ? left : right;
                                }
                            }
                            if (fallback) {
                                for (size_t i = 0; i < segments; ++i)
                                    evaluate(i);
                            }
                            const Vector3D_t tangent =
                                    Vector3D_t(
                                            geometry(best, 0), geometry(best, 1), geometry(best, 2))
                                    / geometry(best, 4);
                            const Vector3D_t normal(0.0, 1.0, 0.0);
                            const Vector3D_t binormal = ippl::cross(normal, tangent);
                            const Vector3D_t offset   = position - closest;
                            projected(n)              = Vector3D_t(
                                    offset.dot(binormal), offset.dot(normal),
                                    geometry(best, 5) + bestFraction * geometry(best, 4));
                            // Ownership wraps the exact end point to the first slab.
                            // Keep the cached interpolation coordinate consistent with it.
                            if (periodicLength > 0. && projected(n)[2] >= periodicLength)
                                projected(n)[2] = 0.;
                            indices(n) = best;
                        });
            }
            Kokkos::fence();
        });
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

        fft2d5_detail::timed("2D5 Initialize", [&] {
            // Publish dependent storage together only after path loading and allocation succeed.
            referencePath_m = std::make_unique<ReferencePath>(
                    ReferencePath::load(config_m.referencePathFile));
            solveReferencePath_m = ReferenceView_t(
                    "FFT2D5SolveReferencePath", referencePath_m->deviceView().extent(0));
            Kokkos::deep_copy(solveReferencePath_m, referencePath_m->deviceView());
            fieldStorage_m =
                    std::make_unique<FFT2D5FieldStorage>(config_m, referencePath_m->length());
            lineDensity_m = LineDensityView_t(
                    "FFT2D5LineDensity", fieldStorage_m->meshSize()[2] + LineDensityGhostCells);
            lineDensityGradient_m =
                    LineDensityView_t("FFT2D5LineDensityGradient", fieldStorage_m->meshSize()[2]);
        });
        fft2d5_detail::timed("2D5 Warmup", [&] {
            fieldStorage_m->warmup();
        });
    }

}  // namespace opalx::spacecharge
