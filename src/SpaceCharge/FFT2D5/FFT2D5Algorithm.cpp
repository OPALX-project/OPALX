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

            enterParticleSolveFrames(context, labToSolve);
            try {
                doRunSolver(context);
            } catch (...) {
                restoreParticleSolveFrames(context, labToSolve);
                throw;
            }
            restoreParticleSolveFrames(context, labToSolve);
            result.backendSolves = fieldStorage_m->backendSolves();
            if (ippl::Comm->size() > 1) {
                const auto local = result.backendSolves;
                ippl::Comm->allreduce(local, result.backendSolves, 1, std::plus<size_t>());
            }
        });
        firstPhysicalSolve_m = false;
        return result;
    }

    void FFT2D5Algorithm::enterParticleSolveFrames(
            const SpaceChargeSolveContext& context, const CoordinateSystemTrafo& labToSolve) {
        for (std::size_t index = 0; index < particles_m.size(); ++index) {
            if (context.trackingActive()[index] != 0) {
                auto& particles           = *particles_m[index];
                const auto trackerToSolve = labToSolve * particles.getToLabTrafo();
                clearSelfFields(particles);
                enterSolveFrame({trackerToSolve, trackerToSolve.inverted()}, particles);
            }
        }
    }

    void FFT2D5Algorithm::restoreParticleSolveFrames(
            const SpaceChargeSolveContext& context, const CoordinateSystemTrafo& labToSolve) {
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
    }

    void FFT2D5Algorithm::migrateToSlabs(const SpaceChargeSolveContext& context) {
        fft2d5_detail::timed("2D5 Migration", [&] {
            for (size_t c = 0; c < particles_m.size(); ++c)
                if (context.trackingActive()[c]) migrateContainer(c);
            preparePathGeometryAndProjections(context, false);
        });
    }

    void FFT2D5Algorithm::checkMigrationDomain(size_t container) {
        const auto origin = fieldStorage_m->origin(), size = fieldStorage_m->size();
        const auto projected = projectedPositions_m[container];
        const auto invalid   = particles_m[container]->InvalidMask.getView();
        auto& pc             = *particles_m[container];
        size_t outside       = 0;
        Kokkos::parallel_reduce(
                "FFT2D5::checkMigrationDomain", pc.getLocalNum(),
                KOKKOS_LAMBDA(size_t i, size_t& count) {
                    if (!invalid(i)
                        && (projected(i)[0] < origin[0] || projected(i)[0] >= origin[0] + size[0]
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
    }

    void FFT2D5Algorithm::saveMigrationPositions(size_t container) {
        auto& pc             = *particles_m[container];
        const auto projected = projectedPositions_m[container];
        const auto invalid   = pc.InvalidMask.getView();
        const auto length    = referencePath_m->length();
        const auto closed    = config_m.closedRing;
        auto saved           = pc.migrationPositions().getView();
        auto r               = pc.R.getView();
        Kokkos::parallel_for(
                "FFT2D5::saveCartesianPositions", pc.getLocalNum(), KOKKOS_LAMBDA(size_t i) {
                    saved(i) = r(i);
                    auto q   = invalid(i) ? Vector3D_t(0.) : projected(i);
                    q[2]     = q[2] >= length ? (closed ? 0. : length * (1. - 1.e-15)) : q[2];
                    r(i)     = q;
                });
        Kokkos::fence();
    }

    void FFT2D5Algorithm::restoreMigrationPositions(ParticleContainer& pc) {
        const auto positions = pc.R.getView();
        const auto cartesian = pc.migrationPositions().getView();
        Kokkos::parallel_for(
                "FFT2D5::restoreCartesianPositions", pc.getLocalNum(),
                KOKKOS_LAMBDA(size_t i) { positions(i) = cartesian(i); });
        Kokkos::fence();
        pc.markMomentsDirty();
    }

    void FFT2D5Algorithm::migrateContainer(size_t container) {
        checkMigrationDomain(container);
        auto& pc = *particles_m[container];
        saveMigrationPositions(container);
        pc.updateLayout(fieldStorage_m->layout(), fieldStorage_m->mesh());
        // Disable spatial periodic remapping of the already-bounded logical coordinates.
        const auto particleBC = pc.getPL().getParticleBC();
        pc.setParticleBC(ippl::BC::NO);
        try {
            pc.update();
        } catch (...) {
            pc.setParticleBC(particleBC);
            restoreMigrationPositions(pc);
            throw;
        }
        pc.setParticleBC(particleBC);
        restoreMigrationPositions(pc);
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
        fft2d5_detail::timed("2D5 Geometry", [&] {
            if (solveReferencePath_m.extent(0) < 2) return;
            const size_t leaves = resizePathGeometry();
            if (rebuildGeometry) {
                buildPathSegments(leaves);
                buildPathHierarchy(leaves);
                Kokkos::fence();
            }
            projectedPositions_m.resize(particles_m.size());
            projectedSegments_m.resize(particles_m.size());
            const T periodicLength =
                    !rebuildGeometry && config_m.closedRing ? referencePath_m->length() : 0.;
            for (size_t c = 0; c < particles_m.size(); ++c)
                if (context.trackingActive()[c]) projectContainer(c, leaves, periodicLength);
            Kokkos::fence();
        });
    }

    size_t FFT2D5Algorithm::resizePathGeometry() {
        const size_t segments = solveReferencePath_m.extent(0) - 1;
        if (segmentGeometry_m.extent(0) != segments)
            segmentGeometry_m = SegmentGeometryView_t("FFT2D5SegmentGeometry", segments);
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
        return leaves;
    }

    void FFT2D5Algorithm::buildPathSegments(size_t leaves) {
        const auto ref      = solveReferencePath_m;
        const auto geometry = segmentGeometry_m, bounds = segmentBounds_m;
        const size_t segments = ref.extent(0) - 1;
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
                        bounds(leaves + i, d) = i < segments
                                                        ? Kokkos::min(ref(i)[d], ref(i + 1)[d])
                                                        : Kokkos::Experimental::finite_max_v<T>;
                        bounds(leaves + i, d + 3) =
                                i < segments ? Kokkos::max(ref(i)[d], ref(i + 1)[d])
                                             : -Kokkos::Experimental::finite_max_v<T>;
                    }
                });
    }

    void FFT2D5Algorithm::buildPathHierarchy(size_t leaves) {
        const auto geometry = segmentGeometry_m, bounds = segmentBounds_m;
        const size_t segments = solveReferencePath_m.extent(0) - 1;
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
                                        bounds(i, d) =
                                                Kokkos::min(bounds(2 * i, d), bounds(2 * i + 1, d));
                                        bounds(i, d + 3) = Kokkos::max(
                                                bounds(2 * i, d + 3), bounds(2 * i + 1, d + 3));
                                    }
                                });
                        team.team_barrier();
                    }
                });
    }

    void FFT2D5Algorithm::projectContainer(size_t container, size_t leaves, T periodicLength) {
        const size_t count    = particles_m[container]->getLocalNum();
        const auto ref        = solveReferencePath_m;
        const size_t segments = ref.extent(0) - 1;
        if (projectedPositions_m[container].extent(0) != count) {
            projectedPositions_m[container] = VectorView_t("FFT2D5ProjectedPositions", count);
            projectedSegments_m[container]  = SegmentIndexView_t("FFT2D5ProjectedSegments", count);
            Kokkos::deep_copy(projectedSegments_m[container], segments);
        }
        const auto positions = particles_m[container]->R.getView();
        const auto invalid   = particles_m[container]->InvalidMask.getView();
        const auto projected = projectedPositions_m[container];
        const auto indices   = projectedSegments_m[container];
        const auto geometry = segmentGeometry_m, bounds = segmentBounds_m;
        Kokkos::parallel_for(
                "FFT2D5::prepareParticleProjections", count, KOKKOS_LAMBDA(const size_t n) {
                    projectPreparedParticle(
                            n, positions, invalid, ref, geometry, bounds, projected, indices,
                            leaves, periodicLength);
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

    void FFT2D5Algorithm::validateDepositionContainer(
            size_t container, bool usePreparedProjections) {
        auto* pc            = particles_m[container];
        const auto ref      = solveReferencePath_m;
        const auto geometry = usePreparedProjections ? segmentGeometry_m : SegmentGeometryView_t();
        const auto projectedR =
                usePreparedProjections ? projectedPositions_m[container] : VectorView_t();
        const auto projectedSegment =
                usePreparedProjections ? projectedSegments_m[container] : SegmentIndexView_t();
        pc->updateMoments();
        if (pc->getTotalNum() > 0 && pc->getChargePerParticle() == 0.0) {
            throw OpalException(
                    "FFT2D5Algorithm::scatterToGrid",
                    "Per-particle charge is zero for an active FFT2D5 container.");
        }
        // A clamped projection would put particles beyond a finite path
        // into its end slice. Fail before mutating the per-particle time steps.
        if (!config_m.closedRing) {
            const auto length    = referencePath_m->length();
            const auto tolerance = 64 * std::numeric_limits<T>::epsilon() * std::max(T(1), length);
            size_t outside       = 0;
            Kokkos::parallel_reduce(
                    "FFT2D5::checkPathCoverage", pc->getLocalNum(),
                    fft2d5_detail::CheckPathCoverageFunctor{
                            pc->R.getView(), pc->P.getView(), ref, geometry, projectedR,
                            projectedSegment, pc->InvalidMask.getView(), length, tolerance},
                    outside);
            ippl::Comm->allreduce(outside, 1, std::plus<size_t>());
            if (outside != 0) {
                throw OpalException(
                                "FFT2D5Algorithm::scatterToGrid",
                                "Open reference path does not cover " + std::to_string(outside)
                                + " active particle(s). Extend REFPATHFNAME to cover the full "
                                  "bunch throughout tracking; endpoint clamping would pile up charge.");
            }
        }
    }

    void FFT2D5Algorithm::sumLineDensityProfile() {
        const auto rho3d             = fieldStorage_m->chargeDensity().getView();
        const auto deviceLineDensity = lineDensity_m;
        const auto hostLineDensity   = Kokkos::create_mirror_view(lineDensity_m);
        const auto numSlices         = fieldStorage_m->meshSize()[2];
        const auto localDomain       = fieldStorage_m->layout().getLocalNDIndex();
        const bool distributed       = ippl::Comm->size() > 1;
        const size_t begin           = distributed ? 0 : LineDensityFirstRealCell;
        const size_t end             = size_t(localDomain[2].length()) + (distributed ? 2 : 1);
        // Reuse the profile as device scratch until the final global profile is
        // assembled. One batched reduction and one small copy replace a synchronized
        // scalar reduction/copy for every slice. Always allocate an independent host
        // mirror: create_mirror_view may alias on HostSpace, and clearing the global
        // profile below must not erase the local slice sums.
        const auto hostSliceSums = Kokkos::create_mirror(deviceLineDensity);
#ifdef KOKKOS_ENABLE_CUDA
        constexpr int vectorLength = 32;
#else
        constexpr int vectorLength = 1;
#endif
        Kokkos::parallel_for(
                "FFT2D5::batchedSliceDensity",
                Kokkos::TeamPolicy<>(end - begin, Kokkos::AUTO, vectorLength),
                fft2d5_detail::SliceDensitySumFunctor{rho3d, begin, deviceLineDensity});
        Kokkos::deep_copy(hostSliceSums, deviceLineDensity);
        Kokkos::deep_copy(hostLineDensity, 0.);
        for (size_t k = begin; k < end; ++k) {
            const T sum = hostSliceSums(k);
            int global  = int(k) + localDomain[2].first() - LineDensityFirstRealCell;
            if (config_m.closedRing) global = (global + int(numSlices)) % int(numSlices);
            if (global >= 0 && global < int(numSlices))
                hostLineDensity(global + LineDensityFirstRealCell) += sum;
        }
        // Only the small 1D profile is replicated; no 3D field all-reduction.
        if (ippl::Comm->size() > 1)
            ippl::Comm->allreduce(
                    hostLineDensity.data(), numSlices + LineDensityGhostCells, std::plus<T>());
        // Set the ghost cells to the boundary conditions
        if (config_m.closedRing) {
            hostLineDensity(0) = hostLineDensity(numSlices + LineDensityGhostCells - 2);
            hostLineDensity(numSlices + LineDensityGhostCells - 1) =
                    hostLineDensity(LineDensityFirstRealCell);
        } else {
            hostLineDensity(0)                                     = 0.0;
            hostLineDensity(numSlices + LineDensityGhostCells - 1) = 0.0;
        }
        Kokkos::deep_copy(deviceLineDensity, hostLineDensity);
    }

    void FFT2D5Algorithm::setTransverseFieldGhosts() {
        using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto e3d     = fieldStorage_m->electricField().getView();
        // Set the ghost slices
        constexpr auto leftGhost = 0;
        const auto rightGhost    = e3d.extent(2) - 1;
        if (ippl::Comm->size() > 1) {
            fillFieldHalos();
        } else if (config_m.closedRing) {
            Kokkos::parallel_for(
                    "Solve2d5::solvePoissons::ghostClosed",
                    Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                    fft2d5_detail::FieldGhostsFunctor<true>{e3d, leftGhost, rightGhost});
            Kokkos::fence();
        } else {
            Kokkos::parallel_for(
                    "Solve2d5::solvePoissons::ghostOpen",
                    Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                    fft2d5_detail::FieldGhostsFunctor<false>{e3d, leftGhost, rightGhost});
            Kokkos::fence();
        }
    }

    FFT2D5Algorithm::T FFT2D5Algorithm::longitudinalCoupling() const {
        const auto pipeRadius = std::min(fieldStorage_m->size()[0], fieldStorage_m->size()[1]);
        T g;
        if (config_m.longitudinalFieldMode == LongitudinalFieldMode::Cylindrical) {
            g = CircularPipeG0 + 2 * Kokkos::log(pipeRadius / config_m.beamRadius);
        } else if (config_m.longitudinalFieldMode == LongitudinalFieldMode::Plates) {
            g = ParallelPlatesG0
                + 2 * Kokkos::log(4 * pipeRadius / Physics::pi / config_m.beamRadius);
        } else {
            g = OpenG0;
        }
        return g / (4 * Physics::pi * Physics::epsilon_0);
    }

}  // namespace opalx::spacecharge
