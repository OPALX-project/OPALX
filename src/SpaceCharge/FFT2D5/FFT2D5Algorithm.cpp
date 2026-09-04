/**
 * @file FFT2D5Algorithm.cpp
 * @brief Implements the independent reference-path FFT2D5 space-charge algorithm.
 */

#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.h"

#include "Physics/Physics.h"
#include "SpaceCharge/SpaceChargeDiagnostics.h"
#include "Utilities/OpalException.h"

#include <Kokkos_NumericTraits.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <utility>

namespace opalx::spacecharge {

    FFT2D5Algorithm::FFT2D5Algorithm(
            FFT2D5Config config, std::span<const ParticleFieldBinding3D> particleBindings)
        : config_m(std::move(config)) {
        particles_m.reserve(particleBindings.size());
        for (const ParticleFieldBinding3D& binding : particleBindings) {
            if (binding.container == nullptr) {
                throw OpalException(
                        "FFT2D5Algorithm::FFT2D5Algorithm",
                        "FFT2D5 cannot borrow a null particle container.");
            }
            particles_m.push_back(binding.container);
        }
        if (particles_m.empty()) {
            throw OpalException(
                    "FFT2D5Algorithm::FFT2D5Algorithm",
                    "FFT2D5 requires at least one particle container.");
        }
    }

    SpaceChargeCapabilities FFT2D5Algorithm::capabilities() const {
        SpaceChargeCapabilities result;
        result.particleSelection       = ParticleSelectionMode::AllTrackingActive;
        result.supportsBinning         = false;
        result.supportsImageCharge     = false;
        result.supportsShiftedGreen    = false;
        result.supportsRedistribution  = false;
        result.supportsPotentialOutput = false;
        return result;
    }

    void FFT2D5Algorithm::solve(
            SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) {
        if (context.stepState().mpiSize != 1) {
            throw OpalException(
                    "FFT2D5Algorithm::solve",
                    "FFT2D5 currently supports only one MPI rank. Distributed fields and ORB "
                    "load balancing are not implemented for this solver.");
        }
        ensureInitialized();

        SpaceChargeFrameTransform<double, 3> frameTransform(
                context.stepState().frames, *particles_m.front());
        frameTransform.enter();
        frameTransform.markComputedFields();
        solveFields(context, diagnostics);
        frameTransform.leave();
    }

    std::size_t FFT2D5Algorithm::sliceCount() const {
        return fieldStorage_m == nullptr ? 0 : fieldStorage_m->slices().size();
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

        // Construct every dependent object locally first. A failed path load or field allocation
        // leaves the algorithm uninitialized, so the next solve can retry without partial state.
        auto referencePath =
                std::make_unique<ReferencePath>(ReferencePath::load(config_m.referencePathFile()));
        auto fieldStorage = std::make_unique<FFT2D5FieldStorage>(config_m, referencePath->length());
        LineDensityView lineDensity(
                "FFT2D5LineDensity", fieldStorage->slices().size() + LineDensityGhostCells);
        LineDensityView lineDensityGradient(
                "FFT2D5LineDensityGradient", fieldStorage->slices().size());

        referencePath_m       = std::move(referencePath);
        fieldStorage_m        = std::move(fieldStorage);
        lineDensity_m         = std::move(lineDensity);
        lineDensityGradient_m = std::move(lineDensityGradient);
    }

    void FFT2D5Algorithm::solveFields(
            SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) {
        scatterToGrid(context);
        solveSlicePoissonProblems(diagnostics);
        calculateLineDensity();
        gatherFromGrid(context);
    }

    bool FFT2D5Algorithm::isSelected(
            const SpaceChargeSolveContext& context, std::size_t index) const {
        return context.particles().bindings()[index].selectedForSolve;
    }

    void FFT2D5Algorithm::scatterToGrid(const SpaceChargeSolveContext& context) {
        using Policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        using Policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

        const auto& reference        = referencePath_m->deviceView();
        const Vector3 spacing        = fieldStorage_m->spacing();
        const Vector3 inverseSpacing = 1.0 / spacing;
        auto& rho                    = fieldStorage_m->chargeDensity();
        const int ghostCells         = rho.getNghost();
        const auto localDomain       = rho.getLayout().getLocalNDIndex();
        const auto rhoView           = rho.getView();
        const Vector3 origin         = fieldStorage_m->origin();
        Kokkos::deep_copy(rhoView, 0.0);

        // FFT2D5 is the multi-container algorithm: all tracking-active bindings selected by
        // SpaceChargeSolver contribute to the same slice field.
        for (std::size_t containerIndex = 0; containerIndex < particles_m.size();
             ++containerIndex) {
            if (!isSelected(context, containerIndex)) {
                continue;
            }
            ParticleContainer& particles = *particles_m[containerIndex];
            particles.updateMoments();
            particles.scaleDtByCharge();
            const auto position       = particles.R.getView();
            const auto momentum       = particles.P.getView();
            const double meanPs       = particles.getMeanP()[2];
            const auto timeStepCharge = particles.dt.getView();
            const auto invalid        = particles.InvalidMask.getView();
            if (config_m.scatterLongitudinally()) {
                Kokkos::parallel_for(
                        "FFT2D5Algorithm::scatterToGrid::3d", particles.getLocalNum(),
                        KOKKOS_LAMBDA(const std::size_t index) {
                            scatterParticle<true>(
                                    index, position, momentum, reference, meanPs, timeStepCharge,
                                    invalid, inverseSpacing, ghostCells, localDomain, rhoView,
                                    origin);
                        });
            } else {
                Kokkos::parallel_for(
                        "FFT2D5Algorithm::scatterToGrid::2d", particles.getLocalNum(),
                        KOKKOS_LAMBDA(const std::size_t index) {
                            scatterParticle<false>(
                                    index, position, momentum, reference, meanPs, timeStepCharge,
                                    invalid, inverseSpacing, ghostCells, localDomain, rhoView,
                                    origin);
                        });
            }
            Kokkos::fence();
            particles.unscaleDtByCharge();
        }

        if (config_m.closedRing()) {
            // Fold periodic ghost charge into the opposite real boundary before slice solves.
            constexpr std::size_t firstReal = LineDensityFirstRealCell;
            const std::size_t lastReal =
                    fieldStorage_m->meshSize()[2] - 1 + LineDensityFirstRealCell;
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::scatterToGrid::closed-boundary",
                    Policy2({0, 0}, {rhoView.extent(0), rhoView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        Kokkos::atomic_add(&rhoView(i, j, firstReal), rhoView(i, j, lastReal + 1));
                        Kokkos::atomic_add(&rhoView(i, j, lastReal), rhoView(i, j, firstReal - 1));
                        rhoView(i, j, lastReal + 1)  = 0.0;
                        rhoView(i, j, firstReal - 1) = 0.0;
                    });
            Kokkos::fence();
        }

        const double cellVolume =
                std::reduce(spacing.begin(), spacing.end(), 1.0, std::multiplies());

        // Deposition used dt*Q, so divide by dt and cell volume to retain volume charge density.
        // The 2D slice solver receives this density unchanged apart from 1/epsilon_0 below.
        const double normalization = context.stepState().timeStep * cellVolume;
        if (normalization == 0.0) {
            throw OpalException(
                    "FFT2D5Algorithm::scatterToGrid",
                    "FFT2D5 charge-density normalization requires a non-zero time step.");
        }
        Kokkos::parallel_for(
                "FFT2D5Algorithm::scatterToGrid::normalize",
                Policy3({0, 0, 0}, {rhoView.extent(0), rhoView.extent(1), rhoView.extent(2)}),
                KOKKOS_LAMBDA(const std::size_t i, const std::size_t j, const std::size_t k) {
                    rhoView(i, j, k) /= normalization;
                });
        Kokkos::fence();
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void FFT2D5Algorithm::scatterParticle(
            std::size_t index, const VectorView& position, const VectorView& momentum,
            const ReferenceView& reference, double meanLongitudinalMomentum,
            const ScalarView& timeStepCharge, const BooleanView& invalid, Vector3 inverseSpacing,
            int ghostCells, ippl::NDIndex<3> localDomain, ScalarGridView3 chargeDensity,
            Vector3 origin) {
        if (invalid(index)) {
            return;
        }
        Vector3 frenetPosition, frenetMomentum, binormal, normal, tangent;
        convertToFrenet(
                index, position, momentum, reference, frenetPosition, frenetMomentum, binormal,
                normal, tangent);
        boostMomentum(meanLongitudinalMomentum, frenetMomentum);
        deposit<ScatterLongitudinally>(
                index, frenetPosition, timeStepCharge, inverseSpacing, ghostCells, localDomain,
                chargeDensity, origin);
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::convertToFrenet(
            std::size_t index, const VectorView& position, const VectorView& momentum,
            const ReferenceView& reference, Vector3& frenetPosition, Vector3& frenetMomentum,
            Vector3& binormal, Vector3& normal, Vector3& tangent) {
        double bestDistanceSquared = Kokkos::Experimental::finite_max_v<double>;
        std::size_t bestSegment    = 0;
        double bestFraction        = 0.0;
        Vector3 bestDirection(0.0);
        Vector3 bestPoint(0.0);
        const Vector3 particlePosition = position(index);
        for (std::size_t segment = 0; segment + 1 < reference.extent(0); ++segment) {
            const Vector3 direction       = reference(segment + 1) - reference(segment);
            const double directionSquared = direction.dot(direction);
            if (directionSquared <= 0.0) {
                continue;
            }
            const double fraction = Kokkos::clamp(
                    direction.dot(particlePosition - reference(segment)) / directionSquared, 0.0,
                    1.0);
            const Vector3 pathPoint      = reference(segment) + fraction * direction;
            const Vector3 difference     = particlePosition - pathPoint;
            const double distanceSquared = difference.dot(difference);
            if (distanceSquared < bestDistanceSquared) {
                bestSegment         = segment;
                bestFraction        = fraction;
                bestDirection       = direction;
                bestDistanceSquared = distanceSquared;
                bestPoint           = pathPoint;
            }
        }

        const double segmentLength = bestDirection.Pnorm();
        tangent                    = bestDirection / segmentLength;
        normal                     = Vector3(0.0);
        normal[1]                  = 1.0;
        binormal                   = ippl::cross(normal, tangent);

        frenetPosition[2] = bestFraction * segmentLength;
        for (std::size_t segment = 1; segment <= bestSegment; ++segment) {
            const Vector3 pathSegment = reference(segment) - reference(segment - 1);
            frenetPosition[2] += pathSegment.Pnorm();
        }
        const Vector3 difference = particlePosition - bestPoint;
        frenetPosition[0]        = difference.dot(binormal);
        frenetPosition[1]        = difference.dot(normal);
        frenetMomentum[0]        = momentum(index).dot(binormal);
        frenetMomentum[1]        = momentum(index).dot(normal);
        frenetMomentum[2]        = momentum(index).dot(tangent);
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::boostMomentum(
            double meanLongitudinalMomentum, Vector3& momentum) {
        const double beamGamma =
                Kokkos::sqrt(1.0 + meanLongitudinalMomentum * meanLongitudinalMomentum);
        const double particleGamma = Kokkos::sqrt(1.0 + momentum[2] * momentum[2]);
        momentum[2] = beamGamma * momentum[2] - meanLongitudinalMomentum * particleGamma;
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void FFT2D5Algorithm::deposit(
            std::size_t index, Vector3 frenetPosition, const ScalarView& timeStepCharge,
            Vector3 inverseSpacing, int ghostCells, const ippl::NDIndex<3>& localDomain,
            ScalarGridView3 chargeDensity, Vector3 origin) {
        Vector3 upperWeights, lowerWeights;
        ippl::Vector<int, 3> indices;
        if (!makeWeights(
                    frenetPosition, origin, inverseSpacing, ghostCells, localDomain, chargeDensity,
                    upperWeights, lowerWeights, indices)) {
            return;
        }
        if constexpr (ScatterLongitudinally) {
            deposit3d(
                    chargeDensity, lowerWeights, upperWeights, indices[0], indices[1], indices[2],
                    timeStepCharge(index));
        } else {
            deposit2d(
                    chargeDensity, lowerWeights, upperWeights, indices[0], indices[1], indices[2],
                    timeStepCharge(index));
        }
    }

    void FFT2D5Algorithm::solveSlicePoissonProblems(SpaceChargeDiagnostics& diagnostics) {
        using Policy         = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto& rho3d          = fieldStorage_m->chargeDensity();
        auto& electric3d     = fieldStorage_m->electricField();
        auto electric3dView  = electric3d.getView();
        const int ghostCells = electric3d.getNghost();

        for (std::size_t z = 0; z < fieldStorage_m->slices().size(); ++z) {
            auto& slice    = fieldStorage_m->slices()[z];
            auto rho2dView = slice.chargeDensity->getView();

            // IPPL fields do not accept a Kokkos subview as owned storage. Copy each z plane into
            // its persistent 2D field, solve it, and copy the transverse result back to staging.
            Kokkos::deep_copy(
                    rho2dView,
                    Kokkos::subview(rho3d.getView(), Kokkos::ALL(), Kokkos::ALL(), z + ghostCells));
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::solveSlicePoissonProblems::coupling",
                    Policy({0, 0}, {rho2dView.extent(0), rho2dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        rho2dView(i, j) /= Physics::epsilon_0;
                    });
            slice.solver->solve();
            diagnostics.recordBackendSolve();
            Kokkos::fence();
            auto electric2dView = slice.electricField->getView();
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::solveSlicePoissonProblems::copy",
                    Policy({0, 0}, {electric2dView.extent(0), electric2dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        electric3dView(i, j, z + ghostCells)[0] = electric2dView(i, j)[0];
                        electric3dView(i, j, z + ghostCells)[1] = electric2dView(i, j)[1];
                        electric3dView(i, j, z + ghostCells)[2] = 0.0;
                    });
            Kokkos::fence();
        }

        constexpr std::size_t leftGhost = 0;
        const std::size_t rightGhost    = electric3dView.extent(2) - 1;
        if (config_m.closedRing()) {
            // Periodic ghosts expose the opposite real slice to the bilateral gather.
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::solveSlicePoissonProblems::closed-boundary",
                    Policy({0, 0}, {electric3dView.extent(0), electric3dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        electric3dView(i, j, leftGhost)  = electric3dView(i, j, rightGhost - 1);
                        electric3dView(i, j, rightGhost) = electric3dView(i, j, leftGhost + 1);
                    });
        } else {
            // Open longitudinal boundaries contribute no transverse field outside the path.
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::solveSlicePoissonProblems::open-boundary",
                    Policy({0, 0}, {electric3dView.extent(0), electric3dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        electric3dView(i, j, leftGhost)  = Vector3(0.0);
                        electric3dView(i, j, rightGhost) = Vector3(0.0);
                    });
        }
        Kokkos::fence();
    }

    void FFT2D5Algorithm::calculateLineDensity() {
        if (config_m.longitudinalFieldMode() == FFT2D5LongitudinalFieldMode::None) {
            Kokkos::deep_copy(lineDensity_m, 0.0);
            Kokkos::deep_copy(lineDensityGradient_m, 0.0);
            return;
        }

        using Policy                 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        const auto rho               = fieldStorage_m->chargeDensity().getView();
        auto hostLineDensity         = Kokkos::create_mirror_view(lineDensity_m);
        const std::size_t sliceCount = fieldStorage_m->slices().size();
        for (std::size_t k = 0; k < sliceCount + LineDensityGhostCells; ++k) {
            double sum = 0.0;
            Kokkos::parallel_reduce(
                    Policy({0, 0}, {rho.extent(0), rho.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j, double& local) {
                        local += rho(i, j, k);
                    },
                    sum);
            hostLineDensity(k) = sum;
        }
        if (config_m.closedRing()) {
            hostLineDensity(0)              = hostLineDensity(sliceCount);
            hostLineDensity(sliceCount + 1) = hostLineDensity(LineDensityFirstRealCell);
        } else {
            hostLineDensity(0)              = 0.0;
            hostLineDensity(sliceCount + 1) = 0.0;
        }
        Kokkos::deep_copy(lineDensity_m, hostLineDensity);

        const Vector3 spacing       = fieldStorage_m->spacing();
        const double transverseArea = spacing[0] * spacing[1];
        const auto lineDensity      = lineDensity_m;
        Kokkos::parallel_for(
                "FFT2D5Algorithm::calculateLineDensity::convert",
                sliceCount + LineDensityGhostCells,
                KOKKOS_LAMBDA(const std::size_t k) { lineDensity(k) *= transverseArea; });

        const auto gradient = lineDensityGradient_m;
        const double dz     = spacing[2];

        // The two explicit ghost entries make this central difference valid for every real slice.
        // With one slice, both neighbors are equal for a ring or zero for open boundaries, giving
        // the required zero longitudinal gradient without an out-of-range special case.
        Kokkos::parallel_for(
                "FFT2D5Algorithm::calculateLineDensity::gradient", sliceCount,
                KOKKOS_LAMBDA(const std::size_t k) {
                    gradient(k) = (lineDensity(k + LineDensityFirstRealCell + 1)
                                   - lineDensity(k + LineDensityFirstRealCell - 1))
                                  / (2.0 * dz);
                });
        Kokkos::fence();
    }

    double FFT2D5Algorithm::longitudinalGeometryFactor() const {
        const double pipeRadius = std::min(fieldStorage_m->size()[0], fieldStorage_m->size()[1]);
        double factor           = OpenG0;
        if (config_m.longitudinalFieldMode() == FFT2D5LongitudinalFieldMode::Cylindrical) {
            factor = CircularPipeG0 + 2.0 * std::log(pipeRadius / config_m.beamRadius());
        } else if (config_m.longitudinalFieldMode() == FFT2D5LongitudinalFieldMode::Plates) {
            factor = ParallelPlatesG0
                     + 2.0 * std::log(4.0 * pipeRadius / (Physics::pi * config_m.beamRadius()));
        }
        return factor / (4.0 * Physics::pi * Physics::epsilon_0);
    }

    void FFT2D5Algorithm::gatherFromGrid(const SpaceChargeSolveContext& context) {
        const auto reference         = referencePath_m->deviceView();
        const Vector3 inverseSpacing = 1.0 / fieldStorage_m->spacing();
        const int ghostCells         = fieldStorage_m->chargeDensity().getNghost();
        const auto localDomain       = fieldStorage_m->layout().getLocalNDIndex();
        const Vector3 origin         = fieldStorage_m->origin();
        const auto electricGrid      = fieldStorage_m->electricField().getView();
        const auto gradient          = lineDensityGradient_m;
        const double geometryFactor  = longitudinalGeometryFactor();

        for (std::size_t containerIndex = 0; containerIndex < particles_m.size();
             ++containerIndex) {
            if (!isSelected(context, containerIndex)) {
                continue;
            }
            ParticleContainer& particles = *particles_m[containerIndex];
            particles.updateMoments();
            const auto position    = particles.R.getView();
            const auto momentum    = particles.P.getView();
            const double meanPs    = particles.getMeanP()[2];
            const double beamGamma = Kokkos::sqrt(1.0 + meanPs * meanPs);
            const double beamBeta  = meanPs / beamGamma;
            const auto electric    = particles.E.getView();
            const auto magnetic    = particles.B.getView();
            const auto invalid     = particles.InvalidMask.getView();
            Kokkos::parallel_for(
                    "FFT2D5Algorithm::gatherFromGrid", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const std::size_t index) {
                        gatherParticle(
                                index, position, momentum, reference, beamGamma, beamBeta, electric,
                                magnetic, invalid, inverseSpacing, ghostCells, localDomain,
                                electricGrid, origin, geometryFactor, gradient);
                    });
            Kokkos::fence();
        }
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::gatherParticle(
            std::size_t index, const VectorView& position, const VectorView& momentum,
            const ReferenceView& reference, double beamGamma, double beamBeta,
            const VectorView& electric, const VectorView& magnetic, const BooleanView& invalid,
            Vector3 inverseSpacing, int ghostCells, ippl::NDIndex<3> localDomain,
            VectorGridView3 electricGrid, Vector3 origin, double geometryFactor,
            LineDensityView lineDensityGradient) {
        if (invalid(index)) {
            return;
        }
        Vector3 frenetPosition, frenetMomentum, binormal, normal, tangent;
        convertToFrenet(
                index, position, momentum, reference, frenetPosition, frenetMomentum, binormal,
                normal, tangent);
        gatherTransverse(
                index, frenetPosition, electric, inverseSpacing, ghostCells, localDomain,
                electricGrid, origin);
        unboostFields(index, beamGamma, beamBeta, electric, magnetic);

        const int longitudinalIndex =
                static_cast<int>((frenetPosition[2] - origin[2]) * inverseSpacing[2] + 0.5);
        double gradient = 0.0;
        if (longitudinalIndex >= 0
            && longitudinalIndex < static_cast<int>(lineDensityGradient.extent(0))) {
            gradient = lineDensityGradient(longitudinalIndex);
        }
        electric(index)[2] = -geometryFactor * gradient / (beamGamma * beamGamma);
        convertFieldsFromFrenet(index, binormal, normal, tangent, electric, magnetic);
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::gatherTransverse(
            std::size_t index, Vector3 frenetPosition, const VectorView& electric,
            Vector3 inverseSpacing, int ghostCells, const ippl::NDIndex<3>& localDomain,
            VectorGridView3 electricGrid, Vector3 origin) {
        Vector3 upperWeights, lowerWeights;
        ippl::Vector<int, 3> indices;
        if (makeWeights(
                    frenetPosition, origin, inverseSpacing, ghostCells, localDomain, electricGrid,
                    upperWeights, lowerWeights, indices)) {
            // Preserve master's bilateral transverse gather across both adjacent longitudinal
            // slices. Longitudinal weighting is intentionally not introduced here.
            electric(index) = gatherSlice(
                                      electricGrid, lowerWeights, upperWeights, indices[0],
                                      indices[1], indices[2] - 1)
                              + gatherSlice(
                                      electricGrid, lowerWeights, upperWeights, indices[0],
                                      indices[1], indices[2]);
        }
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::unboostFields(
            std::size_t index, double beamGamma, double beamBeta, const VectorView& electric,
            const VectorView& magnetic) {
        electric(index)[0] *= beamGamma;
        electric(index)[1] *= beamGamma;
        magnetic(index)[0] = beamBeta * electric(index)[1] / Physics::c;
        magnetic(index)[1] = -beamBeta * electric(index)[0] / Physics::c;
        magnetic(index)[2] = 0.0;
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::convertFieldsFromFrenet(
            std::size_t index, const Vector3& binormal, const Vector3& normal,
            const Vector3& tangent, const VectorView& electric, const VectorView& magnetic) {
        electric(index) = electric(index)[0] * binormal + electric(index)[1] * normal
                          + electric(index)[2] * tangent;
        magnetic(index) = magnetic(index)[0] * binormal + magnetic(index)[1] * normal
                          + magnetic(index)[2] * tangent;
    }

    template <typename View>
    KOKKOS_FUNCTION bool FFT2D5Algorithm::makeWeights(
            Vector3 frenetPosition, Vector3 origin, Vector3 inverseSpacing, int ghostCells,
            const ippl::NDIndex<3>& localDomain, const View& view, Vector3& upperWeights,
            Vector3& lowerWeights, ippl::Vector<int, 3>& indices) {
        const Vector3 location = (frenetPosition - origin) * inverseSpacing + 0.5;
        indices                = location;
        upperWeights           = location - indices;
        lowerWeights           = 1.0 - upperWeights;
        indices                = indices - localDomain.first() + ghostCells;
        bool inBounds          = indices[0] > 0 && indices[0] < static_cast<int>(view.extent(0));
        inBounds = inBounds && indices[1] > 0 && indices[1] < static_cast<int>(view.extent(1));
        inBounds = inBounds && indices[2] > 0 && indices[2] < static_cast<int>(view.extent(2));
        return inBounds;
    }

    KOKKOS_FUNCTION FFT2D5Algorithm::Vector3 FFT2D5Algorithm::gatherSlice(
            VectorGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights, int x,
            int y, int z) {
        Vector3 result = lowerWeights[0] * lowerWeights[1] * field(x - 1, y - 1, z);
        result += upperWeights[0] * lowerWeights[1] * field(x, y - 1, z);
        result += lowerWeights[0] * upperWeights[1] * field(x - 1, y, z);
        result += upperWeights[0] * upperWeights[1] * field(x, y, z);
        return result;
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::deposit2d(
            ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights, int x,
            int y, int z, double charge) {
        Kokkos::atomic_add(&field(x - 1, y - 1, z), lowerWeights[0] * lowerWeights[1] * charge);
        Kokkos::atomic_add(&field(x, y - 1, z), upperWeights[0] * lowerWeights[1] * charge);
        Kokkos::atomic_add(&field(x - 1, y, z), lowerWeights[0] * upperWeights[1] * charge);
        Kokkos::atomic_add(&field(x, y, z), upperWeights[0] * upperWeights[1] * charge);
    }

    KOKKOS_FUNCTION void FFT2D5Algorithm::deposit3d(
            ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights, int x,
            int y, int z, double charge) {
        Kokkos::atomic_add(
                &field(x - 1, y - 1, z - 1),
                lowerWeights[0] * lowerWeights[1] * lowerWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x, y - 1, z - 1),
                upperWeights[0] * lowerWeights[1] * lowerWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x - 1, y, z - 1),
                lowerWeights[0] * upperWeights[1] * lowerWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x, y, z - 1), upperWeights[0] * upperWeights[1] * lowerWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x - 1, y - 1, z),
                lowerWeights[0] * lowerWeights[1] * upperWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x, y - 1, z), upperWeights[0] * lowerWeights[1] * upperWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x - 1, y, z), lowerWeights[0] * upperWeights[1] * upperWeights[2] * charge);
        Kokkos::atomic_add(
                &field(x, y, z), upperWeights[0] * upperWeights[1] * upperWeights[2] * charge);
    }

}  // namespace opalx::spacecharge
