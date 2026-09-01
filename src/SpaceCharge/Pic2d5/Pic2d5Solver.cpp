/**
 * @file Pic2d5Solver.cpp
 * @brief Implements the independent reference-path FFT2D5 self-field algorithm.
 */

#include "SpaceCharge/Pic2d5/Pic2d5Solver.h"

#include "Physics/Physics.h"
#include "SpaceCharge/SelfFieldDiagnostics.h"
#include "Utilities/OpalException.h"

#include <Kokkos_NumericTraits.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <numeric>
#include <typeindex>
#include <utility>

namespace opalx::spacecharge {

    Pic2d5Solver::Pic2d5Solver(
            Pic2d5Config config,
            std::span<const std::shared_ptr<ParticleContainer>> particleContainers)
        : config_m(std::move(config)) {
        particles_m.reserve(particleContainers.size());
        for (const auto& particles : particleContainers) {
            if (particles == nullptr) {
                throw OpalException(
                        "Pic2d5Solver::Pic2d5Solver",
                        "FFT2D5 cannot borrow a null particle container.");
            }
            particles_m.push_back(particles.get());
        }
        if (particles_m.empty()) {
            throw OpalException(
                    "Pic2d5Solver::Pic2d5Solver",
                    "FFT2D5 requires at least one particle container.");
        }
    }

    SolverCapabilities Pic2d5Solver::capabilities() const {
        SolverCapabilities result;
        result.algorithm                  = SelfFieldAlgorithmKind::Pic2d5;
        result.particleSelection          = ParticleSelectionPolicy::AllTrackingActive;
        result.supportsBinning            = false;
        result.supportsImageCharge        = false;
        result.supportsShiftedGreen       = false;
        result.supportsRedistribution     = false;
        result.supportsMultipleContainers = true;
        result.supportsPotentialOutput    = false;
        result.requiredReadableAttributes =
                ParticleAttribute::Position | ParticleAttribute::Momentum
                | ParticleAttribute::Charge | ParticleAttribute::TimeStep
                | ParticleAttribute::InvalidMask;
        result.requiredWritableAttributes = ParticleAttribute::ElectricField
                                            | ParticleAttribute::MagneticField
                                            | ParticleAttribute::TimeStep;
        return result;
    }

    void Pic2d5Solver::execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) {
        if (context.stepState().communicator.size != 1) {
            throw OpalException(
                    "Pic2d5Solver::execute",
                    "FFT2D5 currently supports only one MPI rank. Distributed fields and ORB "
                    "load balancing are not implemented for this solver.");
        }
        validateParticleMembership(context);
        ensureInitialized();
        framePolicy_m.validate(context);

        Pic2d5FramePolicy::State frameState;
        try {
            framePolicy_m.enter(context, *particles_m.front(), frameState);
            framePolicy_m.markComputedFields(frameState);
            run(context, diagnostics);
            framePolicy_m.leave(context, *particles_m.front(), frameState);
        } catch (...) {
            const std::exception_ptr original = std::current_exception();
            framePolicy_m.restore(context, *particles_m.front(), frameState);
            std::rethrow_exception(original);
        }
    }

    std::size_t Pic2d5Solver::sliceCount() const {
        return workspace_m == nullptr ? 0 : workspace_m->slices().size();
    }

    const ReferencePath::View& Pic2d5Solver::referencePathView() const {
        if (referencePath_m == nullptr) {
            throw OpalException(
                    "Pic2d5Solver::referencePathView", "FFT2D5 has not been initialized.");
        }
        return referencePath_m->deviceView();
    }

    void Pic2d5Solver::ensureInitialized() {
        if (workspace_m != nullptr) {
            return;
        }

        auto referencePath =
                std::make_unique<ReferencePath>(ReferencePath::load(config_m.referencePathFile()));
        auto workspace = std::make_unique<SliceWorkspace>(config_m, referencePath->length());
        LineDensityView lineDensity(
                "Pic2d5LineDensity", workspace->slices().size() + LineDensityGhostCells);
        LineDensityView lineDensityGradient(
                "Pic2d5LineDensityGradient", workspace->slices().size());

        referencePath_m       = std::move(referencePath);
        workspace_m           = std::move(workspace);
        lineDensity_m         = std::move(lineDensity);
        lineDensityGradient_m = std::move(lineDensityGradient);
    }

    void Pic2d5Solver::validateParticleMembership(const SolveContext& context) const {
        const auto views = context.particles().containers();
        if (views.size() != particles_m.size() || context.particles().primaryIndex() != 0) {
            throw OpalException(
                    "Pic2d5Solver::execute",
                    "SolveContext particle membership does not match FFT2D5 construction.");
        }
        using Position = typename ParticleContainer::particle_position_type;
        for (std::size_t index = 0; index < particles_m.size(); ++index) {
            const auto* handle = views[index].find(ParticleAttribute::Position);
            if (handle == nullptr || handle->nativeType() != std::type_index(typeid(Position))
                || std::addressof(handle->native<Position>())
                           != std::addressof(particles_m[index]->R)) {
                throw OpalException(
                        "Pic2d5Solver::execute",
                        "SolveContext particle order does not match FFT2D5 construction.");
            }
        }
    }

    void Pic2d5Solver::run(SolveContext& context, SelfFieldDiagnostics& diagnostics) {
        scatterToGrid(context);
        solvePoissons(diagnostics);
        calculateLineDensity();
        gatherFromGrid(context);
    }

    bool Pic2d5Solver::selected(const SolveContext& context, std::size_t index) const {
        return context.particles().containers()[index].selectedForSolve();
    }

    void Pic2d5Solver::scatterToGrid(const SolveContext& context) {
        using Policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        using Policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;

        const auto& reference        = referencePath_m->deviceView();
        const Vector3 spacing        = workspace_m->spacing();
        const Vector3 inverseSpacing = 1.0 / spacing;
        auto& rho                    = workspace_m->chargeDensity();
        const int ghostCells         = rho.getNghost();
        const auto localDomain       = rho.getLayout().getLocalNDIndex();
        const auto rhoView           = rho.getView();
        const Vector3 origin         = workspace_m->origin();
        Kokkos::deep_copy(rhoView, 0.0);

        for (std::size_t containerIndex = 0; containerIndex < particles_m.size();
             ++containerIndex) {
            if (!selected(context, containerIndex)) {
                continue;
            }
            ParticleContainer& particles = *particles_m[containerIndex];
            particles.updateMoments();
            particles.scaleDtByCharge();
            try {
                const auto position       = particles.R.getView();
                const auto momentum       = particles.P.getView();
                const double meanPs       = particles.getMeanP()[2];
                const auto timeStepCharge = particles.dt.getView();
                const auto invalid        = particles.InvalidMask.getView();
                if (config_m.scatterLongitudinally()) {
                    Kokkos::parallel_for(
                            "Pic2d5Solver::scatterToGrid::3d", particles.getLocalNum(),
                            KOKKOS_LAMBDA(const std::size_t index) {
                                scatterParticle<true>(
                                        index, position, momentum, reference, meanPs,
                                        timeStepCharge, invalid, inverseSpacing, ghostCells,
                                        localDomain, rhoView, origin);
                            });
                } else {
                    Kokkos::parallel_for(
                            "Pic2d5Solver::scatterToGrid::2d", particles.getLocalNum(),
                            KOKKOS_LAMBDA(const std::size_t index) {
                                scatterParticle<false>(
                                        index, position, momentum, reference, meanPs,
                                        timeStepCharge, invalid, inverseSpacing, ghostCells,
                                        localDomain, rhoView, origin);
                            });
                }
                Kokkos::fence();
                particles.unscaleDtByCharge();
            } catch (...) {
                particles.unscaleDtByCharge();
                throw;
            }
        }

        if (config_m.closedRing()) {
            constexpr std::size_t firstReal = LineDensityFirstRealCell;
            const std::size_t lastReal = workspace_m->meshSize()[2] - 1 + LineDensityFirstRealCell;
            Kokkos::parallel_for(
                    "Pic2d5Solver::scatterToGrid::closed-boundary",
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
        const double normalization = context.stepState().timeStep * cellVolume;
        if (normalization == 0.0) {
            throw OpalException(
                    "Pic2d5Solver::scatterToGrid",
                    "FFT2D5 charge-density normalization requires a non-zero time step.");
        }
        Kokkos::parallel_for(
                "Pic2d5Solver::scatterToGrid::normalize",
                Policy3({0, 0, 0}, {rhoView.extent(0), rhoView.extent(1), rhoView.extent(2)}),
                KOKKOS_LAMBDA(const std::size_t i, const std::size_t j, const std::size_t k) {
                    rhoView(i, j, k) /= normalization;
                });
        Kokkos::fence();
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void Pic2d5Solver::scatterParticle(
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

    KOKKOS_FUNCTION void Pic2d5Solver::convertToFrenet(
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

    KOKKOS_FUNCTION void Pic2d5Solver::boostMomentum(
            double meanLongitudinalMomentum, Vector3& momentum) {
        const double beamGamma =
                Kokkos::sqrt(1.0 + meanLongitudinalMomentum * meanLongitudinalMomentum);
        const double particleGamma = Kokkos::sqrt(1.0 + momentum[2] * momentum[2]);
        momentum[2] = beamGamma * momentum[2] - meanLongitudinalMomentum * particleGamma;
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void Pic2d5Solver::deposit(
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

    void Pic2d5Solver::solvePoissons(SelfFieldDiagnostics& diagnostics) {
        using Policy         = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto& rho3d          = workspace_m->chargeDensity();
        auto& electric3d     = workspace_m->electricField();
        auto electric3dView  = electric3d.getView();
        const int ghostCells = electric3d.getNghost();

        for (std::size_t z = 0; z < workspace_m->slices().size(); ++z) {
            auto& slice    = workspace_m->slices()[z];
            auto rho2dView = slice.chargeDensity->getView();
            Kokkos::deep_copy(
                    rho2dView,
                    Kokkos::subview(rho3d.getView(), Kokkos::ALL(), Kokkos::ALL(), z + ghostCells));
            Kokkos::parallel_for(
                    "Pic2d5Solver::solvePoissons::coupling",
                    Policy({0, 0}, {rho2dView.extent(0), rho2dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        rho2dView(i, j) /= Physics::epsilon_0;
                    });
            {
                auto backendEvent =
                        diagnostics.scopedEvent(SelfFieldEventKind::BackendSolve, "FFT2D5-slice");
                slice.solver->solve();
                diagnostics.completeBackendSolve();
            }
            Kokkos::fence();
            auto electric2dView = slice.electricField->getView();
            Kokkos::parallel_for(
                    "Pic2d5Solver::solvePoissons::copy",
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
            Kokkos::parallel_for(
                    "Pic2d5Solver::solvePoissons::closed-boundary",
                    Policy({0, 0}, {electric3dView.extent(0), electric3dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        electric3dView(i, j, leftGhost)  = electric3dView(i, j, rightGhost - 1);
                        electric3dView(i, j, rightGhost) = electric3dView(i, j, leftGhost + 1);
                    });
        } else {
            Kokkos::parallel_for(
                    "Pic2d5Solver::solvePoissons::open-boundary",
                    Policy({0, 0}, {electric3dView.extent(0), electric3dView.extent(1)}),
                    KOKKOS_LAMBDA(const std::size_t i, const std::size_t j) {
                        electric3dView(i, j, leftGhost)  = Vector3(0.0);
                        electric3dView(i, j, rightGhost) = Vector3(0.0);
                    });
        }
        Kokkos::fence();
    }

    void Pic2d5Solver::calculateLineDensity() {
        if (config_m.longitudinalFieldMode() == Pic2d5LongitudinalFieldMode::None) {
            Kokkos::deep_copy(lineDensity_m, 0.0);
            Kokkos::deep_copy(lineDensityGradient_m, 0.0);
            return;
        }

        using Policy                 = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        const auto rho               = workspace_m->chargeDensity().getView();
        auto hostLineDensity         = Kokkos::create_mirror_view(lineDensity_m);
        const std::size_t sliceCount = workspace_m->slices().size();
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

        const Vector3 spacing       = workspace_m->spacing();
        const double transverseArea = spacing[0] * spacing[1];
        const auto lineDensity      = lineDensity_m;
        Kokkos::parallel_for(
                "Pic2d5Solver::calculateLineDensity::convert", sliceCount + LineDensityGhostCells,
                KOKKOS_LAMBDA(const std::size_t k) { lineDensity(k) *= transverseArea; });

        const auto gradient = lineDensityGradient_m;
        const double dz     = spacing[2];
        Kokkos::parallel_for(
                "Pic2d5Solver::calculateLineDensity::gradient", sliceCount,
                KOKKOS_LAMBDA(const std::size_t k) {
                    gradient(k) = (lineDensity(k + LineDensityFirstRealCell + 1)
                                   - lineDensity(k + LineDensityFirstRealCell - 1))
                                  / (2.0 * dz);
                });
        Kokkos::fence();
    }

    double Pic2d5Solver::longitudinalGeometryFactor() const {
        const double pipeRadius = std::min(workspace_m->size()[0], workspace_m->size()[1]);
        double factor           = OpenG0;
        if (config_m.longitudinalFieldMode() == Pic2d5LongitudinalFieldMode::Cylindrical) {
            factor = CircularPipeG0 + 2.0 * std::log(pipeRadius / config_m.beamRadius());
        } else if (config_m.longitudinalFieldMode() == Pic2d5LongitudinalFieldMode::Plates) {
            factor = ParallelPlatesG0
                     + 2.0 * std::log(4.0 * pipeRadius / (Physics::pi * config_m.beamRadius()));
        }
        return factor / (4.0 * Physics::pi * Physics::epsilon_0);
    }

    void Pic2d5Solver::gatherFromGrid(const SolveContext& context) {
        const auto reference         = referencePath_m->deviceView();
        const Vector3 inverseSpacing = 1.0 / workspace_m->spacing();
        const int ghostCells         = workspace_m->chargeDensity().getNghost();
        const auto localDomain       = workspace_m->layout().getLocalNDIndex();
        const Vector3 origin         = workspace_m->origin();
        const auto electricGrid      = workspace_m->electricField().getView();
        const auto gradient          = lineDensityGradient_m;
        const double geometryFactor  = longitudinalGeometryFactor();

        for (std::size_t containerIndex = 0; containerIndex < particles_m.size();
             ++containerIndex) {
            if (!selected(context, containerIndex)) {
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
                    "Pic2d5Solver::gatherFromGrid", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const std::size_t index) {
                        gatherParticle(
                                index, position, momentum, reference, beamGamma, beamBeta, electric,
                                magnetic, invalid, inverseSpacing, ghostCells, localDomain,
                                electricGrid, origin, geometryFactor, gradient);
                    });
            Kokkos::fence();
        }
    }

    KOKKOS_FUNCTION void Pic2d5Solver::gatherParticle(
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

    KOKKOS_FUNCTION void Pic2d5Solver::gatherTransverse(
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

    KOKKOS_FUNCTION void Pic2d5Solver::unboostFields(
            std::size_t index, double beamGamma, double beamBeta, const VectorView& electric,
            const VectorView& magnetic) {
        electric(index)[0] *= beamGamma;
        electric(index)[1] *= beamGamma;
        magnetic(index)[0] = beamBeta * electric(index)[1] / Physics::c;
        magnetic(index)[1] = -beamBeta * electric(index)[0] / Physics::c;
        magnetic(index)[2] = 0.0;
    }

    KOKKOS_FUNCTION void Pic2d5Solver::convertFieldsFromFrenet(
            std::size_t index, const Vector3& binormal, const Vector3& normal,
            const Vector3& tangent, const VectorView& electric, const VectorView& magnetic) {
        electric(index) = electric(index)[0] * binormal + electric(index)[1] * normal
                          + electric(index)[2] * tangent;
        magnetic(index) = magnetic(index)[0] * binormal + magnetic(index)[1] * normal
                          + magnetic(index)[2] * tangent;
    }

    template <typename View>
    KOKKOS_FUNCTION bool Pic2d5Solver::makeWeights(
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

    KOKKOS_FUNCTION Pic2d5Solver::Vector3 Pic2d5Solver::gatherSlice(
            VectorGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights, int x,
            int y, int z) {
        Vector3 result = lowerWeights[0] * lowerWeights[1] * field(x - 1, y - 1, z);
        result += upperWeights[0] * lowerWeights[1] * field(x, y - 1, z);
        result += lowerWeights[0] * upperWeights[1] * field(x - 1, y, z);
        result += upperWeights[0] * upperWeights[1] * field(x, y, z);
        return result;
    }

    KOKKOS_FUNCTION void Pic2d5Solver::deposit2d(
            ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights, int x,
            int y, int z, double charge) {
        Kokkos::atomic_add(&field(x - 1, y - 1, z), lowerWeights[0] * lowerWeights[1] * charge);
        Kokkos::atomic_add(&field(x, y - 1, z), upperWeights[0] * lowerWeights[1] * charge);
        Kokkos::atomic_add(&field(x - 1, y, z), lowerWeights[0] * upperWeights[1] * charge);
        Kokkos::atomic_add(&field(x, y, z), upperWeights[0] * upperWeights[1] * charge);
    }

    KOKKOS_FUNCTION void Pic2d5Solver::deposit3d(
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
