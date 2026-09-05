//
// Copyright (c) 2008 - 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
//
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
/**
 * @file FFT2D5Algorithm.hpp
 * @brief Numerical routines retained from master a5c859240, PartBunch/Solve2d5.hpp.
 *
 * Arithmetic, interpolation, field signs, stage ordering, and diagnostic hooks are retained.
 * Integration changes replace PartBunch access with borrowed containers and per-call state,
 * redirect workspace access, select tracking-active containers, and check deposition inputs.
 */
#ifndef OPALX_FFT2D5_ALGORITHM_HPP
#define OPALX_FFT2D5_ALGORITHM_HPP

#include "Physics/Physics.h"
#include "VectorMath.h"
#include "Utilities/OpalException.h"
#include <Kokkos_NumericTraits.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

namespace opalx::spacecharge {

template <typename DiagnosticPolicy>
void FFT2D5Algorithm::doRunSolver(const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic) {
    scatterToGrid<DiagnosticPolicy>(context, diagnostic);
    solvePoissons<DiagnosticPolicy>(diagnostic);
    calculateLineDensity<DiagnosticPolicy>(diagnostic);
    gatherFromGrid<DiagnosticPolicy>(context, diagnostic);
}

template <typename DiagnosticPolicy>
std::unique_ptr<DiagnosticPolicy> FFT2D5Algorithm::createDiagnostic(NullDiagnostic::Kind kind) {
    auto result = std::make_unique<DiagnosticPolicy>(kind);
    result->initialise(particles_m, fieldStorage_m->chargeDensity(), lineDensity_m,
                       lineDensityGradient_m, fieldStorage_m->electricField());
    return result;
}

// Place the charge from all the particles into the charge density grid.
// Multiple particle containers in the bunch are supported.
// The particle container array dt is temporarily used to contain the charge to deposit.
// The real work is handled by the doScatterToGrid Kokkos kernel.
// The periodic boundary conditions are handled when the ring is marked as closed.
template <typename DiagnosticPolicy>
void FFT2D5Algorithm::scatterToGrid(const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic) {
    using Policy2D_t = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    using Policy3D_t = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
    if (referencePath_m->deviceView().extent(0) > 1) {
        const auto& ref    = referencePath_m->deviceView();
        const auto invDr   = 1.0 / fieldStorage_m->spacing();
        const int nghost   = fieldStorage_m->chargeDensity().getNghost();
        const auto& layout = fieldStorage_m->chargeDensity().getLayout();
        const auto& lDom   = layout.getLocalNDIndex();
        const auto rhoView = fieldStorage_m->chargeDensity().getView();
        const auto& origin = fieldStorage_m->origin();
        // Do the scattering for all the particle containers
        Kokkos::deep_copy(rhoView, 0.0);
        for (std::size_t container = 0; container < particles_m.size(); ++container) {
            if (context.trackingActive()[container] == 0) {
                continue;
            }
            auto* pc = particles_m[container];
            pc->updateMoments();
            if (pc->getTotalNum() > 0 && pc->getChargePerParticle() == 0.0) {
                throw OpalException("FFT2D5Algorithm::scatterToGrid",
                                    "Per-particle charge is zero for an active FFT2D5 container.");
            }
            pc->scaleDtByCharge();
            const auto& r       = pc->R.getView();
            const auto& p       = pc->P.getView();
            const auto meanPs   = pc->getMeanP().data_m[2];
            const auto& dt      = pc->dt.getView();
            const auto& invalid = pc->InvalidMask.getView();
            if (config_m.scatterLongitudinally) {
                Kokkos::parallel_for(
                        "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                        KOKKOS_LAMBDA(const size_t n) {
                            doScatterToGrid<true>(
                                    n, r, p, ref, meanPs, dt, invalid, invDr, nghost, lDom, rhoView,
                                    origin, diagnostic);
                        });
            } else {
                Kokkos::parallel_for(
                        "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                        KOKKOS_LAMBDA(const size_t n) {
                            doScatterToGrid<false>(
                                    n, r, p, ref, meanPs, dt, invalid, invDr, nghost, lDom, rhoView,
                                    origin, diagnostic);
                        });
            }
            Kokkos::fence();
            diagnostic.scatterCharge(rhoView);
            pc->unscaleDtByCharge();
        }
        // Handle the closed ring periodic boundary condition
        if (config_m.closedRing) {
            constexpr auto firstRealZ = LineDensityFirstRealCell;
            const auto lastRealZ      = fieldStorage_m->meshSize()[2] - 1 + LineDensityFirstRealCell;
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid::boundaries",
                    Policy2D_t({0, 0}, {rhoView.extent(0), rhoView.extent(1)}),
                    KOKKOS_LAMBDA(const size_t i, const size_t j) {
                        rhoView(i, j, firstRealZ) += rhoView(i, j, lastRealZ + 1);
                        rhoView(i, j, lastRealZ) += rhoView(i, j, firstRealZ - 1);
                        rhoView(i, j, lastRealZ + 1)  = 0;
                        rhoView(i, j, firstRealZ - 1) = 0;
                    });
            Kokkos::fence();
        }
        // Now scale by volume and time step to get the proper charge density.
        const auto cellVolume = std::reduce(fieldStorage_m->spacing().begin(), fieldStorage_m->spacing().end(), 1.0, std::multiplies());
        const auto scale      = context.stepState().timeStep * cellVolume;
        if (!std::isfinite(scale) || scale == 0.0) {
            throw OpalException("FFT2D5Algorithm::scatterToGrid",
                                "FFT2D5 deposition requires a finite nonzero time step and cell volume.");
        }
        Kokkos::parallel_for(
                "Solve2d5::scatterToGrid::scale",
                Policy3D_t({0, 0, 0}, {rhoView.extent(0), rhoView.extent(1), rhoView.extent(2)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j, const size_t k) {
                    rhoView(i, j, k) /= scale;
                });
        Kokkos::fence();
        diagnostic.scatterChargeDensity(rhoView);
    }
}

// This is a Kokkos kernel that deposits the charge for a single particle into the
// charge density grid.  The particle's coordinates are first translated into
// the Frenet-Serret coordinate system using the previously loaded reference path.
// The particle is then boosted into the beam's reference frame before the charge
// is scattered using one of two CiC algorithms depending on the ScatterLongitudinally
// template parameter.   When true, trilinear CiC is used.  When false, the particle
// is assigned to a single z slice and bilinear CiC is used in that slice.
template <bool ScatterLongitudinally, typename DiagnosticPolicy>
KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::doScatterToGrid(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        const T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
        const int nghost, const ippl::NDIndex<3U> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
        DiagnosticPolicy diagnostic) {
    if (!invalid(n)) {
        // Into Frenet-Serret coordinates
        Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
        convertToFrenetSerret(n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit);
        diagnostic.frenetSerretScatter(n, fsR, fsP, invalid(n));
        // Boost into the beam frame
        boostToBeamFrame(meanPs, fsP);
        diagnostic.boostToBeam(n, fsR, fsP, invalid(n));
        // CiC scatter the charge to the 3D rho grid
        scatterToRho<ScatterLongitudinally>(n, fsR, dt, invDr, nghost, lDom, rho, origin);
    }
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::convertToFrenetSerret(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit, Vector3D_t& tUnit) {
    // Find the segment with the shortest normal to the point
    T bestDist2 = Kokkos::Experimental::finite_max_v<T>;
    size_t bestI{};
    T bestU{};
    Vector3D_t bestDi{};
    Vector3D_t bestRc{};
    auto rn       = r(n);
    auto segments = ref.extent(0);
    for (size_t i = 0; i < segments - 1; ++i) {
        auto refi     = ref(i);
        auto refip1   = ref(i + 1);
        Vector3D_t di = refip1 - refi;
        const T di2   = di.dot(di);
        if (di2 > 0.0) {
            const T u       = Kokkos::clamp(di.dot(rn - refi) / di2, 0.0, 1.0);
            Vector3D_t rc   = refi + u * di;
            Vector3D_t diff = rn - rc;
            const T dist2   = diff.dot(diff);
            if (dist2 < bestDist2) {
                bestI     = i;
                bestU     = u;
                bestDi    = di;
                bestDist2 = dist2;
                bestRc    = rc;
            }
        }
    }
    // The basis unit vectors
    tUnit = bestDi / euclidean_norm(bestDi);
    nUnit = {0, 1, 0};
    bUnit = ippl::cross(nUnit, tUnit);
    // Calculate the S coordinate
    fsR.data_m[2] = bestU * euclidean_norm(bestDi);
    for (size_t i = 1; i <= bestI; ++i) {
        Vector3D_t diff = ref(i) - ref(i - 1);
        fsR.data_m[2] += euclidean_norm(diff);
    }
    // Remaining position coordinates
    const Vector3D_t rDiff = r(n) - bestRc;
    fsR.data_m[0]          = rDiff.dot(bUnit);
    fsR.data_m[1]          = rDiff.dot(nUnit);
    // The momentum coordinates
    fsP.data_m[0] = p(n).dot(bUnit);
    fsP.data_m[1] = p(n).dot(nUnit);
    fsP.data_m[2] = p(n).dot(tUnit);
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::boostToBeamFrame(const T meanPs, Vector3D_t& fsP) {
    // Transform the longitudinal momentum coordinate into the beam reference frame
    auto gammaB   = Kokkos::sqrt(1.0 + meanPs * meanPs);
    auto gamma    = Kokkos::sqrt(1.0 + fsP.data_m[2] * fsP.data_m[2]);
    fsP.data_m[2] = gammaB * fsP.data_m[2] - meanPs * gamma;
}

template <typename ViewType>
KOKKOS_INLINE_FUNCTION bool FFT2D5Algorithm::makeWeights(
        Vector3D_t fsR, Vector3D_t origin, Vector3D_t invDr, int nghost,
        const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
        ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args) {
    const auto l                = (fsR - origin) * invDr + 0.5;
    ippl::Vector<int, 3U> index = l;
    whi                         = l - index;
    wlo                         = 1.0 - whi;
    args                        = index - lDom.first() + nghost;
    // CIC touches args[d] and args[d] - 1, so valid args are
    // [1, extent - 1]. Anything outside would underflow or
    // overrun the field view on the device.
    bool inBounds = args[0] > 0 && args[0] < static_cast<int>(view.extent(0));
    inBounds      = inBounds && args[1] > 0 && args[1] < static_cast<int>(view.extent(1));
    inBounds      = inBounds && args[2] > 0 && args[2] < static_cast<int>(view.extent(2));
    return inBounds;
}

template <bool ScatterLongitudinally>
KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::scatterToRho(
        const size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U>& lDom, ScalarGridView3D_t rho, Vector3D_t origin) {
    // CiC scatter the charge to the 3D rho grid
    ippl::Vector<T, Dim> whi, wlo;
    ippl::Vector<int, Dim> args;
    if (makeWeights(fsR, origin, invDr, nghost, lDom, rho, whi, wlo, args)) {
        if constexpr (ScatterLongitudinally) {
            scatter3D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
        } else {
            scatter2D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
        }
    }
}

template <typename DiagnosticPolicy>
void FFT2D5Algorithm::calculateLineDensity(DiagnosticPolicy diagnostic) {
    if (config_m.longitudinalFieldMode != LongitudinalFieldMode::None) {
        using Policy2D_t       = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        const auto rho3d       = fieldStorage_m->chargeDensity().getView();
        auto deviceLineDensity = lineDensity_m;
        auto hostLineDensity   = Kokkos::create_mirror_view(lineDensity_m);
        const auto numSlices   = fieldStorage_m->meshSize()[2];
        // Calculate the total charge density for each z slice
        for (size_t k = 0; k < numSlices + LineDensityGhostCells; ++k) {
            T sum{};
            Kokkos::parallel_reduce(
                    Policy2D_t({0, 0}, {rho3d.extent(0), rho3d.extent(1)}),
                    KOKKOS_LAMBDA(const size_t i, const size_t j, T& localSum) {
                        localSum += rho3d(i, j, k);
                    },
                    sum);
            hostLineDensity(k) = sum;
        }
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
        diagnostic.totalDensity(lineDensity_m);
        // Convert this to line density
        auto dx = fieldStorage_m->spacing()[0];
        auto dy = fieldStorage_m->spacing()[1];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::convert", numSlices + LineDensityGhostCells,
                KOKKOS_LAMBDA(const size_t k) { deviceLineDensity(k) *= dx * dy; });
        Kokkos::fence();
        diagnostic.lineDensity(lineDensity_m);
        // Find the gradient
        auto lineDensityGradient = lineDensityGradient_m;
        const auto dz            = fieldStorage_m->chargeDensity().get_mesh().getMeshSpacing().data_m[2];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::gradient", numSlices,
                KOKKOS_LAMBDA(const size_t k) {
                    lineDensityGradient(k) = (deviceLineDensity(k + LineDensityFirstRealCell + 1)
                                              - deviceLineDensity(k + LineDensityFirstRealCell - 1))
                                             / (2.0 * dz);
                });
        Kokkos::fence();
        diagnostic.lineDensityGradient(lineDensityGradient_m);
    } else {
        Kokkos::deep_copy(lineDensity_m, 0);
        diagnostic.totalDensity(lineDensity_m);
        diagnostic.lineDensity(lineDensity_m);
        Kokkos::deep_copy(lineDensityGradient_m, 0);
        diagnostic.lineDensityGradient(lineDensityGradient_m);
    }
}

template <typename DiagnosticPolicy>
void FFT2D5Algorithm::solvePoissons(DiagnosticPolicy diagnostic) {
    using Policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
    auto e3d     = fieldStorage_m->electricField().getView();
    auto nghost  = fieldStorage_m->electricField().getNghost();
    // Copy the 3D charge density grid into the array of 2D grids, solve the 2D poisson on each,
    // then copy the E field results into the 3D E field grid.
    for (size_t z = 0; z < fieldStorage_m->slices().size(); ++z) {
        auto& s = fieldStorage_m->slices()[z];
        // Do the 2D solve to get Ex and Ey
        auto rho2d = s.chargeDensity->getView();
        Kokkos::deep_copy(
                rho2d, Kokkos::subview(fieldStorage_m->chargeDensity().getView(), Kokkos::ALL(), Kokkos::ALL(), z + nghost));
        // Scale by the coupling constant then solve
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::coupling",
                Policy({0, 0}, {rho2d.extent(0), rho2d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    rho2d(i, j) *= 1 / Physics::epsilon_0;
                });
        s.solver->solve();
        Kokkos::fence();
        diagnostic.potential(rho2d, z + nghost);
        auto e2d = s.electricField->getView();
        // Copy the 2D E field into the 3D E field grid
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::copy", Policy({0, 0}, {e2d.extent(0), e2d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, z + nghost)[0] = e2d(i, j)[0];
                    e3d(i, j, z + nghost)[1] = e2d(i, j)[1];
                    e3d(i, j, z + nghost)[2] = 0.0;
                });
        Kokkos::fence();
    }
    // Set the ghost slices
    constexpr auto leftGhost = 0;
    auto rightGhost          = e3d.extent(2) - 1;
    if (config_m.closedRing) {
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::ghostClosed",
                Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, leftGhost)[0]  = e3d(i, j, rightGhost - 1)[0];
                    e3d(i, j, leftGhost)[1]  = e3d(i, j, rightGhost - 1)[1];
                    e3d(i, j, leftGhost)[2]  = e3d(i, j, rightGhost - 1)[2];
                    e3d(i, j, rightGhost)[0] = e3d(i, j, leftGhost + 1)[0];
                    e3d(i, j, rightGhost)[1] = e3d(i, j, leftGhost + 1)[1];
                    e3d(i, j, rightGhost)[2] = e3d(i, j, leftGhost + 1)[2];
                });
        Kokkos::fence();
    } else {
        Kokkos::parallel_for(
                "Solve2d5::solvePoissons::ghostOpen",
                Policy({0, 0}, {e3d.extent(0), e3d.extent(1)}),
                KOKKOS_LAMBDA(const size_t i, const size_t j) {
                    e3d(i, j, leftGhost)[0]  = 0;
                    e3d(i, j, leftGhost)[1]  = 0;
                    e3d(i, j, leftGhost)[2]  = 0;
                    e3d(i, j, rightGhost)[0] = 0;
                    e3d(i, j, rightGhost)[1] = 0;
                    e3d(i, j, rightGhost)[2] = 0;
                });
        Kokkos::fence();
    }
    diagnostic.eField(e3d);
}

template <typename DiagnosticPolicy>
void FFT2D5Algorithm::gatherFromGrid(const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic) {
    if (referencePath_m->deviceView().extent(0) > 1) {
        for (std::size_t container = 0; container < particles_m.size(); ++container) {
            if (context.trackingActive()[container] == 0) {
                continue;
            }
            auto* pc = particles_m[container];
            const auto& r                  = pc->R.getView();
            const auto& p                  = pc->P.getView();
            const auto& ref                = referencePath_m->deviceView();
            const auto meanPs              = pc->getMeanP().data_m[2];
            const auto& e                  = pc->E.getView();
            const auto& b                  = pc->B.getView();
            const auto& invalid            = pc->InvalidMask.getView();
            const auto& dr                 = fieldStorage_m->spacing();
            const auto invDr               = 1.0 / dr;
            const int nghost               = fieldStorage_m->chargeDensity().getNghost();
            const auto& layout             = fieldStorage_m->chargeDensity().getLayout();
            const auto& lDom               = layout.getLocalNDIndex();
            const auto& origin             = fieldStorage_m->origin();
            const auto gammaB              = Kokkos::sqrt(1.0 + meanPs * meanPs);
            const auto betaB               = meanPs / gammaB;
            const auto lineDensityGradient = lineDensityGradient_m;
            const auto eField              = fieldStorage_m->electricField().getView();
            const auto pipeRadius          = std::min(fieldStorage_m->size()[0], fieldStorage_m->size()[1]);
            T gBy4PiEpsilon0;
            if (config_m.longitudinalFieldMode == LongitudinalFieldMode::Cylindrical) {
                gBy4PiEpsilon0 = CircularPipeG0 + 2 * Kokkos::log(pipeRadius / config_m.beamRadius);
            } else if (config_m.longitudinalFieldMode == LongitudinalFieldMode::Plates) {
                gBy4PiEpsilon0 = ParallelPlatesG0
                                 + 2 * Kokkos::log(4 * pipeRadius / Physics::pi / config_m.beamRadius);
            } else {
                gBy4PiEpsilon0 = OpenG0;
            }
            gBy4PiEpsilon0 /= 4 * Physics::pi * Physics::epsilon_0;
            Kokkos::parallel_for(
                    "Solve2d5::gatherFromGrid", pc->getLocalNum(), KOKKOS_LAMBDA(const size_t n) {
                        doGatherFromGrid(
                                n, r, p, ref, gammaB, betaB, e, b, invalid, invDr, nghost, lDom,
                                eField, origin, gBy4PiEpsilon0, lineDensityGradient, diagnostic);
                    });
            Kokkos::fence();
        }
    }
}

template <typename DiagnosticPolicy>
KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::doGatherFromGrid(
        const size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
        const T beamGamma, const T beamBeta, const VectorView_t& e, const VectorView_t& b,
        const BooleanView_t& invalid, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U> lDom, VectorGridView3D_t eField, Vector3D_t origin,
        T gBy4PiEpsilon0, LineDensityView_t lineDensityGradient, DiagnosticPolicy diagnostic) {
    if (!invalid(n)) {
        // Into Frenet-Serret coordinates
        Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
        convertToFrenetSerret(n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit);
        diagnostic.frenetSerretGather(n, fsR, fsP, invalid(n));
        // CiC Gather the boosted E field
        gatherFromEField(n, fsR, e, invDr, nghost, lDom, eField, origin);
        diagnostic.gatherEField(n, e(n), b(n), invalid(n));
        // Unboost from the beam frame
        unboostFromBeamFrame(n, beamGamma, beamBeta, e, b);
        diagnostic.deboostFromBeam(n, e(n), b(n), invalid(n));
        // Calculate the longitudinal E field in Frenet-Serret coordinates
        const int index = (fsR.data_m[2] - origin.data_m[2]) * invDr.data_m[2] + 0.5;
        T ldg{};
        if (index >= 0 && index < static_cast<int>(lineDensityGradient.extent(0))) {
            ldg = lineDensityGradient(index);
        }
        e(n).data_m[2] = -gBy4PiEpsilon0 * ldg / beamGamma / beamGamma;
        diagnostic.longitudinalField(n, e(n), b(n), invalid(n));
        // And finally back into lab coordinates
        convertFromFrenetSerret(n, bUnit, nUnit, tUnit, e, b);
        diagnostic.labFrameFields(n, e(n), b(n), invalid(n));
    }
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::gatherFromEField(
        const size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, const int nghost,
        const ippl::NDIndex<3U>& lDom, VectorGridView3D_t eField, Vector3D_t origin) {
    // CiC gather the boosted E field to the 3D E field grid
    ippl::Vector<T, Dim> whi, wlo;
    ippl::Vector<int, Dim> args;
    if (makeWeights(fsR, origin, invDr, nghost, lDom, eField, whi, wlo, args)) {
        e(n) = gather2D(eField, wlo, whi, args[0], args[1], args[2] - 1)
               + gather2D(eField, wlo, whi, args[0], args[1], args[2]);
    }
}

KOKKOS_INLINE_FUNCTION FFT2D5Algorithm::Vector3D_t FFT2D5Algorithm::gather2D(
        VectorGridView3D_t eField, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z) {
    Vector3D_t result;
    result = wlo[0] * wlo[1] * eField(x - 1, y - 1, z);
    result += whi[0] * wlo[1] * eField(x, y - 1, z);
    result += wlo[0] * whi[1] * eField(x - 1, y, z);
    result += whi[0] * whi[1] * eField(x, y, z);
    return result;
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::scatter3D(
        ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z, T charge) {
    Kokkos::atomic_add(&rho(x - 1, y - 1, z - 1), wlo[0] * wlo[1] * wlo[2] * charge);
    Kokkos::atomic_add(&rho(x, y - 1, z - 1), whi[0] * wlo[1] * wlo[2] * charge);
    Kokkos::atomic_add(&rho(x - 1, y, z - 1), wlo[0] * whi[1] * wlo[2] * charge);
    Kokkos::atomic_add(&rho(x, y, z - 1), whi[0] * whi[1] * wlo[2] * charge);
    Kokkos::atomic_add(&rho(x - 1, y - 1, z), wlo[0] * wlo[1] * whi[2] * charge);
    Kokkos::atomic_add(&rho(x, y - 1, z), whi[0] * wlo[1] * whi[2] * charge);
    Kokkos::atomic_add(&rho(x - 1, y, z), wlo[0] * whi[1] * whi[2] * charge);
    Kokkos::atomic_add(&rho(x, y, z), whi[0] * whi[1] * whi[2] * charge);
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::scatter2D(
        ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
        int x, int y, int z, T charge) {
    Kokkos::atomic_add(&rho(x - 1, y - 1, z), wlo[0] * wlo[1] * charge);
    Kokkos::atomic_add(&rho(x, y - 1, z), whi[0] * wlo[1] * charge);
    Kokkos::atomic_add(&rho(x - 1, y, z), wlo[0] * whi[1] * charge);
    Kokkos::atomic_add(&rho(x, y, z), whi[0] * whi[1] * charge);
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::unboostFromBeamFrame(
        size_t n, T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b) {
    // Transform the E field from the boosted beam frame into E and B fields
    e(n).data_m[0] *= beamGamma;
    e(n).data_m[1] *= beamGamma;
    b(n).data_m[0] = beamBeta * e(n).data_m[1] / Physics::c;
    b(n).data_m[1] = -beamBeta * e(n).data_m[0] / Physics::c;
    b(n).data_m[2] = 0.0;
}

KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::convertFromFrenetSerret(
        size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit, const Vector3D_t& tUnit,
        const VectorView_t& e, const VectorView_t& b) {
    e(n) = e(n).data_m[0] * bUnit + e(n).data_m[1] * nUnit + e(n).data_m[2] * tUnit;
    b(n) = b(n).data_m[0] * bUnit + b(n).data_m[1] * nUnit + b(n).data_m[2] * tUnit;
}


}  // namespace opalx::spacecharge

#endif
