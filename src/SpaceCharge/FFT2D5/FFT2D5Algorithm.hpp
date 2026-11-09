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
 * @brief FFT2D5 deposition, line-density reduction, transverse solves and field gathering.
 *
 * Uses common-frame path projections, laboratory-normalized charge and the integrated
 * free-space transverse Green function. Diagnostic hooks expose the same numerical stages
 * to independent tests; no extra physical solver mode is selected here.
 */
#ifndef OPALX_FFT2D5_ALGORITHM_HPP
#define OPALX_FFT2D5_ALGORITHM_HPP

#include <Kokkos_NumericTraits.hpp>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "VectorMath.h"

namespace opalx::spacecharge {

    /**
     * @brief Device kernels shared by the FFT2D5 library and its diagnostic callers.
     *
     * Named functors avoid NVCC extended-lambda copy helpers whose private state can disagree
     * when the host linker selects a wrapper from another translation unit or shared object.
     * Views and diagnostic policies are copied by value, as in the original lambda captures.
     */
    namespace fft2d5_detail {
        using Algorithm = FFT2D5Algorithm;

        /** @brief Deposit particle charge using trilinear or single-slice bilinear CIC. */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy, typename InverseSpacing>
        struct ScatterToGridFunctor {
            Algorithm::VectorView_t r_m;
            Algorithm::VectorView_t p_m;
            Algorithm::ReferenceView_t ref_m;
            Algorithm::SegmentGeometryView_t geometry_m;
            Algorithm::VectorView_t projectedR_m;
            Algorithm::SegmentIndexView_t projectedSegment_m;
            Algorithm::T meanPs_m;
            Algorithm::ScalarView_t dt_m;
            Algorithm::BooleanView_t invalid_m;
            InverseSpacing invDr_m;
            int nGhost_m;
            ippl::NDIndex<3U> lDom_m;
            Algorithm::ScalarGridView3D_t rho_m;
            Algorithm::Vector3D_t origin_m;
            DiagnosticPolicy diagnostic_m;

            /** @brief Deposit one local particle through the captured projection/grid policy; n is
             * the local index. */
            KOKKOS_FUNCTION void operator()(const size_t n) const {
                Algorithm::doScatterToGrid<ScatterLongitudinally>(
                        n, r_m, p_m, ref_m, meanPs_m, dt_m, invalid_m, invDr_m, nGhost_m, lDom_m,
                        rho_m, origin_m, diagnostic_m, geometry_m, projectedR_m,
                        projectedSegment_m);
            }
        };

        /** @brief Reject open-path endpoint projection that would collapse source charge. */
        struct CheckPathCoverageFunctor {
            Algorithm::VectorView_t r_m;
            Algorithm::VectorView_t p_m;
            Algorithm::ReferenceView_t ref_m;
            Algorithm::SegmentGeometryView_t geometry_m;
            Algorithm::VectorView_t projectedR_m;
            Algorithm::SegmentIndexView_t projectedSegment_m;
            Algorithm::BooleanView_t invalid_m;
            Algorithm::T length_m;
            Algorithm::T tolerance_m;

            /** @brief Add one to outside for an active particle projected beyond either open
             * endpoint. */
            KOKKOS_FUNCTION void operator()(const size_t n, size_t& outside) const {
                if (invalid_m(n)) return;
                Algorithm::Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
                Algorithm::convertToFrenetSerret(
                        n, r_m, p_m, ref_m, fsR, fsP, bUnit, nUnit, tUnit, geometry_m, projectedR_m,
                        projectedSegment_m);
                const Algorithm::Vector3D_t fromStart = r_m(n) - ref_m(0);
                const auto before                     = fromStart.dot(tUnit);
                const Algorithm::Vector3D_t fromEnd   = r_m(n) - ref_m(ref_m.extent(0) - 1);
                const auto after                      = fromEnd.dot(tUnit);
                if ((fsR[2] <= tolerance_m && before < -tolerance_m)
                    || (fsR[2] >= length_m - tolerance_m && after > tolerance_m)) {
                    ++outside;
                }
            }
        };

        /**
         * @brief Mark slices containing any nonzero deposited charge, including neutral sources.
         * @details This kernel is launched from a diagnostic-policy template instantiated in
         * both the solver library and unit/diagnostic translation units. A named functor keeps
         * its captured views in an ordinary public type instead of NVCC extended-lambda copy
         * helpers whose per-translation-unit state can disagree after host linking.
         */
        struct OccupiedSlicesFunctor {
            Algorithm::ScalarGridView3D_t rho_m;
            Kokkos::View<int*> occupied_m;
            int nGhost_m;

            /** @brief Mark plane z occupied when transverse row x contains any nonzero charge. */
            KOKKOS_FUNCTION void operator()(const size_t z, const size_t x) const {
                for (size_t y = 0; y < rho_m.extent(1); ++y) {
                    // Net charge can vanish while opposite charges still generate a field.
                    if (rho_m(x, y, z + nGhost_m) != Algorithm::T(0)) {
                        Kokkos::atomic_max(&occupied_m(z), 1);
                        break;
                    }
                }
            }
        };

        /** @brief Fold deposited ghost charge into the opposite end of a closed ring. */
        struct ScatterBoundariesFunctor {
            Algorithm::ScalarGridView3D_t rhoView_m;
            size_t firstRealZ_m;
            size_t lastRealZ_m;

            /** @brief Fold the two deposited z ghosts for transverse cell (i,j), then clear them.
             */
            KOKKOS_FUNCTION void operator()(const size_t i, const size_t j) const {
                rhoView_m(i, j, firstRealZ_m) += rhoView_m(i, j, lastRealZ_m + 1);
                rhoView_m(i, j, lastRealZ_m) += rhoView_m(i, j, firstRealZ_m - 1);
                rhoView_m(i, j, lastRealZ_m + 1)  = 0;
                rhoView_m(i, j, firstRealZ_m - 1) = 0;
            }
        };

        /** @brief Convert deposited charge times time step to charge density. */
        struct ScaleChargeDensityFunctor {
            Algorithm::ScalarGridView3D_t rho_m;
            Algorithm::T scale_m;

            /** @brief Divide cell (i,j,k), including source ghosts, by the captured volume/time
             * factor. */
            KOKKOS_FUNCTION void operator()(const size_t i, const size_t j, const size_t k) const {
                rho_m(i, j, k) /= scale_m;
            }
        };

        /** @brief Sum the charge density over one transverse slice. */
        struct SliceDensitySumFunctor {
            Algorithm::ScalarGridView3D_t rho3d_m;
            size_t k_m;
            Algorithm::LineDensityView_t sums_m;

            /** @brief Reduce one local z plane per team, including deposited source ghosts.
             * Vector lanes traverse the faster transverse dimension; team threads cover
             * rows. Each team writes a distinct local slot, never a wrapped global slot:
             * tiny periodic grids can map both ghosts to the same global slice. */
            KOKKOS_FUNCTION void operator()(const Kokkos::TeamPolicy<>::member_type& team) const {
                const size_t k    = k_m + team.league_rank();
                const bool xFast  = rho3d_m.stride(0) < rho3d_m.stride(1);
                const size_t fast = rho3d_m.extent(xFast ? 0 : 1);
                const size_t slow = rho3d_m.extent(xFast ? 1 : 0);
                Algorithm::T sum{};
                Kokkos::parallel_reduce(
                        Kokkos::TeamThreadRange(team, slow),
                        [&](const size_t row, Algorithm::T& rows) {
                            Algorithm::T value{};
                            Kokkos::parallel_reduce(
                                    Kokkos::ThreadVectorRange(team, fast),
                                    [&](const size_t col, Algorithm::T& cells) {
                                        cells += rho3d_m(xFast ? col : row, xFast ? row : col, k);
                                    },
                                    value);
                            rows += value;
                        },
                        sum);
                Kokkos::single(Kokkos::PerTeam(team), [&] {
                    sums_m(k) = sum;
                });
            }
        };

        /** @brief Integrate the transverse cell area into the line-charge density. */
        struct ScaleLineDensityFunctor {
            Algorithm::LineDensityView_t lineDensity_m;
            Algorithm::T dx_m;
            Algorithm::T dy_m;

            /** @brief Multiply global profile entry k, including boundary ghosts, by transverse
             * cell area. */
            KOKKOS_FUNCTION void operator()(const size_t k) const {
                lineDensity_m(k) *= dx_m * dy_m;
            }
        };

        /** @brief Central-difference the line density using its boundary ghost cells. */
        struct LineDensityGradientFunctor {
            Algorithm::LineDensityView_t lineDensity_m;
            Algorithm::LineDensityView_t gradient_m;
            size_t firstRealCell_m;
            Algorithm::T dz_m;

            /** @brief Set real-slice gradient k from its two neighboring line-density entries. */
            KOKKOS_FUNCTION void operator()(const size_t k) const {
                gradient_m(k) = (lineDensity_m(k + firstRealCell_m + 1)
                                 - lineDensity_m(k + firstRealCell_m - 1))
                                / (2.0 * dz_m);
            }
        };

        /** @brief Apply the electrostatic coupling constant before a slice solve. */
        struct PoissonCouplingFunctor {
            Algorithm::ScalarGridView2D_t rho2d_m;

            /** @brief Apply inverse permittivity to source cell (i,j) before a serial slice solve.
             */
            KOKKOS_FUNCTION void operator()(const size_t i, const size_t j) const {
                rho2d_m(i, j) *= 1 / Physics::epsilon_0;
            }
        };

        /** @brief Copy a solved transverse field into its longitudinal grid slice. */
        struct CopySliceFieldFunctor {
            Algorithm::VectorGridView2D_t e2d_m;
            Algorithm::VectorGridView3D_t e3d_m;
            size_t z_m;

            /** @brief Copy transverse field cell (i,j) to the captured z plane and zero its z
             * component. */
            KOKKOS_FUNCTION void operator()(const size_t i, const size_t j) const {
                e3d_m(i, j, z_m)[0] = e2d_m(i, j)[0];
                e3d_m(i, j, z_m)[1] = e2d_m(i, j)[1];
                e3d_m(i, j, z_m)[2] = 0.0;
            }
        };

        /** @brief Set longitudinal field ghosts to periodic copies or zero for an open path. */
        template <bool ClosedRing>
        struct FieldGhostsFunctor {
            Algorithm::VectorGridView3D_t e3d_m;
            size_t leftGhost_m;
            size_t rightGhost_m;

            /** @brief Set both z field ghosts for transverse cell (i,j) using the selected boundary
             * mode. */
            KOKKOS_FUNCTION void operator()(const size_t i, const size_t j) const {
                if constexpr (ClosedRing) {
                    e3d_m(i, j, leftGhost_m)[0]  = e3d_m(i, j, rightGhost_m - 1)[0];
                    e3d_m(i, j, leftGhost_m)[1]  = e3d_m(i, j, rightGhost_m - 1)[1];
                    e3d_m(i, j, leftGhost_m)[2]  = e3d_m(i, j, rightGhost_m - 1)[2];
                    e3d_m(i, j, rightGhost_m)[0] = e3d_m(i, j, leftGhost_m + 1)[0];
                    e3d_m(i, j, rightGhost_m)[1] = e3d_m(i, j, leftGhost_m + 1)[1];
                    e3d_m(i, j, rightGhost_m)[2] = e3d_m(i, j, leftGhost_m + 1)[2];
                } else {
                    e3d_m(i, j, leftGhost_m)[0]  = 0;
                    e3d_m(i, j, leftGhost_m)[1]  = 0;
                    e3d_m(i, j, leftGhost_m)[2]  = 0;
                    e3d_m(i, j, rightGhost_m)[0] = 0;
                    e3d_m(i, j, rightGhost_m)[1] = 0;
                    e3d_m(i, j, rightGhost_m)[2] = 0;
                }
            }
        };

        /** @brief Gather slice fields and transform the particle fields back to lab coordinates. */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy, typename InverseSpacing>
        struct GatherFromGridFunctor {
            Algorithm::VectorView_t r_m;
            Algorithm::VectorView_t p_m;
            Algorithm::ReferenceView_t ref_m;
            Algorithm::SegmentGeometryView_t geometry_m;
            Algorithm::VectorView_t projectedR_m;
            Algorithm::SegmentIndexView_t projectedSegment_m;
            Algorithm::T beamGamma_m;
            Algorithm::T beamBeta_m;
            Algorithm::VectorView_t e_m;
            Algorithm::VectorView_t b_m;
            Algorithm::BooleanView_t invalid_m;
            InverseSpacing invDr_m;
            int nGhost_m;
            ippl::NDIndex<3U> lDom_m;
            Algorithm::VectorGridView3D_t eField_m;
            Algorithm::Vector3D_t origin_m;
            Algorithm::T gBy4PiEpsilon0_m;
            Algorithm::LineDensityView_t lineDensityGradient_m;
            bool closedRing_m;
            DiagnosticPolicy diagnostic_m;

            /** @brief Gather and rotate fields for local particle n using the captured
             * interpolation policy. */
            KOKKOS_FUNCTION void operator()(const size_t n) const {
                Algorithm::doGatherFromGrid<ScatterLongitudinally>(
                        n, r_m, p_m, ref_m, beamGamma_m, beamBeta_m, e_m, b_m, invalid_m, invDr_m,
                        nGhost_m, lDom_m, eField_m, origin_m, gBy4PiEpsilon0_m,
                        lineDensityGradient_m, closedRing_m, diagnostic_m, geometry_m, projectedR_m,
                        projectedSegment_m);
            }
        };
    }  // namespace fft2d5_detail

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::doRunSolver(
            const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic) {
        preparePathGeometryAndProjections(context);
        if (ippl::Comm->size() > 1) migrateToSlabs(context);
        scatterToGrid<DiagnosticPolicy>(context, diagnostic, true);
        solvePoissons<DiagnosticPolicy>(diagnostic);
        if (ippl::Comm->size() == 1) calculateLineDensity<DiagnosticPolicy>(diagnostic);
        gatherFromGrid<DiagnosticPolicy>(context, diagnostic, true);
    }

    template <typename DiagnosticPolicy>
    std::unique_ptr<DiagnosticPolicy> FFT2D5Algorithm::createDiagnostic(NullDiagnostic::Kind kind) {
        auto result = std::make_unique<DiagnosticPolicy>(kind);
        result->initialise(
                particles_m, fieldStorage_m->chargeDensity(), lineDensity_m, lineDensityGradient_m,
                fieldStorage_m->electricField());
        return result;
    }

    // Place the charge from all the particles into the charge density grid.
    // Multiple particle containers in the bunch are supported.
    // The particle container array dt is temporarily used to contain the charge to deposit.
    // The real work is handled by the doScatterToGrid Kokkos kernel.
    // The periodic boundary conditions are handled when the ring is marked as closed.
    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::scatterToGrid(
            const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic,
            bool usePreparedProjections) {
        fft2d5_detail::timed("2D5 Scatter", [&] {
            if (solveReferencePath_m.extent(0) < 2) return;
            Kokkos::deep_copy(fieldStorage_m->chargeDensity().getView(), 0.);
            for (size_t c = 0; c < particles_m.size(); ++c)
                if (context.trackingActive()[c])
                    scatterContainer(c, diagnostic, usePreparedProjections);
            finishChargeDeposition(context, diagnostic);
        });
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::scatterContainer(
            size_t container, DiagnosticPolicy diagnostic, bool usePreparedProjections) {
        const auto ref      = solveReferencePath_m;
        const auto invDr    = 1.0 / fieldStorage_m->spacing();
        const int nGhost    = fieldStorage_m->chargeDensity().getNghost();
        const auto lDom     = fieldStorage_m->chargeDensity().getLayout().getLocalNDIndex();
        const auto rhoView  = fieldStorage_m->chargeDensity().getView();
        const auto origin   = fieldStorage_m->origin();
        auto* pc            = particles_m[container];
        const auto geometry = usePreparedProjections ? segmentGeometry_m : SegmentGeometryView_t();
        const auto projectedR =
                usePreparedProjections ? projectedPositions_m[container] : VectorView_t();
        const auto projectedSegment =
                usePreparedProjections ? projectedSegments_m[container] : SegmentIndexView_t();
        validateDepositionContainer(container, usePreparedProjections);
        pc->scaleDtByCharge();
        const auto r       = pc->R.getView();
        const auto p       = pc->P.getView();
        const auto meanPs  = pc->getMeanP().data_m[2];
        const auto dt      = pc->dt.getView();
        const auto invalid = pc->InvalidMask.getView();
        if (config_m.scatterLongitudinally) {
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                    fft2d5_detail::ScatterToGridFunctor<true, DiagnosticPolicy, decltype(invDr)>{
                            r, p, ref, geometry, projectedR, projectedSegment, meanPs, dt, invalid,
                            invDr, nGhost, lDom, rhoView, origin, diagnostic});
        } else {
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid::scatter", pc->getLocalNum(),
                    fft2d5_detail::ScatterToGridFunctor<false, DiagnosticPolicy, decltype(invDr)>{
                            r, p, ref, geometry, projectedR, projectedSegment, meanPs, dt, invalid,
                            invDr, nGhost, lDom, rhoView, origin, diagnostic});
        }
        Kokkos::fence();
        diagnostic.scatterCharge(rhoView);
        pc->unscaleDtByCharge();
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::finishChargeDeposition(
            const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic) {
        using Policy2D_t   = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        using Policy3D_t   = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
        const auto rhoView = fieldStorage_m->chargeDensity().getView();
        // Handle the closed ring periodic boundary condition
        if (ippl::Comm->size() == 1 && config_m.closedRing) {
            constexpr auto firstRealZ = LineDensityFirstRealCell;
            const auto lastRealZ = fieldStorage_m->meshSize()[2] - 1 + LineDensityFirstRealCell;
            Kokkos::parallel_for(
                    "Solve2d5::scatterToGrid::boundaries",
                    Policy2D_t({0, 0}, {rhoView.extent(0), rhoView.extent(1)}),
                    fft2d5_detail::ScatterBoundariesFunctor{rhoView, firstRealZ, lastRealZ});
            Kokkos::fence();
        }
        // Now scale by volume and time step to get the proper charge density.
        const auto cellVolume = std::reduce(
                fieldStorage_m->spacing().begin(), fieldStorage_m->spacing().end(), 1.0,
                std::multiplies());
        const auto scale = context.stepState().timeStep * cellVolume;
        if (!std::isfinite(scale) || scale == 0.0) {
            throw OpalException(
                    "FFT2D5Algorithm::scatterToGrid",
                    "FFT2D5 deposition requires a finite nonzero time step and cell "
                    "volume.");
        }
        Kokkos::parallel_for(
                "Solve2d5::scatterToGrid::scale",
                Policy3D_t({0, 0, 0}, {rhoView.extent(0), rhoView.extent(1), rhoView.extent(2)}),
                fft2d5_detail::ScaleChargeDensityFunctor{rhoView, scale});
        Kokkos::fence();
        if (ippl::Comm->size() > 1) {
            // Sum each rank's deposits before halo accumulation. This includes
            // transverse CIC spill for the original line-charge convention,
            // while the transverse free-space FFT uses physical x-y cells only.
            calculateLineDensity<DiagnosticPolicy>(diagnostic);
            accumulateDensityHalos();
        }
        diagnostic.scatterChargeDensity(rhoView);
    }

    // This is a Kokkos kernel that deposits the charge for a single particle into the
    // charge density grid.  The particle's coordinates are first translated into
    // the Frenet-Serret coordinate system using the previously loaded reference path.
    // Momentum is transformed for diagnostics; positions and deposited line charge
    // remain in the lab frame. Charge is scattered with a CiC algorithm selected by
    // template parameter.   When true, trilinear CiC is used.  When false, the particle
    // is assigned to a single z slice and bilinear CiC is used in that slice.
    template <bool ScatterLongitudinally, typename DiagnosticPolicy>
    KOKKOS_FUNCTION void FFT2D5Algorithm::doScatterToGrid(
            const size_t n, const VectorView_t& r, const VectorView_t& p,
            const ReferenceView_t& ref, const T meanPs, const ScalarView_t& dt,
            const BooleanView_t& invalid, const Vector3D_t invDr, const int nGhost,
            const ippl::NDIndex<3U> lDom, const ScalarGridView3D_t rho, const Vector3D_t origin,
            DiagnosticPolicy diagnostic, const SegmentGeometryView_t& geometry,
            const VectorView_t& projectedR, const SegmentIndexView_t& projectedSegment) {
        if (!invalid(n)) {
            // Into Frenet-Serret coordinates
            Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
            convertToFrenetSerret(
                    n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit, geometry, projectedR,
                    projectedSegment);
            diagnostic.frenetSerretScatter(n, fsR, fsP, invalid(n));
            // Boost into the beam frame
            boostToBeamFrame(meanPs, fsP);
            diagnostic.boostToBeam(n, fsR, fsP, invalid(n));
            // CiC scatter the charge to the 3D rho grid
            scatterToRho<ScatterLongitudinally>(n, fsR, dt, invDr, nGhost, lDom, rho, origin);
        }
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::convertToFrenetSerret(
            const size_t n, const VectorView_t& r, const VectorView_t& p,
            const ReferenceView_t& ref, Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit,
            Vector3D_t& nUnit, Vector3D_t& tUnit, const SegmentGeometryView_t& geometry,
            const VectorView_t& projectedR, const SegmentIndexView_t& projectedSegment) {
        if (projectedR.extent(0) != 0) {
            loadPreparedCoordinates(
                    n, p, geometry, projectedR, projectedSegment, fsR, fsP, bUnit, nUnit, tUnit);
            return;
        }
        // Find the segment with the shortest normal to the point
        T bestDist2 = Kokkos::Experimental::finite_max_v<T>;
        size_t bestI{};
        T bestU{};
        Vector3D_t bestDi{};
        Vector3D_t bestRc{};
        const auto rN       = r(n);
        const auto segments = ref.extent(0);
        for (size_t i = 0; i < segments - 1; ++i) {
            auto refI     = ref(i);
            auto refIp1   = ref(i + 1);
            Vector3D_t di = refIp1 - refI;
            const T di2   = di.dot(di);
            if (di2 > 0.0) {
                const T u       = Kokkos::clamp(di.dot(rN - refI) / di2, 0.0, 1.0);
                Vector3D_t rc   = refI + u * di;
                Vector3D_t diff = rN - rc;
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
        const auto gammaB = Kokkos::sqrt(1.0 + meanPs * meanPs);
        const auto gamma  = Kokkos::sqrt(1.0 + fsP.data_m[2] * fsP.data_m[2]);
        fsP.data_m[2]     = gammaB * fsP.data_m[2] - meanPs * gamma;
    }

    template <typename ViewType>
    KOKKOS_FUNCTION bool FFT2D5Algorithm::makeWeights(
            const Vector3D_t fsR, const Vector3D_t origin, const Vector3D_t invDr, const int nGhost,
            const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
            ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args) {
        const auto l                      = (fsR - origin) * invDr + 0.5;
        const ippl::Vector<int, 3U> index = l;
        whi                               = l - index;
        wlo                               = 1.0 - whi;
        args                              = index - lDom.first() + nGhost;
        // CIC touches args[d] and args[d] - 1, so valid args are
        // [1, extent - 1]. Anything outside would underflow or
        // overrun the field view on the device.
        bool inBounds = args[0] > 0 && args[0] < static_cast<int>(view.extent(0));
        inBounds      = inBounds && args[1] > 0 && args[1] < static_cast<int>(view.extent(1));
        inBounds      = inBounds && args[2] > 0 && args[2] < static_cast<int>(view.extent(2));
        return inBounds;
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void FFT2D5Algorithm::scatterToRho(
            const size_t n, const Vector3D_t fsR, const ScalarView_t& dt, const Vector3D_t invDr,
            const int nGhost, const ippl::NDIndex<3U>& lDom, const ScalarGridView3D_t rho,
            const Vector3D_t origin) {
        // CiC scatter the charge to the 3D rho grid
        ippl::Vector<T, Dim> whi, wlo;
        ippl::Vector<int, Dim> args;
        if (makeWeights(fsR, origin, invDr, nGhost, lDom, rho, whi, wlo, args)) {
            if constexpr (ScatterLongitudinally) {
                scatter3D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
            } else {
                scatter2D(rho, wlo, whi, args[0], args[1], args[2], dt(n));
            }
        }
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::calculateLineDensity(DiagnosticPolicy diagnostic) {
        fft2d5_detail::timed("2D5 Line density", [&] {
            if (config_m.longitudinalFieldMode == LongitudinalFieldMode::None) {
                Kokkos::deep_copy(lineDensity_m, 0);
                diagnostic.totalDensity(lineDensity_m);
                diagnostic.lineDensity(lineDensity_m);
                Kokkos::deep_copy(lineDensityGradient_m, 0);
                diagnostic.lineDensityGradient(lineDensityGradient_m);
                return;
            }
            sumLineDensityProfile();
            diagnostic.totalDensity(lineDensity_m);
            scaleLineDensityAndGradient(diagnostic);
        });
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::scaleLineDensityAndGradient(DiagnosticPolicy diagnostic) {
        const auto deviceLineDensity = lineDensity_m;
        const auto numSlices         = fieldStorage_m->meshSize()[2];
        // Convert this to line density
        const auto dx = fieldStorage_m->spacing()[0];
        const auto dy = fieldStorage_m->spacing()[1];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::convert", numSlices + LineDensityGhostCells,
                fft2d5_detail::ScaleLineDensityFunctor{deviceLineDensity, dx, dy});
        Kokkos::fence();
        diagnostic.lineDensity(lineDensity_m);
        // Find the gradient
        const auto lineDensityGradient = lineDensityGradient_m;
        const auto dz = fieldStorage_m->chargeDensity().get_mesh().getMeshSpacing().data_m[2];
        Kokkos::parallel_for(
                "Solve2d5::calculateLineDensity::gradient", numSlices,
                fft2d5_detail::LineDensityGradientFunctor{
                        deviceLineDensity, lineDensityGradient, LineDensityFirstRealCell, dz});
        Kokkos::fence();
        diagnostic.lineDensityGradient(lineDensityGradient_m);
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::solvePoissons(DiagnosticPolicy diagnostic) {
        fft2d5_detail::timed("2D5 Poisson", [&] {
            using Policy         = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
            auto e3d             = fieldStorage_m->electricField().getView();
            const auto rho3d     = fieldStorage_m->chargeDensity().getView();
            const auto nGhost    = fieldStorage_m->electricField().getNghost();
            const auto numSlices = fieldStorage_m->layout().getLocalNDIndex()[2].length();
            fieldStorage_m->resetBackendSolves();
            Kokkos::deep_copy(e3d, Vector3D_t(0.));
            Kokkos::View<int*> occupiedSlices("FFT2D5OccupiedSlices", numSlices);
            Kokkos::deep_copy(occupiedSlices, 0);
            Kokkos::parallel_for(
                    "Solve2d5::solvePoissons::occupancy",
                    Policy({0, 0}, {numSlices, rho3d.extent(0)}),
                    fft2d5_detail::OccupiedSlicesFunctor{rho3d, occupiedSlices, nGhost});
            auto occupiedHost =
                    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), occupiedSlices);
            bool batched = false;
            if constexpr (std::is_same_v<DiagnosticPolicy, NullDiagnostic>) {
                if (fieldStorage_m->batchSize() > 0) {
                    std::vector<size_t> occupied;
                    for (size_t z = 0; z < size_t(numSlices); ++z)
                        if (occupiedHost(z)) occupied.push_back(z);
                    fieldStorage_m->solveBatchedSlices(occupied);
                    batched = true;
                }
            }
            if (!batched && ippl::Comm->size() > 1) {
                solveDiagnosticSlabs(occupiedHost, diagnostic);
                batched = true;
            }
            if (!batched) solveSerialSlices(occupiedHost, diagnostic);
            setTransverseFieldGhosts();
            diagnostic.eField(e3d);
        });
    }

    template <typename OccupiedView, typename DiagnosticPolicy>
    void FFT2D5Algorithm::solveDiagnosticSlabs(
            OccupiedView occupiedHost, DiagnosticPolicy diagnostic) {
        const auto nGhost    = fieldStorage_m->electricField().getNghost();
        const auto numSlices = fieldStorage_m->layout().getLocalNDIndex()[2].length();
        std::vector<size_t> occupied;
        for (size_t z = 0; z < size_t(numSlices); ++z)
            if (occupiedHost(z)) occupied.push_back(z);
        FFT2D5FieldStorage::ScalarField3 potential(
                fieldStorage_m->mesh(), fieldStorage_m->layout());
        fieldStorage_m->solveBatchedSlices(occupied, &potential);
        const auto phi   = potential.getView();
        const auto first = fieldStorage_m->layout().getLocalNDIndex()[2].first();
        for (size_t z = 0; z < size_t(numSlices); ++z) {
            const auto slice = Kokkos::subview(phi, Kokkos::ALL(), Kokkos::ALL(), z + nGhost);
            diagnostic.potential(slice, z + first + nGhost);
        }
    }

    template <typename OccupiedView, typename DiagnosticPolicy>
    void FFT2D5Algorithm::solveSerialSlices(
            OccupiedView occupiedHost, DiagnosticPolicy diagnostic) {
        using Policy         = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
        auto e3d             = fieldStorage_m->electricField().getView();
        const auto rho3d     = fieldStorage_m->chargeDensity().getView();
        const auto nGhost    = fieldStorage_m->electricField().getNghost();
        const auto numSlices = fieldStorage_m->layout().getLocalNDIndex()[2].length();
        const auto& s        = fieldStorage_m->workspace();
        auto rho2d           = s.chargeDensity->getView();
        for (size_t z = 0; z < numSlices; ++z) {
            // Do the 2D solve to get Ex and Ey
            Kokkos::deep_copy(
                    rho2d, Kokkos::subview(rho3d, Kokkos::ALL(), Kokkos::ALL(), z + nGhost));
            if (occupiedHost(z) == 0) {
                diagnostic.potential(rho2d, z + nGhost);
                continue;
            }
            // Scale by the coupling constant then solve
            Kokkos::parallel_for(
                    "Solve2d5::solvePoissons::coupling",
                    Policy({0, 0}, {rho2d.extent(0), rho2d.extent(1)}),
                    fft2d5_detail::PoissonCouplingFunctor{rho2d});
            fieldStorage_m->solveSlice();
            Kokkos::fence();
            diagnostic.potential(rho2d, z + nGhost);
            auto e2d = s.electricField->getView();
            // Copy the 2D E field into the 3D E field grid
            Kokkos::parallel_for(
                    "Solve2d5::solvePoissons::copy", Policy({0, 0}, {e2d.extent(0), e2d.extent(1)}),
                    fft2d5_detail::CopySliceFieldFunctor{e2d, e3d, z + nGhost});
            Kokkos::fence();
        }
    }

    template <typename DiagnosticPolicy>
    void FFT2D5Algorithm::gatherFromGrid(
            const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic,
            bool usePreparedProjections) {
        fft2d5_detail::timed("2D5 Gather", [&] {
            if (solveReferencePath_m.extent(0) < 2) return;
            for (size_t c = 0; c < particles_m.size(); ++c) {
                if (!context.trackingActive()[c]) continue;
                if (config_m.scatterLongitudinally)
                    gatherContainer<true>(c, diagnostic, usePreparedProjections);
                else
                    gatherContainer<false>(c, diagnostic, usePreparedProjections);
            }
        });
    }

    template <bool ScatterLongitudinally, typename DiagnosticPolicy>
    void FFT2D5Algorithm::gatherContainer(
            size_t container, DiagnosticPolicy diagnostic, bool usePreparedProjections) {
        auto* pc            = particles_m[container];
        const auto geometry = usePreparedProjections ? segmentGeometry_m : SegmentGeometryView_t();
        const auto projectedR =
                usePreparedProjections ? projectedPositions_m[container] : VectorView_t();
        const auto projectedSegment =
                usePreparedProjections ? projectedSegments_m[container] : SegmentIndexView_t();
        const auto r                   = pc->R.getView();
        const auto p                   = pc->P.getView();
        const auto ref                 = solveReferencePath_m;
        const auto meanPs              = pc->getMeanP().data_m[2];
        const auto e                   = pc->E.getView();
        const auto b                   = pc->B.getView();
        const auto invalid             = pc->InvalidMask.getView();
        const auto dr                  = fieldStorage_m->spacing();
        const auto invDr               = 1.0 / dr;
        const int nGhost               = fieldStorage_m->chargeDensity().getNghost();
        const auto& layout             = fieldStorage_m->chargeDensity().getLayout();
        const auto lDom                = layout.getLocalNDIndex();
        const auto origin              = fieldStorage_m->origin();
        const auto gammaB              = Kokkos::sqrt(1.0 + meanPs * meanPs);
        const auto betaB               = meanPs / gammaB;
        const auto lineDensityGradient = lineDensityGradient_m;
        const auto eField              = fieldStorage_m->electricField().getView();
        const T gBy4PiEpsilon0         = longitudinalCoupling();
        const auto kernel              = fft2d5_detail::GatherFromGridFunctor<
                             ScatterLongitudinally, DiagnosticPolicy, decltype(invDr)>{
                r,
                p,
                ref,
                geometry,
                projectedR,
                projectedSegment,
                gammaB,
                betaB,
                e,
                b,
                invalid,
                invDr,
                nGhost,
                lDom,
                eField,
                origin,
                gBy4PiEpsilon0,
                lineDensityGradient,
                config_m.closedRing,
                diagnostic};
        Kokkos::parallel_for("Solver2d5::GatherFromGrid", pc->getLocalNum(), kernel);
        Kokkos::fence();
    }

    template <bool ScatterLongitudinally, typename DiagnosticPolicy>
    KOKKOS_FUNCTION void FFT2D5Algorithm::doGatherFromGrid(
            const size_t n, const VectorView_t& r, const VectorView_t& p,
            const ReferenceView_t& ref, const T beamGamma, const T beamBeta, const VectorView_t& e,
            const VectorView_t& b, const BooleanView_t& invalid, const Vector3D_t invDr,
            const int nGhost, const ippl::NDIndex<3U> lDom, const VectorGridView3D_t eField,
            const Vector3D_t origin, const T gBy4PiEpsilon0,
            const LineDensityView_t lineDensityGradient, const bool closedRing,
            DiagnosticPolicy diagnostic, const SegmentGeometryView_t& geometry,
            const VectorView_t& projectedR, const SegmentIndexView_t& projectedSegment) {
        if (!invalid(n)) {
            // Into Frenet-Serret coordinates
            Vector3D_t fsR, fsP, bUnit, nUnit, tUnit;
            convertToFrenetSerret(
                    n, r, p, ref, fsR, fsP, bUnit, nUnit, tUnit, geometry, projectedR,
                    projectedSegment);
            diagnostic.frenetSerretGather(n, fsR, fsP, invalid(n));
            // CiC Gather the boosted E field
            gatherFromEField<ScatterLongitudinally>(n, fsR, e, invDr, nGhost, lDom, eField, origin);
            diagnostic.gatherEField(n, e(n), b(n), invalid(n));
            // Complete the lab self-fields from the lab-density transverse solve
            setTransverseLabMagneticField(n, beamBeta, e, b);
            diagnostic.deboostFromBeam(n, e(n), b(n), invalid(n));
            // Calculate the longitudinal E field in Frenet-Serret coordinates
            // Gradients live at cell centers (k+1/2)*dz, with no ghost offset.
            // Use the same longitudinal shape as charge deposition. The old
            // int(s/dz+1/2) shifted the force by half a cell and broke symmetry.
            const T cell = (fsR[2] - origin[2]) * invDr[2] - 0.5;
            int left     = static_cast<int>(Kokkos::floor(cell));
            T fraction   = cell - left;
            if constexpr (!ScatterLongitudinally) {
                left     = static_cast<int>(Kokkos::floor(cell + 0.5));
                fraction = 0;
            }
            int right       = left + 1;
            const int count = static_cast<int>(lineDensityGradient.extent(0));
            if (closedRing) {
                left  = (left % count + count) % count;
                right = (right % count + count) % count;
            } else {
                left  = Kokkos::clamp(left, 0, count - 1);
                right = Kokkos::clamp(right, 0, count - 1);
            }
            const T ldg = (1 - fraction) * lineDensityGradient(left)
                          + fraction * lineDensityGradient(right);
            e(n).data_m[2] = -gBy4PiEpsilon0 * ldg / beamGamma / beamGamma;
            diagnostic.longitudinalField(n, e(n), b(n), invalid(n));
            // And finally back into lab coordinates
            convertFromFrenetSerret(n, bUnit, nUnit, tUnit, e, b);
            diagnostic.labFrameFields(n, e(n), b(n), invalid(n));
        }
    }

    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION void FFT2D5Algorithm::gatherFromEField(
            const size_t n, const Vector3D_t& fsR, const VectorView_t& e, const Vector3D_t& invDr,
            const int nGhost, const ippl::NDIndex<3U>& lDom, const VectorGridView3D_t& eField,
            const Vector3D_t& origin) {
        // CiC gather the boosted E field to the 3D E field grid
        ippl::Vector<T, Dim> whi, wlo;
        ippl::Vector<int, Dim> args;
        if (makeWeights(fsR, origin, invDr, nGhost, lDom, eField, whi, wlo, args)) {
            if constexpr (ScatterLongitudinally) {
                e(n) = gather3D(eField, wlo, whi, args[0], args[1], args[2]);
            } else {
                e(n) = gather2D(eField, wlo, whi, args[0], args[1], args[2]);
            }
        }
    }

    KOKKOS_INLINE_FUNCTION FFT2D5Algorithm::Vector3D_t FFT2D5Algorithm::gather2D(
            const VectorGridView3D_t& eField, const ippl::Vector<T, 3U>& wlo,
            const ippl::Vector<T, 3U>& whi, const int x, const int y, const int z) {
        Vector3D_t result = wlo[0] * wlo[1] * eField(x - 1, y - 1, z);
        result += whi[0] * wlo[1] * eField(x, y - 1, z);
        result += wlo[0] * whi[1] * eField(x - 1, y, z);
        result += whi[0] * whi[1] * eField(x, y, z);
        return result;
    }

    KOKKOS_INLINE_FUNCTION FFT2D5Algorithm::Vector3D_t FFT2D5Algorithm::gather3D(
            const VectorGridView3D_t& eField, const ippl::Vector<T, 3U>& wlo,
            const ippl::Vector<T, 3U>& whi, const int x, const int y, const int z) {
        Vector3D_t result = wlo[0] * wlo[1] * wlo[2] * eField(x - 1, y - 1, z - 1);
        result += whi[0] * wlo[1] * wlo[2] * eField(x, y - 1, z - 1);
        result += wlo[0] * whi[1] * wlo[2] * eField(x - 1, y, z - 1);
        result += whi[0] * whi[1] * wlo[2] * eField(x, y, z - 1);
        result += wlo[0] * wlo[1] * whi[2] * eField(x - 1, y - 1, z);
        result += whi[0] * wlo[1] * whi[2] * eField(x, y - 1, z);
        result += wlo[0] * whi[1] * whi[2] * eField(x - 1, y, z);
        result += whi[0] * whi[1] * whi[2] * eField(x, y, z);
        return result;
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::scatter3D(
            const ScalarGridView3D_t& rho, const ippl::Vector<T, 3U>& wlo,
            const ippl::Vector<T, 3U>& whi, const int x, const int y, const int z, const T charge) {
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
            const ScalarGridView3D_t& rho, const ippl::Vector<T, 3U>& wlo,
            const ippl::Vector<T, 3U>& whi, const int x, const int y, const int z, const T charge) {
        Kokkos::atomic_add(&rho(x - 1, y - 1, z), wlo[0] * wlo[1] * charge);
        Kokkos::atomic_add(&rho(x, y - 1, z), whi[0] * wlo[1] * charge);
        Kokkos::atomic_add(&rho(x - 1, y, z), wlo[0] * whi[1] * charge);
        Kokkos::atomic_add(&rho(x, y, z), whi[0] * whi[1] * charge);
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::setTransverseLabMagneticField(
            const size_t n, const T beamBeta, const VectorView_t& e, const VectorView_t& b) {
        // The slice source uses laboratory positions and charge per laboratory
        // length. Its transverse Poisson field is therefore already the lab E:
        // the rest-density 1/gamma and the E-field boost gamma cancel. Applying
        // another gamma here would overcount that transformation.
        // B = beta * t_hat cross E / c, so v cross B cancels beta^2 of E.
        b(n).data_m[0] = -beamBeta * e(n).data_m[1] / Physics::c;
        b(n).data_m[1] = beamBeta * e(n).data_m[0] / Physics::c;
        b(n).data_m[2] = 0.0;
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::convertFromFrenetSerret(
            const size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit,
            const Vector3D_t& tUnit, const VectorView_t& e, const VectorView_t& b) {
        e(n) = e(n).data_m[0] * bUnit + e(n).data_m[1] * nUnit + e(n).data_m[2] * tUnit;
        b(n) = b(n).data_m[0] * bUnit + b(n).data_m[1] * nUnit + b(n).data_m[2] * tUnit;
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::evaluatePathSegment(
            size_t i, const Vector3D_t& position, const ReferenceView_t& ref,
            const SegmentGeometryView_t& geometry, T& bestDistance2, size_t& best, T& bestFraction,
            Vector3D_t& closest) {
        const T length2 = geometry(i, 3);
        if (length2 <= 0.0) return;
        const Vector3D_t direction(geometry(i, 0), geometry(i, 1), geometry(i, 2));
        const Vector3D_t start = ref(i);
        const T fraction       = Kokkos::clamp(direction.dot(position - start) / length2, 0.0, 1.0);
        const Vector3D_t point = start + fraction * direction;
        const Vector3D_t difference = position - point;
        const T distance2           = difference.dot(difference);
        if (distance2 < bestDistance2 || (distance2 == bestDistance2 && i < best)) {
            best          = i;
            bestFraction  = fraction;
            closest       = point;
            bestDistance2 = distance2;
        }
    }

    KOKKOS_INLINE_FUNCTION FFT2D5Algorithm::T FFT2D5Algorithm::pathDistanceRoundoff(
            const Vector3D_t& position, const SegmentGeometryView_t& bounds) {
        T scale = 0.0;
        for (unsigned d = 0; d < 3; ++d) {
            scale = Kokkos::max(scale, Kokkos::abs(position[d]));
            scale = Kokkos::max(scale, Kokkos::abs(bounds(1, d)));
            scale = Kokkos::max(scale, Kokkos::abs(bounds(1, d + 3)));
        }
        constexpr T epsilon = Kokkos::Experimental::epsilon_v<T>;
        return 256 * epsilon * scale * scale + 256 * Kokkos::Experimental::norm_min_v<T>;
    }

    KOKKOS_INLINE_FUNCTION FFT2D5Algorithm::T FFT2D5Algorithm::pathNodeDistance2(
            size_t node, const Vector3D_t& position, const SegmentGeometryView_t& bounds) {
        T distance2 = 0.0;
        for (unsigned d = 0; d < 3; ++d) {
            const T gap = Kokkos::max(
                    T(0),
                    Kokkos::max(bounds(node, d) - position[d], position[d] - bounds(node, d + 3)));
            distance2 += gap * gap;
        }
        return distance2;
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::findClosestPathSegment(
            const Vector3D_t& position, size_t seed, const ReferenceView_t& ref,
            const SegmentGeometryView_t& geometry, const SegmentGeometryView_t& bounds,
            size_t leaves, size_t& best, T& bestFraction, Vector3D_t& closest) {
        const size_t segments = ref.extent(0) - 1;
        T bestDistance2       = Kokkos::Experimental::finite_max_v<T>;
        const auto evaluate   = [&](size_t i) {
            evaluatePathSegment(
                    i, position, ref, geometry, bestDistance2, best, bestFraction, closest);
        };
        evaluate(seed < segments ? seed : 0);
        const T roundoff      = pathDistanceRoundoff(position, bounds);
        const auto lowerBound = [&](size_t node) {
            return pathNodeDistance2(node, position, bounds);
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
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::projectPreparedParticle(
            size_t n, const VectorView_t& positions, const BooleanView_t& invalid,
            const ReferenceView_t& ref, const SegmentGeometryView_t& geometry,
            const SegmentGeometryView_t& bounds, const VectorView_t& projected,
            const SegmentIndexView_t& indices, size_t leaves, T periodicLength) {
        const size_t segments = ref.extent(0) - 1;
        if (invalid(n)) {
            projected(n) = Vector3D_t(0.);
            indices(n)   = segments;
            return;
        }
        const Vector3D_t position = positions(n);
        size_t best               = 0;
        T bestFraction            = 0.;
        Vector3D_t closest(0.);
        findClosestPathSegment(
                position, indices(n), ref, geometry, bounds, leaves, best, bestFraction, closest);
        const Vector3D_t tangent =
                Vector3D_t(geometry(best, 0), geometry(best, 1), geometry(best, 2))
                / geometry(best, 4);
        const Vector3D_t normal(0.0, 1.0, 0.0);
        const Vector3D_t binormal = ippl::cross(normal, tangent);
        const Vector3D_t offset   = position - closest;
        projected(n)              = Vector3D_t(
                offset.dot(binormal), offset.dot(normal),
                geometry(best, 5) + bestFraction * geometry(best, 4));
        // Ownership wraps the exact end point to the first slab.
        // Keep the cached interpolation coordinate consistent with it.
        if (periodicLength > 0. && projected(n)[2] >= periodicLength) projected(n)[2] = 0.;
        indices(n) = best;
    }

    KOKKOS_INLINE_FUNCTION void FFT2D5Algorithm::loadPreparedCoordinates(
            size_t n, const VectorView_t& p, const SegmentGeometryView_t& geometry,
            const VectorView_t& projectedR, const SegmentIndexView_t& projectedSegment,
            Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit,
            Vector3D_t& tUnit) {
        const auto segment = projectedSegment(n);
        fsR                = projectedR(n);
        tUnit = Vector3D_t(geometry(segment, 0), geometry(segment, 1), geometry(segment, 2))
                / geometry(segment, 4);
        nUnit = {0, 1, 0};
        bUnit = ippl::cross(nUnit, tUnit);
        fsP   = Vector3D_t(p(n).dot(bUnit), p(n).dot(nUnit), p(n).dot(tUnit));
    }

}  // namespace opalx::spacecharge

#endif
