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

/** @file FFT2D5Algorithm.h @brief Adapts the original 2.5D solver to one solve context. */
#ifndef OPALX_SPACE_CHARGE_FFT2D5_ALGORITHM_H
#define OPALX_SPACE_CHARGE_FFT2D5_ALGORITHM_H

#include <memory>
#include <span>
#include <vector>
#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/FFT2D5/FFT2D5FieldStorage.h"
#include "SpaceCharge/FFT2D5/ReferencePath.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeFrames.h"

class BunchStateHandler;

namespace opalx::spacecharge {

    /**
     * @brief A 2D5 solver for electromagnetic field simulations.
     *
     * @details
     * The solver performs the following sequence of operations to account for space charge
     * in the electric and magnetic fields that particles feel.
     * 1. Map the fixed laboratory reference path and the particles' moving reference
     *    coordinates into one aligned spatial frame, then into Frenet-Serret coordinates.
     *    The reference path is provided by a file or the OrbitThreader. Only diagnostic
     *    momenta are boosted; source positions and line charge retain laboratory normalization.
     * 2. Do CiC scatter of charge into the 3D charge density grid.  This spreads the charge out
     * both transversely and (optionally) longitudinally which should be a help in reducing noise.
     * With PARFFTZ, OPEN's spatial layout migrates particles by their projected s
     * coordinate. Registered Cartesian scratch positions migrate with all attributes
     * and are restored before any force calculation. Halo accumulation transfers CIC
     * charge at slab interfaces (and the closed-ring seam); field halos serve gathering.
     * Each rank owns complete transverse planes. Only the small line-density vector
     * is globally reduced, before halo accumulation to count every deposited weight once.
     * No longitudinal FFT or IPPL modification is used.
     *
     * 3. Copy each occupied transverse charge-density slice into one reusable 2D workspace.
     * 4. Perform the 2D Poisson solve to determine the X and Y electric field components;
     *    identically empty slices have zero field and need no FFT. Boundary conditions are open.
     * 5. Calculate the longitudinal S electric field component using the algorithm from PyHEADTAIL.
     *    This uses the local line density along the reference orbit and a geometry constant that
     *    encodes the effect of the beam pipe on the longitudinal component.
     * 6. Do CiC gather of fields for the particle positions.
     * 7. Rotate the local Frenet-Serret fields into the common Cartesian frame and restore
     *    each container's moving reference coordinates and fields.
     */
    class FFT2D5Algorithm final : public SpaceChargeAlgorithm {
    public:
        using T                       = double;
        static constexpr unsigned Dim = 3;
        using ParticleContainer       = ::ParticleContainer<double>;
        using LongitudinalFieldMode   = FFT2D5LongitudinalFieldMode;
        using Mesh3D_t                = ippl::UniformCartesian<T, 3U>;
        using Mesh2D_t                = ippl::UniformCartesian<T, 2U>;
        using Vector2D_t              = Mesh2D_t::vector_type;
        using Vector3D_t              = Mesh3D_t::vector_type;
        using Layout2D_t              = ippl::FieldLayout<2U>;
        using Point_t                 = ippl::Vector<double, 3U>;
        using VectorView_t            = ParticleAttrib<Vector3D_t>::view_type;
        using ScalarView_t            = ParticleAttrib<T>::view_type;
        using BooleanView_t           = ParticleAttrib<bool>::view_type;
        using ReferenceView_t         = Kokkos::View<Vector3D_t*>;
        using LineDensityView_t       = Kokkos::View<T*>;
        using SegmentGeometryView_t   = Kokkos::View<T* [6]>;
        using SegmentIndexView_t      = Kokkos::View<size_t*>;
        using ScalarGridView3D_t      = Field<T, 3U>::view_type;
        using VectorGridView3D_t      = Field<Vector3D_t, 3U>::view_type;
        using ScalarGridView2D_t      = Field<T, 2U>::view_type;
        using VectorGridView2D_t      = Field<Vector2D_t, 2U>::view_type;

        FFT2D5Algorithm(
                FFT2D5Config config, std::span<ParticleContainer* const> particles,
                std::shared_ptr<const BunchStateHandler> bunchState);
        [[nodiscard]] SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override;
        void ensureInitialized();
        [[nodiscard]] bool initialized() const { return fieldStorage_m != nullptr; }
        [[nodiscard]] std::size_t sliceCount() const;
        [[nodiscard]] const ReferencePath::View& referencePathView() const;

        /** @brief Algorithm-local no-op hooks retained for numerical verification; allocate no
         * storage. */
        class NullDiagnostic {
        public:
            enum class Kind {
                FrenetSerretScatter,
                BoostToBeam,
                ScatterCharge,
                ScatterChargeDensity,
                TotalDensity,
                LineDensity,
                LineDensityGradient,
                EField,
                FrenetSerretGather,
                GatherEField,
                Deboosted,
                LongitudinalField,
                LabFrameFields,
                Potential
            };
            KOKKOS_FUNCTION void frenetSerretScatter(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void boostToBeam(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void scatterCharge(const ScalarGridView3D_t&) const {}
            KOKKOS_FUNCTION void scatterChargeDensity(const ScalarGridView3D_t&) const {}
            KOKKOS_FUNCTION void potential(const Field_t<2U>::view_type&, size_t) const {}
            KOKKOS_FUNCTION void eField(const VField_t<T, 3U>::view_type&) const {}
            KOKKOS_FUNCTION void totalDensity(const LineDensityView_t&) const {}
            KOKKOS_FUNCTION void lineDensity(const LineDensityView_t&) const {}
            KOKKOS_FUNCTION void lineDensityGradient(const LineDensityView_t&) const {}
            KOKKOS_FUNCTION void frenetSerretGather(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void gatherEField(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void deboostFromBeam(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void longitudinalField(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            KOKKOS_FUNCTION void labFrameFields(
                    const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
            void initialise(
                    std::span<ParticleContainer* const> /*particles*/, const Field_t<3U>& /*rho*/,
                    const LineDensityView_t& /*lineDensity*/,
                    const LineDensityView_t& /*lineDensityGradient*/,
                    const VField_t<T, 3U>& /*eField*/) {}
        };

        /**
         * @brief Internal implementation functions for the 2.5D space-charge solver.
         *
         * @details
         * These functions implement the individual stages of the 2.5D space-charge
         * calculation. They are templated on a diagnostic policy so that unit tests
         * can supply a diagnostic implementation derived from NullDiagnostic and
         * inspect intermediate results produced by the solver.
         *
         * In normal solver operation, the default NullDiagnostic policy is used,
         * which performs no diagnostic processing. These functions are considered
         * internal implementation details of the solver and its unit tests.
         *
         * @tparam DiagnosticPolicy Diagnostic policy type. The policy must provide
         *                          the diagnostic callback interface defined by
         *                          NullDiagnostic.
         */

        /**
         * @brief Executes the complete 2.5D space-charge calculation.
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect intermediate
         *                           results. Defaults to NullDiagnostic.
         * @param context Per-call activity and tracker state.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void doRunSolver(const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {});

        /**
         * @brief Deposits particle charge onto the 3D grid.
         * @param usePreparedProjections Use the cache prepared by doRunSolver; standalone
         *        diagnostic calls leave this false and calculate fresh projections.
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect the charge
         *                           deposition process. Defaults to NullDiagnostic.
         * @param context Per-call activity and tracker state.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void scatterToGrid(
                const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {},
                bool usePreparedProjections = false);

        /**
         * @brief Solves the two-dimensional Poisson equations on the longitudinal slices.
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect intermediate
         *                           solver results. Defaults to NullDiagnostic.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void solvePoissons(DiagnosticPolicy diagnostic = {});

        /**
         * @brief Calculates the longitudinal line-charge density and its gradient.
         * @details A single team reduction processes all local transverse planes, including
         * deposited source ghosts before charge-halo accumulation. One small host transfer
         * supplies the existing ordered local-to-global accumulation and MPI profile sum.
         * Local scratch slots remain distinct even when periodic ghosts wrap to the same
         * global slice. Host scratch is independent on CPU backends, where a mirror view
         * could alias the profile being cleared. Boundary conditions, cell-area scaling,
         * longitudinal gradient and diagnostic callbacks are unchanged; no longitudinal FFT
         * or IPPL modification is introduced. Scratch storage is O(global z slices).
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect the calculated
         *                           line density and its gradient. Defaults to
         *                           NullDiagnostic.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void calculateLineDensity(DiagnosticPolicy diagnostic = {});

        /**
         * @brief Interpolates the calculated electric field back to the particles.
         * @param usePreparedProjections Use only the current doRunSolver preparation;
         *        standalone probe/diagnostic calls default to fresh conversion.
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect the field
         *                           gathering process. Defaults to NullDiagnostic.
         * @param context Per-call activity and tracker state.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void gatherFromGrid(
                const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {},
                bool usePreparedProjections = false);

        /**
         * @brief Creates a diagnostic policy for use by the solver unit tests.
         *
         * This convenience function constructs a diagnostic policy configured to
         * capture information from the specified stage of the 2.5D space-charge
         * calculation.
         *
         * The returned diagnostic policy is intended to be passed to one of the
         * internal solver functions templated on DiagnosticPolicy. It is primarily
         * used by unit tests to inspect intermediate solver results for verification.
         *
         * @tparam DiagnosticPolicy The diagnostic policy type to construct. The type
         *                          is expected to derive from NullDiagnostic.
         * @param kind Identifies the stage of the solver for which diagnostic
         *             information is required.
         * @return A unique pointer to the newly constructed diagnostic policy.
         *
         * @note This is a test-support utility and is not intended to be used by
         *       normal solver clients.
         */
        template <typename DiagnosticPolicy>
        std::unique_ptr<DiagnosticPolicy> createDiagnostic(NullDiagnostic::Kind kind);

        /**
         * @brief Kokkos kernel that deposits a single particle's charge onto
         * the three-dimensional charge grid.  See implementation for more details.
         *
         * @tparam ScatterLongitudinally If true, scatter particle charge
         *                               longitudinally between adjacent slices.
         *                               If false, deposit the charge onto a single
         *                               longitudinal slice.
         * @tparam DiagnosticPolicy Diagnostic policy.
         *
         * @param n The particle number.
         * @param r Particle positions.
         * @param p Particle momenta.
         * @param ref Reference path used for the Frenet-Serret coordinate
         *            transformation.
         * @param meanPs Mean longitudinal momentum of the reference particle.
         * @param dt Particle time-step values.
         * @param invalid Flags identifying particles that should be excluded
         *                from the charge deposition.
         * @param invDr Inverse grid spacing in each spatial dimension.
         * @param nGhost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param origin Physical origin of the charge-density grid.
         * @param diagnostic Diagnostic policy used to record intermediate results.
         * @param geometry Current solve-frame segment geometry for prepared conversion.
         * @param projectedR Prepared local coordinates; empty uses the original conversion.
         * @param projectedSegment Nearest-segment indices paired with projectedR.
         */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy = NullDiagnostic>
        KOKKOS_FUNCTION static void doScatterToGrid(
                size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
                T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
                int nGhost, ippl::NDIndex<3U> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
                DiagnosticPolicy diagnostic, const SegmentGeometryView_t& geometry = {},
                const VectorView_t& projectedR             = {},
                const SegmentIndexView_t& projectedSegment = {});

        /**
         * @brief Kokkos function that converts from lab to Frenet-Serret coordinates.
         * @details Empty projection views retain the original exhaustive conversion. With
         * prepared views, reuse only position and segment selection; project the current
         * momentum onto the same basis so diagnostic momentum values are never cached.
         * @param geometry Current solve-frame segment geometry (optional cached mode).
         * @param projectedR Prepared local particle coordinates (empty selects uncached mode).
         * @param projectedSegment Prepared nearest-segment indices, paired with projectedR.
         *
         * @param n The particle number.
         * @param r Particle positions.
         * @param p Particle momenta.
         * @param ref Reference path used for the Frenet-Serret coordinate
         *            transformation.
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param fsP The particle momentum in Frenet-Serret coordinates.
         * @param bUnit The binormal unit vector.
         * @param nUnit The normal unit vector.
         * @param tUnit The tangential unit vector.
         */
        KOKKOS_FUNCTION static void convertToFrenetSerret(
                size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
                Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit,
                Vector3D_t& tUnit, const SegmentGeometryView_t& geometry = {},
                const VectorView_t& projectedR             = {},
                const SegmentIndexView_t& projectedSegment = {});

        /**
         * @brief Kokkos function that boosts coordinates into the beam reference frame.
         *
         * @param meanPs The mean momentum of the particles in the beam.
         * @param fsP The particle momentum.
         */
        KOKKOS_FUNCTION static void boostToBeamFrame(T meanPs, Vector3D_t& fsP);

        /**
         * @brief Kokkos function that scatters a particle's charge
         *
         * @tparam ScatterLongitudinally If true, scatter particle charge
         *                               longitudinally between adjacent slices.
         *                               If false, deposit the charge onto a single
         *                               longitudinal slice.
         *
         * @param n The particle number.
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param dt Particle time-step * charge.
         * @param invDr Inverse grid spacing in each spatial dimension.
         * @param nGhost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param origin Physical origin of the charge-density grid.
         */
        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void scatterToRho(
                size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, int nGhost,
                const ippl::NDIndex<3U>& lDom, ScalarGridView3D_t rho, Vector3D_t origin);

        /**
         * @brief Kokkos kernel that gathers the E and B fields for a single particle.
         * See implementation for more details.
         *
         * @tparam DiagnosticPolicy Diagnostic policy.
         * @tparam ScatterLongitudinally If true, scatter particle charge
         *                               longitudinally between adjacent slices.
         *                               If false, deposit the charge onto a single
         *                               longitudinal slice.
         *
         * @param n The particle number.
         * @param r Particle positions.
         * @param p Particle momenta.
         * @param ref Reference path used for the Frenet-Serret coordinate
         *            transformation.
         * @param beamGamma The relativistic gamma for the beam.
         * @param beamBeta The relativistic beta for the beam.
         * @param e Particle E fields.
         * @param b Particle B fields.
         * @param invalid Flags identifying particles that should be excluded
         *                from the charge deposition.
         * @param invDr Inverse grid spacing in each spatial dimension.
         * @param nGhost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param eField Electric field grid from the Poisson solver.
         * @param origin Physical origin of the charge-density grid.
         * @param gBy4PiEpsilon0 Constant for the longitudinal field estimate.
         * @param lineDensityGradient The longitudinal line charge density gradient.
         * @param closedRing Whether interpolation wraps across the path endpoints.
         * @param diagnostic Diagnostic policy used to record intermediate results.
         * @param geometry Current solve-frame segment geometry for prepared conversion.
         * @param projectedR Prepared local coordinates; empty uses the original conversion.
         * @param projectedSegment Nearest-segment indices paired with projectedR.
         */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy = NullDiagnostic>
        KOKKOS_FUNCTION static void doGatherFromGrid(
                size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
                T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b,
                const BooleanView_t& invalid, Vector3D_t invDr, int nGhost, ippl::NDIndex<3U> lDom,
                VectorGridView3D_t eField, Vector3D_t origin, T gBy4PiEpsilon0,
                LineDensityView_t lineDensityGradient, bool closedRing, DiagnosticPolicy diagnostic,
                const SegmentGeometryView_t& geometry = {}, const VectorView_t& projectedR = {},
                const SegmentIndexView_t& projectedSegment = {});

        /**
         * @brief Kokkos function that gathers the boosted E field for a particle from
         *  the electric field grid.
         *
         * @tparam ScatterLongitudinally If true, scatter particle charge
         *                               longitudinally between adjacent slices.
         *                               If false, deposit the charge onto a single
         *                               longitudinal slice.
         *
         * @param n The particle number.
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param e Particle E fields.
         * @param invDr Inverse grid spacing in each spatial dimension.
         * @param nGhost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param eField Electric field grid from the Poisson solver.
         * @param origin Physical origin of the charge-density grid.
         */
        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void gatherFromEField(
                size_t n, const Vector3D_t& fsR, const VectorView_t& e, const Vector3D_t& invDr,
                int nGhost, const ippl::NDIndex<3U>& lDom, const VectorGridView3D_t& eField,
                const Vector3D_t& origin);

        /**
         * @brief Set B = beta * t_hat cross E / c in the local path basis.
         *
         * The transverse slice solve uses lab-frame line charge, so its E is
         * already the lab electric field and must not receive another gamma.
         */
        KOKKOS_FUNCTION static void setTransverseLabMagneticField(
                size_t n, T beamBeta, const VectorView_t& e, const VectorView_t& b);

        /**
         * @brief Kokkos function that converts a particle's e and b fields from
         * Frenet-Serret coordinates into cartesian coordinates using the supplied
         * basis unit vectors calculated during the original translation to Frenet-
         * Serret coordinates.
         *
         * @param n The particle number.
         * @param bUnit The binormal unit vector.
         * @param nUnit The normal unit vector.
         * @param tUnit The tangential unit vector.
         * @param e Particle E fields.
         * @param b Particle B fields.
         */
        KOKKOS_FUNCTION static void convertFromFrenetSerret(
                size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit, const Vector3D_t& tUnit,
                const VectorView_t& e, const VectorView_t& b);

        /**
         * @brief Kokkos function that performs a bilinear gather from a 2D slice of the 3D
         * electric field grid.
         *
         * @param eField Electric field grid from the Poisson solver.
         * @param wlo Lower weighting factors.
         * @param whi Upper weighting factors.
         * @param x Horizontal transverse coordinate.
         * @param y Vertical transverse coordinate.
         * @param z Longtiduninal coordinate.
         */
        KOKKOS_FUNCTION static Vector3D_t gather2D(
                const VectorGridView3D_t& eField, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z);

        /**
         * @brief Kokkos function that performs a bilinear scatter of charge to a 2D slice of the 3D
         * charge density grid.
         *
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param wlo Lower weighting factors.
         * @param whi Upper weighting factors.
         * @param x Horizontal transverse coordinate.
         * @param y Vertical transverse coordinate.
         * @param z Longitudinal coordinate.
         * @param charge The charge to deposit.
         */
        KOKKOS_FUNCTION static void scatter2D(
                const ScalarGridView3D_t& rho, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z, T charge);

        /**
         * @brief Kokkos function that performs a trilinear gather from the 3D
         * electric field grid.
         *
         * @param eField Electric field grid from the Poisson solver.
         * @param wlo Lower weighting factors.
         * @param whi Upper weighting factors.
         * @param x Horizontal transverse coordinate.
         * @param y Vertical transverse coordinate.
         * @param z Longtiduninal coordinate.
         */
        KOKKOS_FUNCTION static Vector3D_t gather3D(
                const VectorGridView3D_t& eField, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z);

        /**
         * @brief Kokkos function that performs a trilinear scatter of charge to the 3D
         * charge density grid.
         *
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param wlo Lower weighting factors.
         * @param whi Upper weighting factors.
         * @param x Horizontal transverse coordinate.
         * @param y Vertical transverse coordinate.
         * @param z Longitudinal coordinate.
         * @param charge The charge to deposit.
         */
        KOKKOS_FUNCTION static void scatter3D(
                const ScalarGridView3D_t& rho, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z, T charge);

        /**
         * @brief Kokkos function that generates the weights for scatter/gather.
         *
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param origin Physical origin of the charge-density grid.
         * @param invDr Inverse of the grid spacing.
         * @param nGhost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param view The view being indexed.
         * @param whi Upper weighting factors.
         * @param wlo Lower weighting factors.
         * @param args The index coordinates of the cell containing fsR.
         */
        template <typename ViewType>
        KOKKOS_FUNCTION static bool makeWeights(
                Vector3D_t fsR, Vector3D_t origin, Vector3D_t invDr, int nGhost,
                const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
                ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args);

        /**
         * @return Returns the number of slices in the Frenet-Serret domain.
         * For use only by test cases.
         */
        size_t getNumSlices() const { return fieldStorage_m->meshSize()[2]; }

        /**
         * @return Returns the mesh object for the Frenet-Serret domain.
         * For use only by test cases.
         */
        Mesh2D_t* getSliceMesh() const { return &fieldStorage_m->sliceMesh(); }

        /**
         * @return Returns the slice layout object for the Frenet-Serret domain.
         * For use only by test cases.
         */
        Layout2D_t* getSliceLayout() const { return &fieldStorage_m->sliceLayout(); }

        /**
         * @return Returns the reference path view for the Frenet-Serret domain.
         * For use only by test cases.
         */
        const ReferenceView_t& getReferencePath() const { return referencePathView(); }
        /// Return segment data for independent numerical verification; handle constness does not
        /// protect elements.
        const SegmentGeometryView_t& getSegmentGeometry() const { return segmentGeometry_m; }
        /// Return the complete bounding-box heap, including padded leaves, for verification.
        const SegmentGeometryView_t& getSegmentBounds() const { return segmentBounds_m; }

        /**
         * @return Returns the 3D charge density view for the Frenet-Serret domain.
         * For use only by test cases.
         */
        Field_t<3>* getRho() const { return &fieldStorage_m->chargeDensity(); }

        /**
         * @return Returns the 3D electric field view for the Frenet-Serret domain.
         * For use only by test cases.
         */
        VField_t<T, 3>* getEField() const { return &fieldStorage_m->electricField(); }

        /**
         * @return Returns the logntudinal charge density view for the Frenet-Serret domain.
         * For use only by test cases.
         */
        const LineDensityView_t& getLineDensity() const { return lineDensity_m; }

        /**
         * @return Returns the logntudinal charge density gradient view for the Frenet-Serret
         * domain. For use only by test cases.
         */
        const LineDensityView_t& getLineDensityGradient() const { return lineDensityGradient_m; }

        /** @brief Current per-solve projection caches, exposed for numerical verification. */
        const VectorView_t& getProjectedPositions(size_t container = 0) const {
            return projectedPositions_m.at(container);
        }
        /// Return winning segment indices for a container; invalid particles carry the
        /// segment-count sentinel.
        const SegmentIndexView_t& getProjectedSegments(size_t container = 0) const {
            return projectedSegments_m.at(container);
        }

        /**
         * @brief Prepare segment geometry and particle projections for one complete field solve.
         *
         * @details
         * This single preparation call eliminates repeated path-only calculations and the
         * second nearest-segment search during gather. It is called by doRunSolver after
         * solve() has put the reference path and every active particle container into the
         * same Cartesian solve frame, and before coverage checking or charge deposition.
         * Two ordered preparation kernels precede the projection kernel for each active
         * container: independent segments/leaves, then one team for prefix and hierarchy. Geometry
         * must be complete before particles search it; no global inter-block barrier is assumed
         * inside a kernel.
         *
         * For path points R_i, segment i stores six doubles:
         * - columns 0..2: d_i = R_(i+1) - R_i;
         * - column 3: d_i dot d_i;
         * - column 4: l_i = |d_i|;
         * - column 5: S_i = sum_(j<i) l_j, the arc length at the segment start.
         * The squared length, rather than its reciprocal, retains the original division
         * in the closest-point calculation. One team thread
         * accumulates S_i in increasing segment order from the stored lengths; zero-length segments
         * contribute zero length and remain excluded from the nearest-segment search.
         *
         * For each non-invalid particle x, the unchanged segment projection is:
         * @f[
         * u_i = clamp((x-R_i)\cdot d_i / |d_i|^2,0,1),\qquad
         * q_i = R_i + u_i d_i.
         * @f]
         * A balanced implicit binary hierarchy encloses contiguous segments in
         * axis-aligned endpoint boxes. Its power-of-two leaves retain segment indices;
         * padded leaves are empty. Team threads merge one level at a time, with a team
         * barrier between levels. No inter-block barrier or parallel prefix scan is used.
         * Rebuild all boxes each solve, so frame
         * or path changes cannot leave stale bounds. Search nearer boxes first, using
         * the previous segment only as an initial upper bound (or segment zero for a
         * new/invalid cache). Movement and reordering cannot restrict the search.
         * Retain the smaller computed squared distance, breaking exact ties by the
         * lowest segment index, independently of traversal/seed order.
         *
         * A box is rejected only if its computed squared lower distance exceeds the
         * incumbent by 256*epsilon*L^2 + 256*min_normal, where L bounds absolute
         * coordinates of both the particle and all path endpoints. This conservative
         * allowance covers rounded segment directions, closest-point interpolation,
         * box subtraction and squared-distance arithmetic in three dimensions under
         * ordinary IEEE round-to-nearest arithmetic (no fast-math reassociation).
         * It also protects subnormal distances. Borderline, overflow and unordered
         * comparisons descend rather than prune. A 64-entry per-particle stack bounds
         * traversal storage; exhaustion falls back to the complete segment scan.
         * No nearest-neighbour approximation or geometry tolerance changes the winner.
         * Store its index and local coordinates (x_perp,y_perp,s), with
         * @f$ s=S_i+u_i l_i @f$. The normal remains the existing (0,1,0), the tangent is
         * d_i/l_i and the binormal is normal cross tangent. This preserves the present
         * path/basis approximation; it introduces no smoothing or curvature correction.
         * Pre-summing S_i changes floating-point association relative to adding preceding
         * lengths to u_i*l_i for every particle; small rounding differences are expected.
         *
         * The resulting cache is read by open-path coverage checking, deposition and
         * gather in this doRunSolver call. Those stages do not change particle positions,
         * momenta, ordering or the solve frame. Momentum is deliberately not cached: its
         * local components are recomputed by the conversion used for diagnostics. Invalid
         * particles receive a sentinel segment index and are skipped by all consumers.
         * Inactive containers are neither projected nor modified.
         *
         * @par Lifetime and storage
         * Recompute geometry and every active particle projection on every complete solve,
         * including after frame changes, particle movement, reordering or count changes.
         * Reuse allocated views when their extents match; otherwise resize them before
         * launching the particle kernel. Storage is O(M + N): 48*M bytes for M segments
         * plus 48*(2*P) bytes for the bounds heap, where P is the least power of two
         * not smaller than M, and (3*sizeof(double)+sizeof(size_t))*N bytes for caches.
         * Buffers of temporarily inactive containers may retain their previous allocation.
         * Views are shared read-only by OpenMP workers; no N-by-M distance matrix exists.
         * Building geometry/bounds costs O(M). The search typically rejects most
         * segments, but overlapping boxes or conservative roundoff margins retain
         * worst-case O(M*N) work. There is no assumption that path index is spatial
         * distance. This preparation is an internal stage, also exposed to numerical
         * tests just like scatterToGrid; callers must satisfy its frame preconditions.
         * Standalone scatterToGrid/gatherFromGrid diagnostic calls default to the uncached
         * conversion, so a later probe edit or exception cannot reuse a stale projection.
         *
         * @param context Per-call container activity and tracker state.
         * @param rebuildGeometry False only after migration within the same solve;
         *        reuse the unchanged path hierarchy and refresh particle projections.
         * @pre The initialized path and particles use the same aligned solve frame.
         * @pre Inspection-only diagnostic callbacks do not move or reorder particles
         *      between preparation and the end of doRunSolver.
         * @post Every active non-invalid particle has fresh local coordinates and a
         *       nearest-segment index. No particle position, momentum, field or dt is changed.
         */
        void preparePathGeometryAndProjections(
                const SpaceChargeSolveContext& context, bool rebuildGeometry = true);

        /** @brief Reuse OPEN spatial migration with logical (u,v,s) ownership.
         * Cartesian positions are saved in a registered particle attribute so IPPL
         * migrates them, IDs, momenta and charges together. R is temporarily replaced
         * by the prepared path projection, then restored exactly before any physics
         * operation. Invalid particles remain invalid. All active ranks participate,
         * including ranks with no local particles. Projection caches are rebuilt after
         * migration because local particle order and count may change.
         * @pre Particles and path are in the common solve frame; projections are current.
         * @post Tracking coordinates are Cartesian, and each particle belongs to its s slab.
         */
        void migrateToSlabs(const SpaceChargeSolveContext& context);
        /// Exchange slab halos while preserving open transverse boundaries.
        void accumulateDensityHalos();
        /// Exchange z-slab field halos, then zero transverse ghosts to keep open x/y boundaries.
        void fillFieldHalos();

        /** @internal @brief Resize segment and complete binary-heap bounds storage.
         * @return Power-of-two leaf count. @throws OpalException On index overflow. */
        size_t resizePathGeometry();
        /** @internal @brief Compute independent segment data and endpoint boxes in one kernel.
         * @param leaves Complete heap leaf count; padded leaves receive empty boxes. */
        void buildPathSegments(size_t leaves);
        /** @internal @brief Preserve ordered arc lengths and merge bounds level by level.
         * @param leaves Complete heap leaf count; preceding segment kernel shares the stream. */
        void buildPathHierarchy(size_t leaves);
        /** @internal @brief Refresh projections for one active container without moving particles.
         * @param container Container index. @param leaves Heap leaf count.
         * @param periodicLength Positive only for post-migration closed-ring seam wrapping.
         * Kernel-launching helpers are public to satisfy NVCC extended-lambda accessibility. */
        void projectContainer(size_t container, size_t leaves, T periodicLength);
        /** @internal @brief Score one nondegenerate segment with the unchanged lowest-index tie.
         * @param i Segment index. @param position Cartesian query. @param ref Path vertices.
         * @param geometry Segment data. @param bestDistance2 Current squared distance, updated.
         * @param best Winning index, updated. @param bestFraction Winning clamped fraction.
         * @param closest Winning Cartesian closest point; updated with the score. */
        KOKKOS_FUNCTION static void evaluatePathSegment(
                size_t i, const Vector3D_t& position, const ReferenceView_t& ref,
                const SegmentGeometryView_t& geometry, T& bestDistance2, size_t& best,
                T& bestFraction, Vector3D_t& closest);
        /** @internal @brief Conservative pruning allowance for rounded Cartesian arithmetic.
         * @param position Query point. @param bounds Complete heap. @return Squared-distance
         * margin. */
        KOKKOS_FUNCTION static T pathDistanceRoundoff(
                const Vector3D_t& position, const SegmentGeometryView_t& bounds);
        /** @internal @brief Compute the unchanged squared point-to-box lower bound.
         * @param node Heap index. @param position Query point. @param bounds Heap boxes.
         * @return Rounded squared lower bound, used only with the conservative allowance. */
        KOKKOS_FUNCTION static T pathNodeDistance2(
                size_t node, const Vector3D_t& position, const SegmentGeometryView_t& bounds);
        /** @internal @brief Search the full hierarchy, falling back to exhaustive scoring on stack
         * exhaustion.
         * @param position Query point. @param seed Previous segment, used only as an upper bound.
         * @param ref Vertices. @param geometry Segment data. @param bounds Heap boxes.
         * @param leaves Power-of-two leaves. @param best Winning index, initialized by caller.
         * @param bestFraction Winning fraction. @param closest Winning Cartesian point. */
        KOKKOS_FUNCTION static void findClosestPathSegment(
                const Vector3D_t& position, size_t seed, const ReferenceView_t& ref,
                const SegmentGeometryView_t& geometry, const SegmentGeometryView_t& bounds,
                size_t leaves, size_t& best, T& bestFraction, Vector3D_t& closest);
        /** @internal @brief Project one particle using exact segment scoring and the existing
         * basis.
         * @param n Particle index. @param positions Cartesian source. @param invalid Invalid mask.
         * @param ref Vertices. @param geometry Segment data. @param bounds Heap boxes.
         * @param projected Output local coordinates. @param indices Output winning indices.
         * @param leaves Heap leaf count. @param periodicLength Post-migration seam-wrap length. */
        KOKKOS_FUNCTION static void projectPreparedParticle(
                size_t n, const VectorView_t& positions, const BooleanView_t& invalid,
                const ReferenceView_t& ref, const SegmentGeometryView_t& geometry,
                const SegmentGeometryView_t& bounds, const VectorView_t& projected,
                const SegmentIndexView_t& indices, size_t leaves, T periodicLength);

        /** @internal @brief Clear self fields and enter the common solve frame for active
         * containers.
         * @param context Activity flags. @param labToSolve Current laboratory-to-solve transform.
         */
        void enterParticleSolveFrames(
                const SpaceChargeSolveContext& context, const CoordinateSystemTrafo& labToSolve);
        /** @internal @brief Restore tracker frames, refreshing moments unless the tracker defers.
         * @param context Activity flags. @param labToSolve Transform used for entry.
         * @param refreshMoments Always true when restoring after an exception. */
        void restoreParticleSolveFrames(
                const SpaceChargeSolveContext& context, const CoordinateSystemTrafo& labToSolve,
                bool refreshMoments = true);
        /** @internal @brief Reject transverse ownership outside the fixed open FFT domain.
         * @param container Active container index. All ranks participate in the error reduction. */
        void checkMigrationDomain(size_t container);
        /** @internal @brief Save registered Cartesian coordinates and expose logical slab
         * coordinates.
         * @param container Active container with current projections; particle attributes migrate
         * together. */
        void saveMigrationPositions(size_t container);
        /** @internal @brief Restore Cartesian coordinates from the migrated registered attribute.
         * @param pc Container on its current ownership; invalidates its moment cache. */
        void restoreMigrationPositions(ParticleContainer& pc);
        /** @internal @brief Migrate one container through OPEN spatial infrastructure and restore
         * its BC.
         * @param container Active container index. Restoration also runs if particle migration
         * throws. */
        void migrateContainer(size_t container);

        /** @internal @brief Check charge and open-path coverage before
         * deposition.
         * @param container Active source index. @param usePreparedProjections Use the per-solve
         * caches. @param meansOnly Leave full moments dirty; scatter computes mean pz directly. */
        void validateDepositionContainer(size_t container, bool usePreparedProjections, bool meansOnly);
        /// Number-weighted global mean pz in the current solve frame; includes all stored particles.
        T meanLongitudinalMomentum(size_t container) const;
        /** @internal @brief Deposit one container and restore its charge-scaled time steps.
         * @tparam DiagnosticPolicy Stage observer. @param container Active source index.
         * @param diagnostic Observer. @param usePreparedProjections Use current projection caches.
         * @param meansOnly Use the direct global mean pz instead of full distribution moments.
         */
        template <typename DiagnosticPolicy>
        void scatterContainer(
                size_t container, DiagnosticPolicy diagnostic, bool usePreparedProjections,
                bool meansOnly);
        /** @internal @brief Fold serial periodic ghosts, scale charge and accumulate MPI source
         * halos.
         * @tparam DiagnosticPolicy Stage observer. @param context Current dt/activity.
         * @param diagnostic Receives normalized density after pre-halo line-density assembly. */
        template <typename DiagnosticPolicy>
        void finishChargeDeposition(
                const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic);
        /** @internal @brief Sum transverse planes and assemble the replicated ordered 1D profile.
         * Source ghosts stay separate until host accumulation; boundary ghosts are filled before
         * return. */
        void sumLineDensityProfile();
        /** @internal @brief Convert raw plane sums to line charge and centered longitudinal
         * differences.
         * @tparam DiagnosticPolicy Observer. @param diagnostic Receives area-scaled density and
         * gradient. */
        template <typename DiagnosticPolicy>
        void scaleLineDensityAndGradient(DiagnosticPolicy diagnostic);
        /** @internal @brief Run MPI-local batches while exposing each slice potential to
         * diagnostics.
         * @tparam OccupiedView Host mask type. @tparam DiagnosticPolicy Observer.
         * @param occupiedHost Mask of nonzero local planes. @param diagnostic Potential observer.
         */
        template <typename OccupiedView, typename DiagnosticPolicy>
        void solveDiagnosticSlabs(OccupiedView occupiedHost, DiagnosticPolicy diagnostic);
        /** @internal @brief Solve occupied one-rank slices with a shared serial workspace.
         * @tparam OccupiedView Host mask type. @tparam DiagnosticPolicy Observer.
         * @param occupiedHost Nonzero source mask; neutral nonzero sources remain active.
         * @param diagnostic Receives potentials including empty planes and halos. */
        template <typename OccupiedView, typename DiagnosticPolicy>
        void solveSerialSlices(OccupiedView occupiedHost, DiagnosticPolicy diagnostic);
        /** @internal @brief Fill z field ghosts and retain open transverse boundaries. */
        void setTransverseFieldGhosts();
        /** @internal @brief Return the unchanged local longitudinal closure coefficient.
         * @return Geometry factor divided by 4*pi*epsilon0; OPEN uses its fixed model constant. */
        T longitudinalCoupling() const;
        /** @internal @brief Gather one active container with the selected z interpolation.
         * @tparam ScatterLongitudinally Select linear z interpolation. @tparam DiagnosticPolicy
         * Observer.
         * @param container Source index. @param diagnostic Field observer.
         * @param usePreparedProjections Use current per-solve path caches.
         * @param meansOnly Reuse the prepared scalar mean, or compute it for a standalone gather. */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy>
        void gatherContainer(
                size_t container, DiagnosticPolicy diagnostic, bool usePreparedProjections,
                bool meansOnly);

        /** @internal @brief Recover the unchanged path basis and momentum components from a cached
         * winner.
         * @param n Particle index. @param p Cartesian momenta. @param geometry Segment data.
         * @param projectedR Cached local positions. @param projectedSegment Cached winners.
         * @param fsR Output position. @param fsP Output momentum. @param bUnit Output binormal.
         * @param nUnit Output fixed normal. @param tUnit Output tangent. */
        KOKKOS_FUNCTION static void loadPreparedCoordinates(
                size_t n, const VectorView_t& p, const SegmentGeometryView_t& geometry,
                const VectorView_t& projectedR, const SegmentIndexView_t& projectedSegment,
                Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit,
                Vector3D_t& tUnit);

    private:
        // The first physical step is retained, but its cache/migration setup cost is
        // reported separately from subsequent (steady-state) space-charge solves.
        bool firstPhysicalSolve_m = true;
        FFT2D5Config config_m;
        std::vector<ParticleContainer*> particles_m;
        /// Prepared with the projection cache; reused by scatter/gather in one complete solve.
        std::vector<T> preparedMeanPs_m;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<ReferencePath> referencePath_m;
        ReferenceView_t solveReferencePath_m;
        SegmentGeometryView_t segmentGeometry_m;
        SegmentGeometryView_t segmentBounds_m;
        std::vector<VectorView_t> projectedPositions_m;
        std::vector<SegmentIndexView_t> projectedSegments_m;
        std::unique_ptr<FFT2D5FieldStorage> fieldStorage_m;
        LineDensityView_t lineDensity_m;
        LineDensityView_t lineDensityGradient_m;
        static constexpr size_t LineDensityGhostCells    = 2;
        static constexpr size_t LineDensityFirstRealCell = 1;

        // Longitudinal field constants from
        // Space charge modules for PyHEADTAIL, A Oeftiger, SE Hegglin, 2016
        // https://proceedings.jacow.org/hb2016/papers/mopr025.pdf
        static constexpr T CircularPipeG0   = 0.67;
        static constexpr T ParallelPlatesG0 = 0.67;
        static constexpr T OpenG0           = 6.36;
    };

}  // namespace opalx::spacecharge
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.hpp"
#endif
