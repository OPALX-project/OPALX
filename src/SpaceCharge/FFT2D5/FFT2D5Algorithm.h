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
     * @brief Original reference-path 2.5D algorithm with solver-owned persistent storage.
     *
     * Numerical routines and their no-op diagnostic policy are retained from master Solve2d5.
     * The adapter borrows particle containers, initializes slices after the design path exists,
     * and applies the common tracker-frame and output contract around the original calculation.
     * PIPEMODE affects only longitudinal fields; transverse slices use the open Poisson solver.
     */
    class FFT2D5Algorithm final : public SpaceChargeAlgorithm {
    public:
        using T                       = double;
        static constexpr unsigned Dim = 3;
        using ParticleContainer       = ::ParticleContainer<double, 3>;
        using LongitudinalFieldMode   = FFT2D5LongitudinalFieldMode;
        using OpenSolver2D_t          = ippl::FFTOpenPoissonSolver<VField_t<T, 2U>, Field_t<2U>>;
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
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void doRunSolver(const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {});

        /**
         * @brief Deposits particle charge onto the 3D grid.
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect the charge
         *                           deposition process. Defaults to NullDiagnostic.
         * @param context Per-call activity and tracker state.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void scatterToGrid(
                const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {});

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
         *
         * @tparam DiagnosticPolicy Diagnostic policy used to inspect the field
         *                           gathering process. Defaults to NullDiagnostic.
         * @param context Per-call activity and tracker state.
         * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        void gatherFromGrid(
                const SpaceChargeSolveContext& context, DiagnosticPolicy diagnostic = {});

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

        // CUDA requires enclosing functions of device lambdas to be public.

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
         * @param nghost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param origin Physical origin of the charge-density grid.
         * @param diagnostic Diagnostic policy used to record intermediate results.
         */
        template <bool ScatterLongitudinally, typename DiagnosticPolicy = NullDiagnostic>
        KOKKOS_FUNCTION static void doScatterToGrid(
                size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
                T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
                int nghost, ippl::NDIndex<3U> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
                DiagnosticPolicy diagnostic);

        /**
         * @brief Kokkos function that converts from lab to Frenet-Serret coordinates.
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
                Vector3D_t& tUnit);

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
         * @param nghost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param rho Charge-density grid onto which particle charge is deposited.
         * @param origin Physical origin of the charge-density grid.
         */
        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void scatterToRho(
                size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, int nghost,
                const ippl::NDIndex<3U>& lDom, ScalarGridView3D_t rho, Vector3D_t origin);

        /**
         * @brief Kokkos kernel that gathers the E and B fields for a single particle.
         * See implementation for more details.
         *
         * @tparam DiagnosticPolicy Diagnostic policy.
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
         * @param nghost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param eField Electric field grid from the Poisson solver.
         * @param origin Physical origin of the charge-density grid.
         * @param gBy4PiEpsilon0 Constant for the longitudinal field estimate.
         * @param lineDensityGradient The longitudinal line charge density gradient.
         * @param diagnostic Diagnostic policy used to record intermediate results.
         */
        template <typename DiagnosticPolicy = NullDiagnostic>
        KOKKOS_FUNCTION static void doGatherFromGrid(
                size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
                T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b,
                const BooleanView_t& invalid, Vector3D_t invDr, int nghost, ippl::NDIndex<3U> lDom,
                VectorGridView3D_t eField, Vector3D_t origin, T gBy4PiEpsilon0,
                LineDensityView_t lineDensityGradient, DiagnosticPolicy diagnostic);

        /**
         * @brief Kokkos function that gathers the boosted E field for a particle from
         *  the electric field grid.
         *
         * @param n The particle number.
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param e Particle E fields.
         * @param invDr Inverse grid spacing in each spatial dimension.
         * @param nghost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param eField Electric field grid from the Poisson solver.
         * @param origin Physical origin of the charge-density grid.
         */
        KOKKOS_FUNCTION static void gatherFromEField(
                size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, int nghost,
                const ippl::NDIndex<3U>& lDom, VectorGridView3D_t eField, Vector3D_t origin);

        /**
         * @brief Kokkos function that unboosts a particle's e efield from the beam
         * reference frame into th elab frame, modifying the e field and creating the b field.
         *
         * @param n The particle number.
         * @param beamGamma The relativistic gamma for the beam.
         * @param beamBeta The relativistic beta for the beam.
         * @param e Particle E fields.
         * @param b Particle B fields.
         */
        KOKKOS_FUNCTION static void unboostFromBeamFrame(
                size_t n, T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b);

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
                VectorGridView3D_t eField, const ippl::Vector<T, 3U>& wlo,
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
                ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z, T charge);

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
                ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo,
                const ippl::Vector<T, 3U>& whi, int x, int y, int z, T charge);

        /**
         * @brief Kokkos function that generates the weights for scatter/gather.
         *
         * @param fsR The particle position in Frenet-Serret coordinates
         * @param origin Physical origin of the charge-density grid.
         * @param nghost Number of ghost cells surrounding the local grid domain.
         * @param lDom Local domain of the charge-density grid.
         * @param view The view being indexed.
         * @param whi Upper weighting factors.
         * @param wlo Lower weighting factors.
         * @param args The index coordinates of the cell containing fsR.
         */
        template <typename ViewType>
        KOKKOS_FUNCTION static bool makeWeights(
                Vector3D_t fsR, Vector3D_t origin, Vector3D_t invDr, int nghost,
                const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
                ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args);

        size_t getNumSlices() const { return fieldStorage_m->slices().size(); }
        Mesh2D_t* getSliceMesh() const {
            return &fieldStorage_m->slices().front().chargeDensity->get_mesh();
        }
        Layout2D_t* getSliceLayout() const {
            return &fieldStorage_m->slices().front().chargeDensity->getLayout();
        }
        const ReferenceView_t& getReferencePath() const { return referencePathView(); }
        Field_t<3>* getRho() const { return &fieldStorage_m->chargeDensity(); }
        VField_t<T, 3>* getEField() const { return &fieldStorage_m->electricField(); }
        const LineDensityView_t& getLineDensity() const { return lineDensity_m; }
        const LineDensityView_t& getLineDensityGradient() const { return lineDensityGradient_m; }

    private:
        FFT2D5Config config_m;
        std::vector<ParticleContainer*> particles_m;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<ReferencePath> referencePath_m;
        std::unique_ptr<FFT2D5FieldStorage> fieldStorage_m;
        LineDensityView_t lineDensity_m;
        LineDensityView_t lineDensityGradient_m;
        static constexpr size_t LineDensityGhostCells    = 2;
        static constexpr size_t LineDensityFirstRealCell = 1;
        // Space charge modules for PyHEADTAIL, A Oeftiger, SE Hegglin, 2016.
        static constexpr T CircularPipeG0   = 0.67;
        static constexpr T ParallelPlatesG0 = 0.67;
        static constexpr T OpenG0           = 6.36;
    };

}  // namespace opalx::spacecharge
#include "SpaceCharge/FFT2D5/FFT2D5Algorithm.hpp"
#endif
