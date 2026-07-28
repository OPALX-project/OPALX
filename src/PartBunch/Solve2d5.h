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

#ifndef OPALX_SOLVE2D5_H
#define OPALX_SOLVE2D5_H

#include "BinnedFieldSolver.h"

/**
 * @brief A 2D5 solver for electromagnetic field simulations.
 *
 * @details
 * The solver performs the following sequence of operations to account for space charge
 * in the electric and magnetic fields that particles feel.
 * 1. Translate the particle cartesian coordinates into Frenet-Serret coordinates and
 *    then boost into the beam reference frame.  The reference orbit can be provided from a file or
 *    acquired from the OrbitThreader.
 * 2. Do CiC scatter of charge into the 3D charge density grid.  This spreads the charge out both
 *    transversely and (optionally) longitudinally which should be a help in reducing noise.
 * 3. Copy the 3D charge density  grid into a stack of 2D grids along the S coordinate.
 * 4. Perform the 2D Poisson solve on each of the 2D grids in the stack to determine the
 *    X and Y electric field components.  Initially, the boundary conditions used will be open.
 * 5. Calculate the longitudinal S electric field component using the algorithm from PyHEADTAIL.
 *    This uses the local line density along the reference orbit and a geometry constant that
 *    encodes the effect of the beam pipe on the longitudinal component.
 * 6. Do CiC gather of fields for the particle positions.
 * 7. Translate from the beam reference frame into the lab frame in Frenet-Serret coodinates and
 *    then back into cartesian coordinates.
 *
 * @tparam T The type of the floating point values (float or double)
 */
template <typename T>
class Solve2d5 : public BinnedFieldSolver<T, 3U> {
public:
    using BCHandler_t        = BinnedFieldSolver<T, 3U>::BCHandler_t;
    using Base               = BinnedFieldSolver<T, 3U>;
    using OpenSolver2D_t     = ippl::FFTOpenPoissonSolver<VField_t<T, 2U>, Field_t<2U>>;
    using Mesh3D_t           = ippl::UniformCartesian<T, 3U>;
    using Mesh2D_t           = ippl::UniformCartesian<T, 2U>;
    using Vector2D_t         = Mesh2D_t::vector_type;
    using Vector3D_t         = Mesh3D_t::vector_type;
    using Layout2D_t         = ippl::FieldLayout<2U>;
    using Point_t            = ippl::Vector<double, 3U>;
    using VectorView_t       = ParticleAttrib<Vector3D_t>::view_type;
    using ScalarView_t       = ParticleAttrib<T>::view_type;
    using BooleanView_t      = ParticleAttrib<bool>::view_type;
    using PartBunch_t        = PartBunch<double, 3U>;
    using ReferenceView_t    = Kokkos::View<Vector3D_t*>;
    using LineDensityView_t  = Kokkos::View<T*>;
    using ScalarGridView3D_t = Field<T, 3U>::view_type;
    using VectorGridView3D_t = Field<Vector3D_t, 3U>::view_type;
    using ScalarGridView2D_t = Field<T, 2U>::view_type;
    using VectorGridView2D_t = Field<Vector2D_t, 2U>::view_type;
    using FieldContainer_t   = FieldContainer<T, 3U>;

    /**
     * @brief Specifies the model used to calculate the longitudinal space-charge field.
     *
     * The longitudinal electric field is calculated from the longitudinal
     * variation of the line charge density, with the selected mode determining
     * the geometric factor used in the calculation.
     */
    enum class LongitudinalFieldMode {
        /**
         * @brief Cylindrical beam pipe geometry. Uses min(pipeSizeX, pipeSizeY) as the radius.
         */
        Cylindrical,
        /**
         * @brief Parallel plates beam pipe geometry. Uses min(pipeSizeX, pipeSizeY) as the
         * plate separation (use for rectangular pipes).
         */
        Plates,
        /**
         * @brief Open-boundary geometry.  The pipe walls are not nearby.
         */
        Open,
        /**
         * @brief Do not calculate the longitudinal electric field.
         */
        None
    };

    /**
     * @brief Constructs a 2.5D space-charge solver.
     *
     * @param partBunch Pointer to the particle bunch on which the space-charge
     *                  calculation will be performed.
     * @param solver Name of this solver.
     * @param rho Pointer to the three-dimensional charge-density field.
     * @param E Pointer to the three-dimensional electric-field field.
     * @param phi Pointer to the three-dimensional electrostatic potential field.
     * @param bcHandler Shared pointer to the boundary-condition handler used
     *                  by the Poisson solver.
     * @param nR Number of grid cells in each of the three spatial dimensions.
     * @param longitudinalFieldMode Specifies the model used to calculate the
     *                              longitudinal space-charge field.
     * @param pipeSizeX Horizontal beam-pipe dimension, defines the horizontal domain extent
     *                  and used when calculating the longitudinal fields.
     * @param pipeSizeY Vertical beam-pipe dimension, defines the horizontal domain extent
     *                  and used when calculating the longitudinal fields.
     * @param beamRadius Beam radius used when calculating the longitudinal
     *                   space-charge field.
     * @param closedRing Indicates whether the simulation represents a closed
     *                   ring. If true, the longitudinal coordinate is treated
     *                   as periodic.
     * @param scatterChargeLongitudinally If true, particle charge is scattered
     *                                    between adjacent longitudinal slices.
     *                                    If false, the full particle charge is
     *                                    deposited onto a single longitudinal
     *                                    slice.
     * @param refPathFileName Name of the file containing the reference path
     *                        used to transform particle coordinates into the
     *                        local Frenet-Serret coordinate system.
     */
    Solve2d5(
            PartBunch_t* partBunch, std::string solver, Field_t<3U>* rho, VField_t<T, 3U>* E,
            Field_t<3U>* phi, std::shared_ptr<BCHandler_t> bcHandler, const Vector<int, 3U>& nR,
            LongitudinalFieldMode longitudinalFieldMode, T pipeSizeX, T pipeSizeY, T beamRadius,
            bool closedRing, bool scatterChargeLongitudinally, const std::string& refPathFileName);

    /**
     * @brief Initialise the solver override, currently does nothing
     * @see orbitThreadersReady()
     */
    void initSolver() override {}

    /**
     * @brief Execute the solver using its attached particle bunch.
     * At the moment I think this interface is unused.
     */
    void runSolver() override { doRunSolver(); }

    /**
     * @brief Execute the solver using its attached particle bunch.
     *
     * @param force_skip_field_dump Set to true if this is an initialising call
     */
    void runSolver(const bool force_skip_field_dump) override {
        if (!force_skip_field_dump) {
            doRunSolver();
        }
    }

    /**
     * @brief Called by the parallel tracker to indicate that the orbit threaders have completed.
     *        The solver can now be initialised as the reference path will now exist.
     */
    void orbitThreadersReady() override;

    /**
     * @brief Called by the parallel tracker to run the solver.
     *
     * @param bunch The particle bunch, parameter unused as it is the same bunch already known.
     */
    void computeSelfFields(PartBunch_t& /*bunch*/) override { doRunSolver(); }

    /**
     * @brief Default no-op diagnostic interface used by the 2.5D space-charge solver.
     *
     * @details
     * This class defines the diagnostic interface used internally by the
     * 2.5D space-charge solver. All methods are no-ops, so the default
     * implementation introduces no diagnostic processing or storage.
     *
     * The class is primarily intended to be used as the default diagnostic
     * type when no diagnostic information is required. Unit tests may derive
     * a diagnostic class from NullDiagnostic and override selected methods
     * to inspect intermediate quantities produced during the space-charge
     * calculation.
     *
     * The solver's step methods are templated on a diagnostic type. The
     * diagnostic type must provide an interface compatible with this class.
     *
     * The diagnostic callbacks correspond to key stages of the solver,
     * including coordinate transformations, charge deposition, charge-density
     * calculation, field calculation, field gathering, and transformation
     * back to the laboratory frame.
     *
     * @note This class is an implementation detail of the 2.5D space-charge
     *       solver and its associated unit tests. It is not intended to be
     *       used as a general-purpose diagnostic interface by external code.
     */
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
            LabFrameFields
        };
        KOKKOS_FUNCTION void frenetSerretScatter(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void boostToBeam(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void scatterCharge(const ScalarGridView3D_t&) const {}
        KOKKOS_FUNCTION void scatterChargeDensity(const ScalarGridView3D_t&) const {}
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
                const PartBunch_t& /*partBunch*/, const Field_t<3U>& /*rho*/,
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
    void doRunSolver(DiagnosticPolicy diagnostic = {});

    /**
     * @brief Deposits particle charge onto the 3D grid.
     *
     * @tparam DiagnosticPolicy Diagnostic policy used to inspect the charge
     *                           deposition process. Defaults to NullDiagnostic.
     * @param bunch Particle bunch from which charge is deposited.
     * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
     */
    template <typename DiagnosticPolicy = NullDiagnostic>
    void scatterToGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic = {});

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
     * @param bunch Particle bunch to which the electric field is gathered.
     * @param diagnostic Diagnostic policy instance used for diagnostic callbacks.
     */
    template <typename DiagnosticPolicy = NullDiagnostic>
    void gatherFromGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic = {});

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

private:
    /**
     * @brief Loads the reference path from the file specified in the constructor.
     *
     * @return The total length of the path
     */
    T loadReferencePath();

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
            ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
            int x, int y, int z, T charge);

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
            ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
            int x, int y, int z, T charge);

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

public:
    /**
     * @return Returns the number of slices in the Frenet-Serret domain.
     * For use only by test cases.
     */
    size_t getNumSlices() const { return twoDSolvers_m.size(); }

    /**
     * @return Returns the mesh object for the Frenet-Serret domain.
     * For use only by test cases.
     */
    Mesh2D_t* getSliceMesh() const { return sliceMesh_m.get(); }

    /**
     * @return Returns the slice layout object for the Frenet-Serret domain.
     * For use only by test cases.
     */
    Layout2D_t* getSliceLayout() const { return sliceLayout_m.get(); }

    /**
     * @return Returns the reference path view for the Frenet-Serret domain.
     * For use only by test cases.
     */
    const ReferenceView_t& getReferencePath() const { return referencePath_m; }

    /**
     * @return Returns the 3D charge density view for the Frenet-Serret domain.
     * For use only by test cases.
     */
    Field_t<3U>* getRho() const { return rho_m; }

    /**
     * @return Returns the 3D electric field view for the Frenet-Serret domain.
     * For use only by test cases.
     */
    VField_t<T, 3U>* getEField() const { return E_m; }

    /**
     * @return Returns the logntudinal charge density view for the Frenet-Serret domain.
     * For use only by test cases.
     */
    const LineDensityView_t& getLineDensity() const { return lineDensity_m; }

    /**
     * @return Returns the logntudinal charge density gradient view for the Frenet-Serret domain.
     * For use only by test cases.
     */
    const LineDensityView_t& getLineDensityGradient() const { return lineDensityGradient_m; }

private:
    /**
     * @brief One of these is created for each slice of the Frenet-Serret domain.
     * Contains the solver object and the 2D charge density and electric field views.
     */
    struct Solver {
        std::shared_ptr<OpenSolver2D_t> solver_m{};
        std::shared_ptr<Field_t<2U>> rho_m{};
        std::shared_ptr<VField_t<T, 2U>> E_m{};
    };

    PartBunch_t* partBunch_m;
    std::vector<Solver> twoDSolvers_m;
    Field_t<3U>* rho_m{};
    VField_t<T, 3U>* E_m{};
    std::shared_ptr<Mesh2D_t> sliceMesh_m;
    std::shared_ptr<Layout2D_t> sliceLayout_m;
    ReferenceView_t referencePath_m;
    LineDensityView_t lineDensity_m;
    LineDensityView_t lineDensityGradient_m;
    Vector3D_t hr_m{1};
    Vector3D_t sizer_m{};
    Vector3D_t originr_m{};
    ippl::NDIndex<3U> domain_m;
    std::unique_ptr<ippl::ParameterList> solverParams_m;

    // Configuration
    T beamRadius_m{1};
    LongitudinalFieldMode longitudinalFieldMode_m{LongitudinalFieldMode::Open};
    bool closedRing_m{false};
    bool scatterChargeLongitudinally_m{true};
    Vector<unsigned int, 3U> nR_m{10};
    std::string referencePathFileName_m{};
    std::string solver_m;
    T pipeSizeX_m;
    T pipeSizeY_m;

    // Constants
    static constexpr size_t LineDensityGhostCells    = 2;
    static constexpr size_t LineDensityFirstRealCell = 1;
};

#include "Solve2d5.hpp"

#endif  // OPALX_SOLVE2D5_H
