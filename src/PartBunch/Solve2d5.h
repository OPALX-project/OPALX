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

// This class provides a 2D5 solver for electromagnetic field simulations.
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

    enum class LongitudinalFieldMode { Cylindrical, Plates, Open, None };

    struct Solver {
        std::shared_ptr<OpenSolver2D_t> solver_m{};
        std::shared_ptr<Field_t<2U>> rho_m{};
        std::shared_ptr<VField_t<T, 2U>> E_m{};
    };

    Solve2d5(
            PartBunch_t* partBunch, std::string solver, Field_t<3U>* rho, VField_t<T, 3U>* E,
            Field_t<3U>* phi, std::shared_ptr<BCHandler_t> bcHandler, const Vector<int, 3U>& nR,
            LongitudinalFieldMode longitudinalFieldMode, T pipeSizeX, T pipeSizeY, T beamRadius,
            bool closedRing, bool scatterChargeLongitudinally, const std::string& refPathFileName);

    void initSolver() override;
    void runSolver() override { doRunSolver(); }
    void runSolver(const bool force_skip_field_dump) override {
        if (!force_skip_field_dump) {
            doRunSolver();
        }
    }
    void orbitThreadersReady() override;
    void computeSelfFields(PartBunch_t& /*bunch*/) override { doRunSolver(); }

    // Algorithm steps, public for testability
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
    template <typename DiagnosticPolicy = NullDiagnostic>
    void doRunSolver(DiagnosticPolicy diagnostic = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void scatterToGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void solvePoissons(DiagnosticPolicy diagnostic = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void calculateLineDensity(DiagnosticPolicy diagnostic = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void gatherFromGrid(const PartBunch_t& bunch, DiagnosticPolicy diagnostic = {});
    template <typename DiagnosticPolicy>
    std::unique_ptr<DiagnosticPolicy> createDiagnostic(NullDiagnostic::Kind kind);

private:
    // Helper functions
    T loadReferencePath();
    template <bool ScatterLongitudinally, typename DiagnosticPolicy = NullDiagnostic>
    KOKKOS_FUNCTION static void doScatterToGrid(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
            int nghost, ippl::NDIndex<3U> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
            DiagnosticPolicy diagnostic);
    KOKKOS_FUNCTION static void convertToFrenetSerret(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit,
            Vector3D_t& tUnit);
    KOKKOS_FUNCTION static void boostToBeamFrame(T meanPs, Vector3D_t& fsP);
    template <bool ScatterLongitudinally>
    KOKKOS_FUNCTION static void scatterToRho(
            size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3U>& lDom, ScalarGridView3D_t rho, Vector3D_t origin);
    template <typename DiagnosticPolicy = NullDiagnostic>
    KOKKOS_FUNCTION static void doGatherFromGrid(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b,
            const BooleanView_t& invalid, Vector3D_t invDr, int nghost, ippl::NDIndex<3U> lDom,
            VectorGridView3D_t eField, Vector3D_t origin, T gBy4PiEpsilon0,
            LineDensityView_t lineDensityGradient, DiagnosticPolicy diagnostic);
    KOKKOS_FUNCTION static void gatherFromEField(
            size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3U>& lDom, VectorGridView3D_t eField, Vector3D_t origin);
    KOKKOS_FUNCTION static void unboostFromBeamFrame(
            size_t n, T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b);
    KOKKOS_FUNCTION static void convertFromFrenetSerret(
            size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit, const Vector3D_t& tUnit,
            const VectorView_t& e, const VectorView_t& b);
    KOKKOS_FUNCTION static Vector3D_t gather2D(
            VectorGridView3D_t eField, const ippl::Vector<T, 3U>& wlo,
            const ippl::Vector<T, 3U>& whi, int x, int y, int z);
    KOKKOS_FUNCTION static void scatter2D(
            ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
            int x, int y, int z, T charge);
    KOKKOS_FUNCTION static void scatter3D(
            ScalarGridView3D_t rho, const ippl::Vector<T, 3U>& wlo, const ippl::Vector<T, 3U>& whi,
            int x, int y, int z, T charge);
    template <typename ViewType>
    KOKKOS_FUNCTION static bool makeWeights(
            Vector3D_t fsR, Vector3D_t origin, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3U>& lDom, const ViewType& view, ippl::Vector<T, 3U>& whi,
            ippl::Vector<T, 3U>& wlo, ippl::Vector<int, 3U>& args);

public:
    // Test case API
    size_t getNumSlices() const { return twoDSolvers_m.size(); }
    Mesh2D_t* getSliceMesh() const { return sliceMesh_m.get(); }
    Layout2D_t* getSliceLayout() const { return sliceLayout_m.get(); }
    const ReferenceView_t& getReferencePath() const { return referencePath_m; }
    Field_t<3U>* getRho() const { return rho_m; }
    VField_t<T, 3U>* getEField() const { return E_m; }
    const LineDensityView_t& getLineDensity() const { return lineDensity_m; }
    const LineDensityView_t& getLineDensityGradient() const { return lineDensityGradient_m; }

private:
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
