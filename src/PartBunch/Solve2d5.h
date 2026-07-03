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
class Solve2d5 : public BinnedFieldSolver<T, 3> {
public:
    using BCHandler_t        = BinnedFieldSolver<T, 3>::BCHandler_t;
    using Base               = BinnedFieldSolver<T, 3>;
    using OpenSolver2D_t     = ippl::FFTOpenPoissonSolver<VField_t<T, 2>, Field_t<2>>;
    using Mesh3D_t           = ippl::UniformCartesian<T, 3>;
    using Mesh2D_t           = ippl::UniformCartesian<T, 2>;
    using Vector2D_t         = Mesh2D_t::vector_type;
    using Vector3D_t         = Mesh3D_t::vector_type;
    using Layout2D_t         = ippl::FieldLayout<2>;
    using Point_t            = ippl::Vector<double, 3>;
    using VectorView_t       = ParticleAttrib<Vector3D_t>::view_type;
    using ScalarView_t       = ParticleAttrib<T>::view_type;
    using BooleanView_t      = ParticleAttrib<bool>::view_type;
    using PartBunch_t        = PartBunch<double, 3>;
    using ReferenceView_t    = Kokkos::View<Vector3D_t*>;
    using LineDensityView_t  = Kokkos::View<T*>;
    using ScalarGridView3D_t = Field<T, 3>::view_type;
    using VectorGridView3D_t = Field<Vector3D_t, 3>::view_type;
    using ScalarGridView2D_t = Field<T, 2>::view_type;
    using VectorGridView2D_t = Field<Vector2D_t, 2>::view_type;
    using FieldContainer_t   = FieldContainer<T, 3>;

    enum class LongitudinalFieldMode { Cylindrical, Plates, Open };

    struct Solver {
        std::shared_ptr<OpenSolver2D_t> solver_m{};
        std::shared_ptr<Field_t<2>> rho_m{};
        std::shared_ptr<VField_t<T, 2>> E_m{};
    };

    Solve2d5(
            PartBunch_t* partBunch, std::string solver, Field_t<3>* rho, VField_t<T, 3>* E,
            Field_t<3>* phi, std::shared_ptr<BCHandler_t> bcHandler, const Vector<int, 3>& nR,
            LongitudinalFieldMode longitudinalFieldMode, T pipeSizeX, T pipeSizeY, T beamRadius,
            bool closedRing);

    void initSolver() override;
    void runSolver() override { doRunSolver(); }
    void runSolver(bool /*force_skip_field_dump*/) override { doRunSolver(); }

    // Algorithm steps, public for testability
    class NullDiagnostic {
    public:
        enum class Kind {
            FrenetSerretScatter, BoostToBeam, ScatterCharge, ScatterChargeDensity,
            TotalDensity, LineDensity, LineDensityGradient, EField, FrenetSerretGather,
            GatherEField, Deboosted, LongitudinalField, LabFrameFields
        };
        KOKKOS_FUNCTION void frenetSerretScatter(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void boostToBeam(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void scatterCharge(const ScalarGridView3D_t&) const {}
        KOKKOS_FUNCTION void scatterChargeDensity(const ScalarGridView3D_t&) const {}
        KOKKOS_FUNCTION void eField(const VField_t<T, 3>::view_type&) const {}
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
                const PartBunch_t& /*partBunch*/, const Field_t<3>& /*rho*/,
                const LineDensityView_t& /*lineDensity*/,
                const LineDensityView_t& /*lineDensityGradient*/,
                const VField_t<T, 3>& /*eField*/) {}
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
    template <typename DiagnosticPolicy = NullDiagnostic>
    KOKKOS_FUNCTION static void doScatterToGrid(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            T meanPs, const ScalarView_t& dt, const BooleanView_t& invalid, Vector3D_t invDr,
            int nghost, ippl::NDIndex<3> lDom, ScalarGridView3D_t rho, Vector3D_t origin,
            DiagnosticPolicy diagnostic);
    KOKKOS_FUNCTION static void convertToFrenetSerret(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            Vector3D_t& fsR, Vector3D_t& fsP, Vector3D_t& bUnit, Vector3D_t& nUnit,
            Vector3D_t& tUnit);
    KOKKOS_FUNCTION static void boostToBeamFrame(T meanPs, Vector3D_t& fsP);
    KOKKOS_FUNCTION static void scatterToRho(
            size_t n, Vector3D_t fsR, const ScalarView_t& dt, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3>& lDom, ScalarGridView3D_t rho, Vector3D_t origin);
    template <typename DiagnosticPolicy = NullDiagnostic>
    KOKKOS_FUNCTION static void doGatherFromGrid(
            size_t n, const VectorView_t& r, const VectorView_t& p, const ReferenceView_t& ref,
            T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b,
            const BooleanView_t& invalid, Vector3D_t invDr, int nghost, ippl::NDIndex<3> lDom,
            VectorGridView3D_t eField, Vector3D_t origin, T gBy4PiEpsilon0,
            LineDensityView_t lineDensityGradient, DiagnosticPolicy diagnostic);
    KOKKOS_FUNCTION static void gatherFromEField(
            size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3>& lDom, VectorGridView3D_t eField, Vector3D_t origin);
    KOKKOS_FUNCTION static void unboostFromBeamFrame(
            size_t n, T beamGamma, T beamBeta, const VectorView_t& e, const VectorView_t& b);
    KOKKOS_FUNCTION static void convertFromFrenetSerret(
            size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit, const Vector3D_t& tUnit,
            const VectorView_t& e, const VectorView_t& b);

public:
    // Test case API
    size_t getNumSlices() const { return twoDSolvers_m.size(); }
    Mesh2D_t* getSliceMesh() const { return sliceMesh_m.get(); }
    Layout2D_t* getSliceLayout() const { return sliceLayout_m.get(); }
    const ReferenceView_t& getReferencePath() const { return referencePath_m; }
    Field_t<3>* getRho() const { return rho_m; }
    VField_t<T, 3>* getEField() const { return E_m; }
    const LineDensityView_t& getLineDensity() const { return lineDensity_m; }
    const LineDensityView_t& getLineDensityGradient() const { return lineDensityGradient_m; }

private:
    PartBunch_t* partBunch_m;
    std::vector<Solver> twoDSolvers_m;
    Field_t<3>* rho_m{};
    VField_t<T, 3>* E_m{};
    std::shared_ptr<Mesh2D_t> sliceMesh_m;
    std::shared_ptr<Layout2D_t> sliceLayout_m;
    ReferenceView_t referencePath_m;
    LineDensityView_t lineDensity_m;
    LineDensityView_t lineDensityGradient_m;
    Vector3D_t hr_m{1};
    Vector3D_t sizer_m{};
    Vector3D_t originr_m{};
    ippl::NDIndex<3> domain_m;

    // Configuration
    T beamRadius_m{1};
    LongitudinalFieldMode longitudinalFieldMode_m{LongitudinalFieldMode::Open};
    bool closedRing_m{false};
    Vector<unsigned int, 3> nR_m{10};

    // Constants
    static constexpr size_t LineDensityGhostCells    = 2;
    static constexpr size_t LineDensityFirstRealCell = 1;
};

#include "Solve2d5.tpp"

#endif  // OPALX_SOLVE2D5_H
