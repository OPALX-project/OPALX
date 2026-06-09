/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

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
    using VField2D_t         = Field<Vector_t<T, 2>, 2>;

    struct Solver {
        std::shared_ptr<OpenSolver2D_t> solver_m{};
        std::shared_ptr<Field_t<2>> rho_m{};
        std::shared_ptr<VField_t<T, 2>> E_m{};
    };

    Solve2d5(
            std::string solver, Field_t<3>* rho, VField_t<T, 3>* E, Field_t<3>* phi,
            std::shared_ptr<BCHandler_t> bcHandler);

    void initSolver() override;

    void runSolver() override;
    void runSolver(bool /*force_skip_field_dump*/) override { Solve2d5::runSolver(); }

    // Algorithm steps, public for testability
    class NullDiagnostic {
    public:
        KOKKOS_FUNCTION void frenetSerret(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void boostToBeam(
                const size_t, const Vector3D_t&, const Vector3D_t&, bool) const {}
        KOKKOS_FUNCTION void scatterCharge(const ScalarGridView3D_t&) const {}
        KOKKOS_FUNCTION void eField(const VField_t<T, 3>::view_type&) const {}
        KOKKOS_FUNCTION void totalDensity(const LineDensityView_t&) const {}
        KOKKOS_FUNCTION void lineDensity(const LineDensityView_t&) const {}
        KOKKOS_FUNCTION void lineDensityGradient(const LineDensityView_t&) const {}
    };
    void loadReferencePath();
    template <typename DiagnosticPolicy = NullDiagnostic>
    void scatterToGrid(const PartBunch_t& bunch, DiagnosticPolicy = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void solvePoissons(DiagnosticPolicy = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void calculateLineDensity(DiagnosticPolicy = {});
    template <typename DiagnosticPolicy = NullDiagnostic>
    void gatherFromGrid(const PartBunch_t& bunch, DiagnosticPolicy = {});

private:
    // Helper functions
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
            T meanPs, const VectorView_t& e, const VectorView_t& b, const BooleanView_t& invalid,
            Vector3D_t invDr, int nghost, ippl::NDIndex<3> lDom, VectorGridView3D_t eField,
            Vector3D_t origin, DiagnosticPolicy diagnostic);
    KOKKOS_FUNCTION static void gatherFromEField(
            size_t n, Vector3D_t fsR, const VectorView_t& e, Vector3D_t invDr, int nghost,
            const ippl::NDIndex<3>& lDom, VectorGridView3D_t eField, Vector3D_t origin);
    KOKKOS_FUNCTION static void unboostFromBeamFrame(
            size_t n, T meanPs, VectorView_t& e, VectorView_t& b);
    KOKKOS_FUNCTION static void convertFromFrenetSerret(
            size_t n, const Vector3D_t& bUnit, const Vector3D_t& nUnit,
            const Vector3D_t& tUnit, VectorView_t& e, VectorView_t& b);

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
    std::vector<Solver> twoDSolvers_m;
    Field_t<3>* phi_m;
    Field_t<3>* rho_m;
    VField_t<T, 3>* E_m;
    std::shared_ptr<Mesh2D_t> sliceMesh_m;
    std::shared_ptr<Layout2D_t> sliceLayout_m;
    ReferenceView_t referencePath_m;
    LineDensityView_t lineDensity_m;
    LineDensityView_t lineDensityGradient_m;
};

#include "Solve2d5.tpp"

#endif  // OPALX_SOLVE2D5_H
