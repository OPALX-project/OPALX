/**
 * @file Pic2d5Solver.h
 * @brief Declares the independent reference-path FFT2D5 self-field algorithm.
 */

#ifndef OPALX_SPACE_CHARGE_PIC2D5_PIC2D5_SOLVER_H
#define OPALX_SPACE_CHARGE_PIC2D5_PIC2D5_SOLVER_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/ParticleFrameGuard.h"
#include "SpaceCharge/Pic2d5/ReferencePath.h"
#include "SpaceCharge/Pic2d5/SliceWorkspace.h"
#include "SpaceCharge/SelfFieldAlgorithm.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <Kokkos_Core.hpp>

#include <memory>
#include <span>
#include <vector>

namespace opalx::spacecharge {

    /**
     * @brief Owns the full production FFT2D5 algorithm without legacy solver orchestration.
     *
     * Stable particle containers are borrowed for the run lifetime. Reference-path data, 3D
     * staging fields, slice fields, and every 2D Poisson backend are solver-owned and created
     * atomically on the first execute() after orbit threading has written the design path.
     */
    class Pic2d5Solver final : public SelfFieldAlgorithm {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using Vector3           = Vector_t<double, 3>;
        using VectorView        = ParticleAttrib<Vector3>::view_type;
        using ScalarView        = ParticleAttrib<double>::view_type;
        using BooleanView       = ParticleAttrib<bool>::view_type;
        using ReferenceView     = ReferencePath::View;
        using LineDensityView   = Kokkos::View<double*>;
        using ScalarGridView3   = Field_t<3>::view_type;
        using VectorGridView3   = VField_t<double, 3>::view_type;

        Pic2d5Solver(Pic2d5Config config, std::span<const ParticleFieldBinding3d> particleBindings);

        [[nodiscard]] SolverCapabilities capabilities() const override;
        void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) override;

        [[nodiscard]] bool initialized() const { return workspace_m != nullptr; }
        [[nodiscard]] std::size_t sliceCount() const;
        [[nodiscard]] const ReferencePath::View& referencePathView() const;

    public:
        void ensureInitialized();
        void run(SolveContext& context, SelfFieldDiagnostics& diagnostics);
        void scatterToGrid(const SolveContext& context);
        void solvePoissons(SelfFieldDiagnostics& diagnostics);
        void calculateLineDensity();
        void gatherFromGrid(const SolveContext& context);

        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void scatterParticle(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, double meanLongitudinalMomentum,
                const ScalarView& timeStepCharge, const BooleanView& invalid,
                Vector3 inverseSpacing, int ghostCells, ippl::NDIndex<3> localDomain,
                ScalarGridView3 chargeDensity, Vector3 origin);

        KOKKOS_FUNCTION static void convertToFrenet(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, Vector3& frenetPosition, Vector3& frenetMomentum,
                Vector3& binormal, Vector3& normal, Vector3& tangent);
        KOKKOS_FUNCTION static void boostMomentum(
                double meanLongitudinalMomentum, Vector3& momentum);

        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void deposit(
                std::size_t index, Vector3 frenetPosition, const ScalarView& timeStepCharge,
                Vector3 inverseSpacing, int ghostCells, const ippl::NDIndex<3>& localDomain,
                ScalarGridView3 chargeDensity, Vector3 origin);

        KOKKOS_FUNCTION static void gatherParticle(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, double beamGamma, double beamBeta,
                const VectorView& electric, const VectorView& magnetic, const BooleanView& invalid,
                Vector3 inverseSpacing, int ghostCells, ippl::NDIndex<3> localDomain,
                VectorGridView3 electricGrid, Vector3 origin, double geometryFactor,
                LineDensityView lineDensityGradient);
        KOKKOS_FUNCTION static void gatherTransverse(
                std::size_t index, Vector3 frenetPosition, const VectorView& electric,
                Vector3 inverseSpacing, int ghostCells, const ippl::NDIndex<3>& localDomain,
                VectorGridView3 electricGrid, Vector3 origin);
        KOKKOS_FUNCTION static void unboostFields(
                std::size_t index, double beamGamma, double beamBeta, const VectorView& electric,
                const VectorView& magnetic);
        KOKKOS_FUNCTION static void convertFieldsFromFrenet(
                std::size_t index, const Vector3& binormal, const Vector3& normal,
                const Vector3& tangent, const VectorView& electric, const VectorView& magnetic);

        template <typename View>
        KOKKOS_FUNCTION static bool makeWeights(
                Vector3 frenetPosition, Vector3 origin, Vector3 inverseSpacing, int ghostCells,
                const ippl::NDIndex<3>& localDomain, const View& view, Vector3& upperWeights,
                Vector3& lowerWeights, ippl::Vector<int, 3>& indices);
        KOKKOS_FUNCTION static Vector3 gatherSlice(
                VectorGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z);
        KOKKOS_FUNCTION static void deposit2d(
                ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z, double charge);
        KOKKOS_FUNCTION static void deposit3d(
                ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z, double charge);

        [[nodiscard]] bool selected(const SolveContext& context, std::size_t index) const;
        [[nodiscard]] double longitudinalGeometryFactor() const;

    private:
        Pic2d5Config config_m;
        std::vector<ParticleContainer*> particles_m;
        std::unique_ptr<ReferencePath> referencePath_m;
        std::unique_ptr<SliceWorkspace> workspace_m;
        LineDensityView lineDensity_m;
        LineDensityView lineDensityGradient_m;

        static constexpr std::size_t LineDensityGhostCells    = 2;
        static constexpr std::size_t LineDensityFirstRealCell = 1;
        static constexpr double CircularPipeG0                = 0.67;
        static constexpr double ParallelPlatesG0              = 0.67;
        static constexpr double OpenG0                        = 6.36;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC2D5_PIC2D5_SOLVER_H
