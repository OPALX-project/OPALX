/**
 * @file FFT2D5Algorithm.h
 * @brief Declares the independent reference-path FFT2D5 space-charge algorithm.
 */

#ifndef OPALX_SPACE_CHARGE_FFT2D5_ALGORITHM_H
#define OPALX_SPACE_CHARGE_FFT2D5_ALGORITHM_H

#include "PartBunch/ParticleContainer.hpp"
#include "SpaceCharge/FFT2D5/FFT2D5FieldStorage.h"
#include "SpaceCharge/FFT2D5/ReferencePath.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeFrameTransform.h"

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
     * atomically on the first solve() after orbit threading has written the design path.
     */
    class FFT2D5Algorithm final : public SpaceChargeAlgorithm {
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

        FFT2D5Algorithm(
                FFT2D5Config config, std::span<const ParticleFieldBinding3D> particleBindings);

        [[nodiscard]] SpaceChargeCapabilities capabilities() const override;
        void solve(SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) override;

        [[nodiscard]] bool initialized() const { return fieldStorage_m != nullptr; }
        [[nodiscard]] std::size_t sliceCount() const;
        [[nodiscard]] const ReferencePath::View& referencePathView() const;

    public:
        /** @brief Lazily construct the reference path, fields, and per-slice Poisson solvers. */
        void ensureInitialized();

        /**
         * @brief Execute the FFT2D5 stages in their compatibility order.
         *
         * Particles are mapped to Frenet coordinates and boosted during deposition. Transverse
         * slice fields and the longitudinal line-density field are then gathered, unboosted, and
         * transformed back to Cartesian coordinates.
         */
        void solveFields(SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics);

        /**
         * @brief Deposit every selected container into the shared 3D charge-density staging field.
         *
         * Charge is scattered as dt*Q and normalized by the global time step and cell volume.
         * Closed rings fold longitudinal ghost charge into the opposite real boundary slice.
         */
        void scatterToGrid(const SpaceChargeSolveContext& context);

        /**
         * @brief Solve one transverse Poisson problem per longitudinal slice.
         *
         * Each 2D IPPL field owns its own allocation, so charge is copied out of the 3D staging
         * field and the solved transverse field is copied back after each solve.
         */
        void solveSlicePoissonProblems(SpaceChargeDiagnostics& diagnostics);

        /** @brief Integrate transverse density and form the guarded longitudinal gradient. */
        void calculateLineDensity();

        /** @brief Gather fields to selected particles and restore Cartesian lab-frame values. */
        void gatherFromGrid(const SpaceChargeSolveContext& context);

        /**
         * @brief Map one particle to Frenet coordinates, boost its momentum, and deposit charge.
         *
         * Bilinear deposition assigns all charge to one slice. Trilinear deposition distributes
         * charge between adjacent slices. Both paths use atomic additions for particle collisions.
         */
        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void scatterParticle(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, double meanLongitudinalMomentum,
                const ScalarView& timeStepCharge, const BooleanView& invalid,
                Vector3 inverseSpacing, int ghostCells, ippl::NDIndex<3> localDomain,
                ScalarGridView3 chargeDensity, Vector3 origin);

        /** @brief Project particle position and momentum onto the nearest path-segment basis. */
        KOKKOS_FUNCTION static void convertToFrenet(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, Vector3& frenetPosition, Vector3& frenetMomentum,
                Vector3& binormal, Vector3& normal, Vector3& tangent);
        /** @brief Lorentz-transform longitudinal momentum into the beam frame. */
        KOKKOS_FUNCTION static void boostMomentum(
                double meanLongitudinalMomentum, Vector3& momentum);

        template <bool ScatterLongitudinally>
        KOKKOS_FUNCTION static void deposit(
                std::size_t index, Vector3 frenetPosition, const ScalarView& timeStepCharge,
                Vector3 inverseSpacing, int ghostCells, const ippl::NDIndex<3>& localDomain,
                ScalarGridView3 chargeDensity, Vector3 origin);

        /**
         * @brief Gather one particle's transverse and longitudinal fields.
         *
         * The result is unboosted and transformed from the Frenet basis back to Cartesian axes.
         */
        KOKKOS_FUNCTION static void gatherParticle(
                std::size_t index, const VectorView& position, const VectorView& momentum,
                const ReferenceView& reference, double beamGamma, double beamBeta,
                const VectorView& electric, const VectorView& magnetic, const BooleanView& invalid,
                Vector3 inverseSpacing, int ghostCells, ippl::NDIndex<3> localDomain,
                VectorGridView3 electricGrid, Vector3 origin, double geometryFactor,
                LineDensityView lineDensityGradient);
        /** @brief Bilinearly gather both neighboring transverse slices without z weighting. */
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
        /** @brief Bilinearly interpolate one transverse slice. */
        KOKKOS_FUNCTION static Vector3 gatherSlice(
                VectorGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z);
        /** @brief Atomically deposit charge bilinearly onto one longitudinal slice. */
        KOKKOS_FUNCTION static void deposit2d(
                ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z, double charge);
        /** @brief Atomically deposit charge trilinearly across adjacent slices. */
        KOKKOS_FUNCTION static void deposit3d(
                ScalarGridView3 field, const Vector3& lowerWeights, const Vector3& upperWeights,
                int x, int y, int z, double charge);

        [[nodiscard]] bool isSelected(
                const SpaceChargeSolveContext& context, std::size_t index) const;
        [[nodiscard]] double longitudinalGeometryFactor() const;

    private:
        FFT2D5Config config_m;
        std::vector<ParticleContainer*> particles_m;
        std::unique_ptr<ReferencePath> referencePath_m;
        std::unique_ptr<FFT2D5FieldStorage> fieldStorage_m;
        LineDensityView lineDensity_m;
        LineDensityView lineDensityGradient_m;

        static constexpr std::size_t LineDensityGhostCells    = 2;
        static constexpr std::size_t LineDensityFirstRealCell = 1;
        static constexpr double CircularPipeG0                = 0.67;
        static constexpr double ParallelPlatesG0              = 0.67;
        static constexpr double OpenG0                        = 6.36;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_FFT2D5_ALGORITHM_H
