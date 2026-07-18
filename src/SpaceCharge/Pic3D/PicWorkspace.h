/**
 * @file PicWorkspace.h
 * @brief Declares persistent field and Cartesian-domain storage for a PIC solve.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_H
#define OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_H

#include <array>
#include <cstddef>
#include <memory>
#include <string_view>

#include "Manager/BaseManager.h"

namespace opalx::spacecharge {

    /**
     * @brief Owns stable mesh, layout, field, and scratch objects for Cartesian PIC.
     *
     * Particle layouts and IPPL fields borrow the addresses of @c mesh_m and @c layout_m.
     * Consequently this object is neither copyable nor movable. Layout changes mutate those
     * objects in place and refresh all persistent fields; no solve-pass allocation is required.
     */
    template <typename T, unsigned Dim>
    class PicWorkspace final {
    public:
        using Mesh        = ippl::UniformCartesian<T, Dim>;
        using Layout      = ippl::FieldLayout<Dim>;
        using Vector      = ippl::Vector<T, Dim>;
        using Extents     = std::array<std::size_t, Dim>;
        using ScalarField = ippl::Field<T, Dim, Mesh, typename Mesh::DefaultCentering>;
        using VectorField =
                ippl::Field<ippl::Vector<T, Dim>, Dim, Mesh, typename Mesh::DefaultCentering>;

        PicWorkspace(
                Vector spacing, Vector lower, Vector upper, std::array<bool, Dim> decomposition,
                ippl::NDIndex<Dim> domain, Vector origin, bool allPeriodic);

        PicWorkspace(const PicWorkspace&)            = delete;
        PicWorkspace& operator=(const PicWorkspace&) = delete;
        PicWorkspace(PicWorkspace&&)                 = delete;
        PicWorkspace& operator=(PicWorkspace&&)      = delete;

        /** @brief Backend output in the active mesh frame; valid on the current IPPL layout. */
        [[nodiscard]] VectorField& electricField() { return electricField_m; }
        [[nodiscard]] const VectorField& electricField() const { return electricField_m; }
        /** @brief Deposited charge in the active mesh frame, including IPPL ghost storage. */
        [[nodiscard]] ScalarField& chargeDensity() { return chargeDensity_m; }
        [[nodiscard]] const ScalarField& chargeDensity() const { return chargeDensity_m; }
        /** @brief Optional CG potential on the same mesh/layout; uninitialized for other backends.
         */
        [[nodiscard]] ScalarField& potential() { return potential_m; }
        [[nodiscard]] const ScalarField& potential() const { return potential_m; }

        [[nodiscard]] Vector& spacing() { return spacing_m; }
        [[nodiscard]] const Vector& spacing() const { return spacing_m; }
        [[nodiscard]] Vector& lower() { return lower_m; }
        [[nodiscard]] const Vector& lower() const { return lower_m; }
        [[nodiscard]] Vector& upper() { return upper_m; }
        [[nodiscard]] const Vector& upper() const { return upper_m; }

        [[nodiscard]] Mesh& mesh() { return mesh_m; }
        [[nodiscard]] const Mesh& mesh() const { return mesh_m; }
        [[nodiscard]] Layout& layout() { return layout_m; }
        [[nodiscard]] const Layout& layout() const { return layout_m; }
        [[nodiscard]] bool isAllPeriodic() const { return allPeriodic_m; }

        /** @brief Return the global extent represented by the mutable field layout. */
        [[nodiscard]] Extents layoutExtents() const;

        /**
         * @brief Rebuild the global field layout without replacing the borrowed layout object.
         *
         * The Cartesian mesh keeps its construction-time index domain. Only the global field
         * layout and persistent field allocations change here, matching the legacy longitudinal
         * image-domain resize behavior while keeping mesh/layout addresses stable.
         *
         * @param decomposition Parallel dimensions to apply on the next extent change.
         * @return true when the extent changed and the fields were refreshed.
         */
        bool rebuildGlobalLayoutInPlace(
                const Extents& extents, const std::array<bool, Dim>& decomposition);

        /** @brief Update the active physical bounds and mesh geometry in place. */
        void setGeometry(Vector lower, Vector upper, Vector spacing, Vector origin);

        /** @brief Per-unit electric accumulation in mesh axes after Lorentz conversion. */
        [[nodiscard]] VectorField& accumulatedElectricField() {
            return *accumulatedElectricField_m;
        }
        [[nodiscard]] const VectorField& accumulatedElectricField() const {
            return *accumulatedElectricField_m;
        }
        /** @brief Per-unit magnetic accumulation in mesh axes after Lorentz conversion. */
        [[nodiscard]] VectorField& accumulatedMagneticField() {
            return *accumulatedMagneticField_m;
        }
        [[nodiscard]] const VectorField& accumulatedMagneticField() const {
            return *accumulatedMagneticField_m;
        }
        /** @brief Out-of-place global z mirror of backend E; ghost values are not physical. */
        [[nodiscard]] VectorField& flippedZSlabField() { return *flippedZSlabField_m; }
        [[nodiscard]] const VectorField& flippedZSlabField() const { return *flippedZSlabField_m; }

        /** @brief Initialize every persistent field once the backend type is known. */
        void initializeFields(std::string_view solverType);

        /** @brief Refresh all initialized fields after an in-place layout mutation. */
        void updateFieldLayoutsAfterLayoutChange();

        /**
         * @brief Return persistent mirror scratch compatible with @p source.
         *
         * Normal solve passes reuse the preinitialized field. A configuration/layout mismatch is
         * rejected so the hot path never allocates or replaces scratch storage.
         */
        [[nodiscard]] VectorField& mirrorScratchFor(const VectorField& source);

    private:
        Vector spacing_m;
        Vector lower_m;
        Vector upper_m;
        std::array<bool, Dim> decomposition_m;
        Mesh mesh_m;
        Layout layout_m;
        bool allPeriodic_m          = false;
        bool fieldsInitialized_m    = false;
        bool potentialInitialized_m = false;

        VectorField electricField_m;
        ScalarField chargeDensity_m;
        ScalarField potential_m;
        std::shared_ptr<VectorField> accumulatedElectricField_m;
        std::shared_ptr<VectorField> accumulatedMagneticField_m;
        std::shared_ptr<VectorField> flippedZSlabField_m;
    };

    extern template class PicWorkspace<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/Pic3D/PicWorkspace.tpp"

#endif  // OPALX_SPACE_CHARGE_PIC3D_WORKSPACE_H
