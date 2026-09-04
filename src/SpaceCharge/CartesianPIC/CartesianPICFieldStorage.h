/**
 * @file CartesianPICFieldStorage.h
 * @brief Declares persistent field and Cartesian-domain storage for a PIC solve.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_H

#include <array>
#include <cstddef>
#include <string_view>

#include "Manager/BaseManager.h"
#include "PartBunch/CartesianDomain.h"

namespace opalx::spacecharge {

    /**
     * @brief Owns persistent PIC fields and borrows the PartBunch Cartesian domain.
     *
     * PartBunch owns the mesh and field-layout lifetime. This object owns only fields and
     * scratch bound to those stable objects and refreshes them after solver-directed mutations.
     */
    template <typename T, unsigned Dim>
    class CartesianPICFieldStorage final {
    public:
        using Mesh        = ippl::UniformCartesian<T, Dim>;
        using Layout      = ippl::FieldLayout<Dim>;
        using Vector      = ippl::Vector<T, Dim>;
        using Domain      = CartesianDomain<T, Dim>;
        using Extents     = typename Domain::Extents;
        using ScalarField = ippl::Field<T, Dim, Mesh, typename Mesh::DefaultCentering>;
        using VectorField =
                ippl::Field<ippl::Vector<T, Dim>, Dim, Mesh, typename Mesh::DefaultCentering>;

        explicit CartesianPICFieldStorage(Domain& domain);

        CartesianPICFieldStorage(const CartesianPICFieldStorage&)            = delete;
        CartesianPICFieldStorage& operator=(const CartesianPICFieldStorage&) = delete;
        CartesianPICFieldStorage(CartesianPICFieldStorage&&)                 = delete;
        CartesianPICFieldStorage& operator=(CartesianPICFieldStorage&&)      = delete;

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

        [[nodiscard]] Vector& spacing() { return domain_m.spacing(); }
        [[nodiscard]] const Vector& spacing() const { return domain_m.spacing(); }
        [[nodiscard]] Vector& lower() { return domain_m.lower(); }
        [[nodiscard]] const Vector& lower() const { return domain_m.lower(); }
        [[nodiscard]] Vector& upper() { return domain_m.upper(); }
        [[nodiscard]] const Vector& upper() const { return domain_m.upper(); }

        [[nodiscard]] Mesh& mesh() { return domain_m.mesh(); }
        [[nodiscard]] const Mesh& mesh() const { return domain_m.mesh(); }
        [[nodiscard]] Layout& layout() { return domain_m.layout(); }
        [[nodiscard]] const Layout& layout() const { return domain_m.layout(); }
        [[nodiscard]] bool isAllPeriodic() const { return domain_m.periodic(); }
        [[nodiscard]] Domain& domain() { return domain_m; }
        [[nodiscard]] const Domain& domain() const { return domain_m; }

        [[nodiscard]] Extents layoutExtents() const { return domain_m.layoutExtents(); }

        /** @brief Per-unit electric accumulation in mesh axes after Lorentz conversion. */
        [[nodiscard]] VectorField& accumulatedElectricField() { return accumulatedElectricField_m; }
        [[nodiscard]] const VectorField& accumulatedElectricField() const {
            return accumulatedElectricField_m;
        }
        /** @brief Per-unit magnetic accumulation in mesh axes after Lorentz conversion. */
        [[nodiscard]] VectorField& accumulatedMagneticField() { return accumulatedMagneticField_m; }
        [[nodiscard]] const VectorField& accumulatedMagneticField() const {
            return accumulatedMagneticField_m;
        }
        /** @brief Out-of-place global z mirror of backend E; ghost values are not physical. */
        [[nodiscard]] VectorField& flippedZSlabField() { return flippedZSlabField_m; }
        [[nodiscard]] const VectorField& flippedZSlabField() const { return flippedZSlabField_m; }

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
        Domain& domain_m;
        bool fieldsInitialized_m    = false;
        bool potentialInitialized_m = false;

        VectorField electricField_m;
        ScalarField chargeDensity_m;
        ScalarField potential_m;
        VectorField accumulatedElectricField_m;
        VectorField accumulatedMagneticField_m;
        VectorField flippedZSlabField_m;
    };

    extern template class CartesianPICFieldStorage<double, 3>;

}  // namespace opalx::spacecharge

#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.tpp"

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_FIELD_STORAGE_H
