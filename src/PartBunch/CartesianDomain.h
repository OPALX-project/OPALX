/**
 * @file CartesianDomain.h
 * @brief Owns the mesh and field-layout lifetime shared by particle storage and PIC fields.
 */

#ifndef OPALX_PART_BUNCH_CARTESIAN_DOMAIN_H
#define OPALX_PART_BUNCH_CARTESIAN_DOMAIN_H

#include "Manager/BaseManager.h"
#include "PartBunch/CartesianDomainConfig.h"
#include "PartBunch/ParticleContainerTypes.h"

#include <array>
#include <cstddef>

namespace opalx::spacecharge {

    /**
     * @brief Stable Cartesian mesh and mutable field layout owned by PartBunch.
     *
     * Particle layouts and PIC fields borrow mesh() and layout(). Their addresses remain stable
     * while geometry and decomposition are updated in place by solver-owned policy.
     *
     * @note PartBunch declares this domain before its particle containers. Reverse member
     * destruction therefore releases every borrowing particle layout before the domain.
     */
    template <typename T, unsigned Dim>
    class CartesianDomain final {
    public:
        using Mesh    = ippl::UniformCartesian<T, Dim>;
        using Layout  = ippl::FieldLayout<Dim>;
        using Vector  = ippl::Vector<T, Dim>;
        using Extents = std::array<std::size_t, Dim>;

        explicit CartesianDomain(const CartesianDomainConfig<T, Dim>& config);

        CartesianDomain(const CartesianDomain&)            = delete;
        CartesianDomain& operator=(const CartesianDomain&) = delete;
        CartesianDomain(CartesianDomain&&)                 = delete;
        CartesianDomain& operator=(CartesianDomain&&)      = delete;

        [[nodiscard]] Mesh& mesh() { return mesh_m; }
        [[nodiscard]] const Mesh& mesh() const { return mesh_m; }
        [[nodiscard]] Layout& layout() { return layout_m; }
        [[nodiscard]] const Layout& layout() const { return layout_m; }

        [[nodiscard]] Vector& spacing() { return spacing_m; }
        [[nodiscard]] const Vector& spacing() const { return spacing_m; }
        [[nodiscard]] Vector& lower() { return lower_m; }
        [[nodiscard]] const Vector& lower() const { return lower_m; }
        [[nodiscard]] Vector& upper() { return upper_m; }
        [[nodiscard]] const Vector& upper() const { return upper_m; }
        [[nodiscard]] const Vector& origin() const { return origin_m; }
        [[nodiscard]] const std::array<bool, Dim>& decomposition() const { return decomposition_m; }
        [[nodiscard]] bool periodic() const { return periodic_m; }
        [[nodiscard]] const ippl::NDIndex<Dim>& indexDomain() const { return indexDomain_m; }
        [[nodiscard]] Extents layoutExtents() const;

        bool rebuildGlobalLayoutInPlace(
                const Extents& extents, const std::array<bool, Dim>& decomposition);
        void setGeometry(Vector lower, Vector upper, Vector spacing, Vector origin);

    private:
        [[nodiscard]] static ippl::NDIndex<Dim> makeIndexDomain(const Extents& extents);
        [[nodiscard]] static Vector initialLength(const CartesianDomainConfig<T, Dim>& config);
        [[nodiscard]] static Vector initialSpacing(const CartesianDomainConfig<T, Dim>& config);
        [[nodiscard]] static Vector initialOrigin(const CartesianDomainConfig<T, Dim>& config);

        Extents extents_m;
        std::array<bool, Dim> decomposition_m;
        bool periodic_m;
        ippl::NDIndex<Dim> indexDomain_m;
        Vector spacing_m;
        Vector lower_m;
        Vector upper_m;
        Vector origin_m;
        Mesh mesh_m;
        Layout layout_m;
    };

    extern template class CartesianDomain<double, 3>;

}  // namespace opalx::spacecharge

#include "PartBunch/CartesianDomain.tpp"

#endif  // OPALX_PART_BUNCH_CARTESIAN_DOMAIN_H
