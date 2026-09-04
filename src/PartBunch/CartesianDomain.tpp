/**
 * @file CartesianDomain.tpp
 * @brief Implements stable Cartesian domain ownership.
 */

#ifndef OPALX_PART_BUNCH_CARTESIAN_DOMAIN_TPP
#define OPALX_PART_BUNCH_CARTESIAN_DOMAIN_TPP

#include "Utilities/OpalException.h"

#include <Kokkos_Core.hpp>

#include <limits>

namespace opalx::spacecharge {

    template <typename T, unsigned Dim>
    CartesianDomain<T, Dim>::CartesianDomain(const ParticleStorageConfig<T, Dim>& config)
        : extents_m(config.meshSize),
          decomposition_m(config.decomposition),
          periodic_m(config.periodicParticleBoundary),
          indexDomain_m(makeIndexDomain(config.meshSize)),
          spacing_m(initialSpacing(config)),
          lower_m(initialOrigin(config)),
          upper_m(lower_m + initialLength(config)),
          origin_m(lower_m),
          mesh_m(indexDomain_m, spacing_m, origin_m),
          layout_m(MPI_COMM_WORLD, indexDomain_m, decomposition_m, periodic_m) {
        if (config.boundingBoxIncreasePercent < T(0)) {
            throw OpalException(
                    "CartesianDomain::CartesianDomain",
                    "The bounding-box increase must not be negative.");
        }
    }

    template <typename T, unsigned Dim>
    typename CartesianDomain<T, Dim>::Extents CartesianDomain<T, Dim>::layoutExtents() const {
        Extents extents{};
        const auto& domain = layout_m.getDomain();
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            extents[dimension] = static_cast<std::size_t>(domain[dimension].length());
        }
        return extents;
    }

    template <typename T, unsigned Dim>
    bool CartesianDomain<T, Dim>::rebuildGlobalLayoutInPlace(
            const Extents& extents, const std::array<bool, Dim>& decomposition) {
        decomposition_m = decomposition;
        if (layoutExtents() == extents) {
            return false;
        }
        // Every particle layout, field, and FFT plan borrows this stable layout object. Device
        // work using the old decomposition must finish before initialize() mutates it in place.
        Kokkos::fence();
        const ippl::NDIndex<Dim> domain = makeIndexDomain(extents);
        layout_m.initialize(domain, decomposition_m, periodic_m);
        extents_m = extents;
        return true;
    }

    template <typename T, unsigned Dim>
    void CartesianDomain<T, Dim>::setGeometry(
            Vector lower, Vector upper, Vector spacing, Vector origin) {
        lower_m   = lower;
        upper_m   = upper;
        spacing_m = spacing;
        origin_m  = origin;
        mesh_m.setMeshSpacing(spacing_m);
        mesh_m.setOrigin(origin_m);
    }

    template <typename T, unsigned Dim>
    ippl::NDIndex<Dim> CartesianDomain<T, Dim>::makeIndexDomain(const Extents& extents) {
        ippl::NDIndex<Dim> domain;
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            if (extents[dimension] == 0
                || extents[dimension] > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw OpalException(
                        "CartesianDomain::makeIndexDomain",
                        "Every Cartesian extent must fit in a positive IPPL Index.");
            }
            domain[dimension] = ippl::Index(static_cast<int>(extents[dimension]));
        }
        return domain;
    }

    template <typename T, unsigned Dim>
    typename CartesianDomain<T, Dim>::Vector CartesianDomain<T, Dim>::initialLength(
            const ParticleStorageConfig<T, Dim>& config) {
        Vector length(T(6));
        if (config.layoutKind == ParticleLayoutKind::SpatialOverlap) {
            if (!(config.overlapCutoff > T(0))) {
                throw OpalException(
                        "CartesianDomain::initialLength",
                        "A spatial-overlap particle layout requires a positive cutoff.");
            }
            for (unsigned dimension = 0; dimension < Dim; ++dimension) {
                length[dimension] =
                        static_cast<T>(config.meshSize[dimension]) * config.overlapCutoff;
            }
        }
        return length;
    }

    template <typename T, unsigned Dim>
    typename CartesianDomain<T, Dim>::Vector CartesianDomain<T, Dim>::initialSpacing(
            const ParticleStorageConfig<T, Dim>& config) {
        const Vector length = initialLength(config);
        Vector spacing(T(0));
        for (unsigned dimension = 0; dimension < Dim; ++dimension) {
            spacing[dimension] = length[dimension] / static_cast<T>(config.meshSize[dimension]);
        }
        return spacing;
    }

    template <typename T, unsigned Dim>
    typename CartesianDomain<T, Dim>::Vector CartesianDomain<T, Dim>::initialOrigin(
            const ParticleStorageConfig<T, Dim>& config) {
        return config.layoutKind == ParticleLayoutKind::SpatialOverlap
                       ? T(-0.5) * initialLength(config)
                       : Vector(T(-3));
    }

}  // namespace opalx::spacecharge

#endif  // OPALX_PART_BUNCH_CARTESIAN_DOMAIN_TPP
