//
// Helper opalx::detail::mirrorField
//
// Produce a spatially-mirrored copy of an ippl::Field along a single axis: the
// cell at global index (i0, i1, i2) in `dst` receives the value from
// (i0, i1, i2) with the coordinate on `axis` replaced by (N_axis - 1 - idx).
//
// GPU + MPI aware: packs and unpacks on device via Kokkos::parallel_for, and
// uses ippl::Comm->isend / ippl::Comm->recv which hand device pointers to the
// MPI implementation when CUDA-aware MPI is available. No host staging in the
// common path.
//
// Intended use: after a shifted-Green's-function Poisson solve, mirror the
// resulting potential / E-field across the cathode plane to obtain the image
// contribution. Sign handling for the vector components is left to the caller
// (e.g. BinnedFieldSolver::accumulateFieldToTemp flips the perpendicular E
// components after the spatial reflection).
//
#ifndef OPALX_FIELD_MIRROR_HPP
#define OPALX_FIELD_MIRROR_HPP

#include <cassert>
#include <vector>

#include "Ippl.h"

#include "Communicate/Communicator.h"
#include "Communicate/Tags.h"
#include "Field/FieldBufferOps.hpp"
#include "FieldLayout/FieldLayout.h"
#include "Index/NDIndex.h"

namespace opalx {
    namespace detail {

        // Spatially mirror `src` into `dst` along `axis`.
        //
        // Requirements:
        // - src and dst have identical layout, mesh, and ghost count.
        // - axis < Field::dim.
        // - Dim == 3 (only 3D supported for now).
        //
        // In-place (same Field for src and dst) is supported via a scratch copy.
        template <typename Field>
        void mirrorField(const Field& src, Field& dst, unsigned axis) {
            constexpr unsigned Dim = Field::dim;
            static_assert(Dim == 3, "mirrorField: only Dim == 3 is implemented");
            assert(axis < Dim);
            assert(&src.getLayout() == &dst.getLayout()
                   && "mirrorField: src and dst must share a layout");
            assert(src.getNghost() == dst.getNghost()
                   && "mirrorField: src and dst must have equal ghost counts");

            using view_type  = typename Field::view_type;
            using value_type = typename view_type::value_type;

            // In-place safety: read-before-write via a scratch copy of src.
            if (static_cast<const void*>(&src) == static_cast<const void*>(&dst)) {
                Field scratch = src.deepCopy();
                mirrorField(scratch, dst, axis);
                return;
            }

            const auto& layout   = src.getLayout();
            const auto& ldom_r   = layout.getLocalNDIndex();
            const auto& gdom     = layout.getDomain();
            const auto& hDomains = layout.getHostLocalDomains();
            const int nghost     = src.getNghost();
            const int N_glob     = static_cast<int>(gdom[axis].length());
            const int P          = ippl::Comm->size();
            const int r          = ippl::Comm->rank();

            auto srcView = src.getView();
            auto dstView = dst.getView();

            // Deterministic zero-init: only cells inside a mirror-intersect region
            // are overwritten below, so ghost cells and any residual interior stay 0.
            Kokkos::deep_copy(dstView, value_type{});

            // Reflect p's local domain on the mirror axis (other axes untouched).
            auto reflectOnAxis = [&](const ippl::NDIndex<Dim>& ldom_p) {
                ippl::NDIndex<Dim> out;
                for (unsigned d = 0; d < Dim; ++d) {
                    if (d == axis) {
                        const int k0_p = ldom_p[d].first();
                        const int k1_p = k0_p + static_cast<int>(ldom_p[d].length());
                        out[d]         = ippl::Index(N_glob - k1_p, N_glob - k0_p - 1);
                    } else {
                        out[d] = ldom_p[d];
                    }
                }
                return out;
            };

            std::vector<MPI_Request> requests;
            ippl::detail::FieldBufferData<value_type> fd_send;
            ippl::detail::FieldBufferData<value_type> fd_recv;
            const int base_tag = ippl::mpi::tag::FIELD_MIRROR;

            const bool flipX = (axis == 0);
            const bool flipY = (axis == 1);
            const bool flipZ = (axis == 2);

            // Phase 1: post all non-blocking sends (avoids deadlock with phase-2 recvs).
            for (int p = 0; p < P; ++p) {
                if (p == r) {
                    continue;
                }
                ippl::NDIndex<Dim> ldom_p_mirror = reflectOnAxis(hDomains(p));
                if (!ldom_r.touches(ldom_p_mirror)) {
                    continue;
                }
                ippl::NDIndex<Dim> intersect = ldom_r.intersect(ldom_p_mirror);

                // Pack r's src data at `intersect` (r's src global coords) and send to p.
                // p will unpack the buffer with the mirror-axis flip flag, placing
                // the data at the reflected dst range.
                ippl::detail::solver_send_field<value_type, value_type, Dim>(
                        base_tag, 0, p, intersect, ldom_r, nghost, srcView, fd_send, requests);
            }

            // Phase 2: blocking recvs + unpack-with-flip on the mirror axis.
            for (int p = 0; p < P; ++p) {
                if (p == r) {
                    continue;
                }
                ippl::NDIndex<Dim> ldom_p_mirror = reflectOnAxis(hDomains(p));
                if (!ldom_r.touches(ldom_p_mirror)) {
                    continue;
                }
                ippl::NDIndex<Dim> intersect = ldom_r.intersect(ldom_p_mirror);

                ippl::detail::solver_recv<value_type, value_type, Dim>(
                        base_tag, 0, p, intersect, ldom_r, nghost, dstView, fd_recv, flipX, flipY,
                        flipZ);
            }

            // Self overlap: r holds both src and dst cells that participate in its
            // own mirror. Common (and only populated) case on single-rank runs.
            {
                ippl::NDIndex<Dim> ldom_r_mirror = reflectOnAxis(ldom_r);
                if (ldom_r.touches(ldom_r_mirror)) {
                    ippl::NDIndex<Dim> intersect = ldom_r.intersect(ldom_r_mirror);

                    const int first0 = intersect[0].first() + nghost - ldom_r[0].first();
                    const int first1 = intersect[1].first() + nghost - ldom_r[1].first();
                    const int first2 = intersect[2].first() + nghost - ldom_r[2].first();
                    const int last0  = intersect[0].last() + nghost - ldom_r[0].first() + 1;
                    const int last1  = intersect[1].last() + nghost - ldom_r[1].first() + 1;
                    const int last2  = intersect[2].last() + nghost - ldom_r[2].first() + 1;

                    // Map dst local index to src local index on the mirror axis:
                    //   dst_global[a] = dst_local[a] - nghost + ldom_r[a].first()
                    //   src_global[a] = N_glob - 1 - dst_global[a]
                    //   src_local[a]  = src_global[a] - ldom_r[a].first() + nghost
                    //                 = N_glob - 1 - dst_local[a] + 2 * (nghost -
                    //                 ldom_r[a].first())
                    const int k_shift = (N_glob - 1) + 2 * (nghost - ldom_r[axis].first());
                    const int axis_c  = static_cast<int>(axis);

                    using mdrange_type = Kokkos::MDRangePolicy<Kokkos::Rank<3>>;
                    Kokkos::parallel_for(
                            "opalx::mirrorField::self",
                            mdrange_type({first0, first1, first2}, {last0, last1, last2}),
                            KOKKOS_LAMBDA(const size_t i, const size_t j, const size_t k) {
                                size_t idx[3]    = {i, j, k};
                                idx[axis_c]      = static_cast<size_t>(k_shift) - idx[axis_c];
                                dstView(i, j, k) = srcView(idx[0], idx[1], idx[2]);
                            });
                    Kokkos::fence();
                }
            }

            if (!requests.empty()) {
                MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
            }
            ippl::Comm->freeAllBuffers();
        }

    }  // namespace detail
}  // namespace opalx

#endif  // OPALX_FIELD_MIRROR_HPP
