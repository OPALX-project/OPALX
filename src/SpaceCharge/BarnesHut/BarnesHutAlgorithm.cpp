/**
 * @file BarnesHutAlgorithm.cpp
 * @brief Implements the Barnes-Hut space-charge algorithm on top of IPPL's NBody module.
 *
 * This is the only OPALX translation unit that includes NBody/cstone/ryoanji headers.
 */

#include "SpaceCharge/BarnesHut/BarnesHutAlgorithm.h"

#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/RestFrame.h"
#include "SpaceCharge/SpaceChargeFrames.h"
#include "Utilities/OpalException.h"

#include "NBody/NBodyParticleContainer.hpp"
#include "NBody/NBodySolver.hpp"
#include "NBody/core/BHFieldLists.hpp"
#include "NBody/wrappers/NBodyKokkosView.hpp"
#include "ryoanji/interface/multipole_holder.cuh"
#include "ryoanji/nbody/heaviside_qpole.hpp"
#include "ryoanji/nbody/traversal_cpu.hpp"

#include "Physics/Physics.h"
#include "Utility/IpplTimings.h"

#include <Kokkos_Core.hpp>
#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <mpi.h>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        using size_type = ippl::detail::size_type;
        using Conserved = ippl::nbody::fields::AuxConserved;
        using Dependent = ippl::nbody::fields::StdDependent;
        /// With per-particle velocities the source rapidity travels with the particles and the
        /// second field channel (Cx, Cy) is an output.
        using RelConserved = ippl::nbody::fields::RelConserved;
        using RelDependent = ippl::nbody::fields::RelDependent;
        /// SFC key type of every precision policy's container.
        using BHKeyType =
                ippl::nbody::NBodyParticleContainer<ippl::nbody::DoublePrecision, 3>::KeyType;

        /// Attributes the copy-back rewrites: R, ID (ParticleBase) and dt, Phi, Bin, P, E, B,
        /// InvalidMask (ParticleContainer). Any other attribute would be lost on migration.
        constexpr unsigned kSupportedAttributeCount = 9;

        constexpr unsigned kMinimumBucketSize = 64;

        /// Below this many particles per rank (on average) cstone cannot give every rank a share
        /// of the octree; such bunches use ryoanji's direct sum instead (cost O(N^2), N small).
        constexpr unsigned long long kDirectSumParticlesPerRank = 64;

        /// Image charges of a grounded plane travel through the Barnes-Hut container as sources
        /// and are told apart from real particles by this bit in the 64-bit ID payload (OPALX IDs
        /// are non-negative).
        constexpr std::uint64_t kImageFlag = std::uint64_t(1) << 63;

        /**
         * @brief Rest-frame z of a lab-frame z: z' = reference + factor (z - reference) with
         * factor = gamma. The identity (factor 1) leaves z bit-for-bit unchanged.
         */
        struct ZStretch {
            double reference = 0.0;
            double factor    = 1.0;

            KOKKOS_INLINE_FUNCTION double apply(double z) const {
                return factor == 1.0 ? z : reference + factor * (z - reference);
            }
            KOKKOS_INLINE_FUNCTION double invert(double z) const {
                return factor == 1.0 ? z : reference + (z - reference) / factor;
            }
        };

        /// Wait for raw CUDA (cstone/ryoanji, default stream) and Kokkos work alike.
        void synchronizeDevice() {
            Kokkos::fence();
            ippl::nbody::syncExec();
        }

        BarnesHutAlgorithm::ParticleContainer& requirePrimary(
                std::span<BarnesHutAlgorithm::ParticleContainer* const> particles) {
            if (particles.empty() || particles.front() == nullptr) {
                throw OpalException(
                        "BarnesHutAlgorithm::BarnesHutAlgorithm",
                        "The primary particle container is not available.");
            }
            return *particles.front();
        }

        /** @brief Global axis-aligned bounds of the local R views in the current frame. */
        std::array<double, 6> globalBounds(BarnesHutAlgorithm::ParticleContainer& particles) {
            const auto R           = particles.R.getView();
            const size_type nLocal = particles.getLocalNum();
            using MinMax           = Kokkos::MinMax<double>;
            std::array<double, 6> bounds{};
            for (unsigned d = 0; d < 3; ++d) {
                typename MinMax::value_type extrema;
                Kokkos::parallel_reduce(
                        "BarnesHutAlgorithm::globalBounds", nLocal,
                        KOKKOS_LAMBDA(const size_type i, typename MinMax::value_type& update) {
                            const double v = R(i)[d];
                            update.min_val = v < update.min_val ? v : update.min_val;
                            update.max_val = v > update.max_val ? v : update.max_val;
                        },
                        MinMax(extrema));
                bounds[2 * d]     = extrema.min_val;
                bounds[2 * d + 1] = -extrema.max_val;
            }
            // One MIN reduction over (min, -max) pairs; empty ranks contribute +/-max().
            std::array<double, 6> global{};
            ippl::Comm->allreduce(bounds.data(), global.data(), 6, std::less<double>());
            for (unsigned d = 0; d < 3; ++d) {
                global[2 * d + 1] = -global[2 * d + 1];
            }
            return global;
        }

        /**
         * @brief Barnes-Hut domain for the axis-aligned bounds @p limits (fixed boundaries).
         *
         * The box is fitted to the sources instead of being a cube, so cstone's mixed-dimension
         * SFC keys give axes shorter than the longest one fewer key bits and the octree cells of
         * flat or long bunches stay close to cubic; a near-cubic box keeps all bits. Axes shorter
         * than 2^-(maxTreeLevel - 1) of the longest one, including the zero extent of a flat
         * bunch, are widened about their centre: the keys cannot resolve them anyway, and every
         * axis must keep at least one key bit.
         */
        cstone::Box<double> fittedBox(const std::array<double, 6>& limits) {
            constexpr unsigned maxLevel = cstone::maxTreeLevel<BHKeyType>{};
            const double longest =
                    std::max({limits[1] - limits[0], limits[3] - limits[2], limits[5] - limits[4]});
            const double shortest = std::ldexp(longest, -static_cast<int>(maxLevel - 1));
            std::array<double, 6> box = limits;
            for (unsigned d = 0; d < 3; ++d) {
                if (box[2 * d + 1] - box[2 * d] < shortest) {
                    const double centre = 0.5 * (box[2 * d] + box[2 * d + 1]);
                    box[2 * d]          = centre - 0.5 * shortest;
                    box[2 * d + 1]      = centre + 0.5 * shortest;
                }
            }
            using cstone::BoundaryType;
            return cstone::Box<double>(
                    box[0], box[1], box[2], box[3], box[4], box[5], BoundaryType::fixed,
                    BoundaryType::fixed, BoundaryType::fixed);
        }

        /// fittedBox() of the source bounds widened by @p boxMargin (BHBOXMARGIN) times the extent
        /// on every side: room for the bunch to grow before the box has to be resized.
        cstone::Box<double> domainBox(const std::array<double, 6>& bounds, double boxMargin) {
            std::array<double, 6> limits = bounds;
            for (unsigned d = 0; d < 3; ++d) {
                const double margin = boxMargin * (bounds[2 * d + 1] - bounds[2 * d]);
                limits[2 * d] -= margin;
                limits[2 * d + 1] += margin;
            }
            return fittedBox(limits);
        }

        size_type countMarkedForDeletion(BarnesHutAlgorithm::ParticleContainer& particles) {
            const auto invalid     = particles.InvalidMask.getView();
            const size_type nLocal = particles.getLocalNum();
            size_type marked       = 0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::countInvalid", nLocal,
                    KOKKOS_LAMBDA(const size_type i, size_type& count) {
                        count += invalid(i) ? 1 : 0;
                    },
                    marked);
            return marked;
        }

        /**
         * @brief Write the local OPALX particles into BH slots [offset, offset + nLocal) with
         * rest-frame z; with images, their mirrors (2 planeZ - z, opposite charge) follow in
         * [offset + nLocal, offset + 2 nLocal). @p planeZ is the rest-frame plane.
         */
        template <class BHContainer>
        void packParticles(
                BHContainer& bh, BarnesHutAlgorithm::ParticleContainer& particles,
                typename BHContainer::LocalIndex offset, size_type nLocal, double chargeScale,
                double softening,
                bool images, double planeZ, ZStretch stretch,
                std::optional<double> referenceRapidity = std::nullopt) {
            using ippl::nbody::getView;
            auto x  = getView<"Rx">(bh);
            auto y  = getView<"Ry">(bh);
            auto z  = getView<"Rz">(bh);
            auto h  = getView<"h">(bh);
            auto m  = getView<"charge">(bh);
            auto px = getView<"Px">(bh);
            auto py = getView<"Py">(bh);
            auto pz = getView<"Pz">(bh);
            auto dt = getView<"aux0">(bh);
            auto id = getView<"aux64">(bh);

            const auto R     = particles.R.getView();
            const auto P     = particles.P.getView();
            const auto dtOut = particles.dt.getView();
            const auto ID    = particles.ID.getView();
            const size_type off = offset;
            const std::uint64_t imageFlag = kImageFlag;  // captured by value in device code
            // Per-particle velocities: source rapidity relative to the frame (images move the
            // other way). The view is empty otherwise and never written.
            const bool rapidity = referenceRapidity.has_value();
            const double wRef   = referenceRapidity.value_or(0.0);
            auto w              = rapidity ? getView<"w">(bh) : decltype(getView<"w">(bh))();

            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::pack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const size_type j = off + i;
                        x(j)              = R(i)[0];
                        y(j)              = R(i)[1];
                        z(j)              = stretch.apply(R(i)[2]);
                        px(j)             = P(i)[0];
                        py(j)             = P(i)[1];
                        pz(j)             = P(i)[2];
                        m(j)              = chargeScale * dtOut(i);
                        h(j)              = softening;
                        dt(j)             = dtOut(i);
                        id(j)             = static_cast<std::uint64_t>(ID(i));
                        const double wi   = rapidity ? longitudinalRapidity(P(i)) : 0.0;
                        if (rapidity) {
                            w(j) = wi - wRef;
                        }
                        if (images) {
                            const size_type k = j + nLocal;
                            x(k)              = R(i)[0];
                            y(k)              = R(i)[1];
                            z(k)              = 2.0 * planeZ - z(j);
                            px(k)             = P(i)[0];
                            py(k)             = P(i)[1];
                            pz(k)             = -P(i)[2];
                            m(k)              = -chargeScale * dtOut(i);
                            h(k)              = softening;
                            dt(k)             = dtOut(i);
                            id(k)             = static_cast<std::uint64_t>(ID(i)) | imageFlag;
                            if (rapidity) {
                                w(k) = -wi - wRef;
                            }
                        }
                    });
            Kokkos::fence();
        }

        /** @brief Number of real (non-image) particles among the BH-owned slots. */
        template <class BHContainer>
        size_type countOwnedReal(BHContainer& bh) {
            const auto id         = ippl::nbody::getView<"aux64">(bh);
            const size_type first = bh.startIndex();
            const std::uint64_t imageFlag = kImageFlag;
            size_type real        = 0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::countReal", static_cast<size_type>(bh.getLocalNum()),
                    KOKKOS_LAMBDA(const size_type k, size_type& count) {
                        count += (id(first + k) & imageFlag) == 0 ? 1 : 0;
                    },
                    real);
            return real;
        }

        /// Field of all sources per owned slot, saved before an image-only traversal.
        using SlotField = Kokkos::View<double* [3]>;

        /**
         * @brief Overwrite the (already resized) OPALX particles with the real particles of the
         * BH-owned set, in slot order, z back in the lab frame; image charges are skipped.
         *
         * E receives the field of all sources. If @p allSources is not empty it holds that field
         * and the slot accelerations are the image charges' field alone, which goes to B as
         * scratch for the boost; otherwise B = 0.
         */
        template <class BHContainer>
        void copyBackParticles(
                BHContainer& bh, BarnesHutAlgorithm::ParticleContainer& particles, double sign,
                ZStretch stretch, const SlotField& allSources, bool secondChannel = false) {
            using ippl::nbody::getView;
            const auto x  = getView<"Rx">(bh);
            const auto y  = getView<"Ry">(bh);
            const auto z  = getView<"Rz">(bh);
            const auto px = getView<"Px">(bh);
            const auto py = getView<"Py">(bh);
            const auto pz = getView<"Pz">(bh);
            const auto dt = getView<"aux0">(bh);
            const auto id = getView<"aux64">(bh);
            const auto ax = getView<"Ex">(bh);
            const auto ay = getView<"Ey">(bh);
            const auto az = getView<"Ez">(bh);
            // Per-particle velocities: the second channel goes to B as scratch for the boost.
            using ChannelView = decltype(getView<"Cx">(bh));
            const auto cx     = secondChannel ? getView<"Cx">(bh) : ChannelView();
            const auto cy     = secondChannel ? getView<"Cy">(bh) : ChannelView();

            auto R       = particles.R.getView();
            auto P       = particles.P.getView();
            auto dtOut   = particles.dt.getView();
            auto ID      = particles.ID.getView();
            auto E       = particles.E.getView();
            auto B       = particles.B.getView();
            auto Phi     = particles.Phi.getView();
            auto Bin     = particles.Bin.getView();
            auto invalid = particles.InvalidMask.getView();

            using index_type      = typename std::decay_t<decltype(ID)>::value_type;
            const size_type first = bh.startIndex();
            const size_type n     = bh.getLocalNum();

            const std::uint64_t imageFlag = kImageFlag;
            const bool imageField         = allSources.extent(0) > 0;
            // Destination of every real particle among the owned slots.
            Kokkos::View<size_type*> destination("BarnesHutAlgorithm::copyBackDestination", n);
            Kokkos::parallel_scan(
                    "BarnesHutAlgorithm::compactReal", n,
                    KOKKOS_LAMBDA(const size_type k, size_type& offset, const bool final) {
                        const bool real = (id(first + k) & imageFlag) == 0;
                        if (final) {
                            destination(k) = real ? offset : n;
                        }
                        offset += real ? 1 : 0;
                    });

            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::copyBack", n, KOKKOS_LAMBDA(const size_type slot) {
                        const size_type k = destination(slot);
                        if (k == n) {
                            return;  // image charge
                        }
                        const size_type j = first + slot;
                        R(k)     = Vector_t<double, 3>(x(j), y(j), stretch.invert(z(j)));
                        P(k)     = Vector_t<double, 3>(px(j), py(j), pz(j));
                        dtOut(k) = dt(j);
                        ID(k)    = static_cast<index_type>(id(j));
                        const Vector_t<double, 3> slotField(
                                sign * ax(j), sign * ay(j), sign * az(j));
                        if (secondChannel) {
                            E(k) = slotField;
                            B(k) = Vector_t<double, 3>(sign * cx(j), sign * cy(j), 0.0);
                        } else if (imageField) {
                            E(k) = Vector_t<double, 3>(
                                    sign * allSources(slot, 0), sign * allSources(slot, 1),
                                    sign * allSources(slot, 2));
                            B(k) = slotField;
                        } else {
                            E(k) = slotField;
                            B(k) = Vector_t<double, 3>(0.0);
                        }
                        Phi(k)     = 0.0;
                        Bin(k)     = 0;
                        invalid(k) = false;
                    });
            Kokkos::fence();
        }

        /**
         * @brief Copy the slot accelerations (field of all sources) of the owned slots and set the
         * charge of the owned real particles to zero, leaving the image charges as sources.
         */
        template <class BHContainer>
        SlotField keepFieldAndClearRealCharges(BHContainer& bh) {
            using ippl::nbody::getView;
            const auto ax                 = getView<"Ex">(bh);
            const auto ay                 = getView<"Ey">(bh);
            const auto az                 = getView<"Ez">(bh);
            const auto id                 = getView<"aux64">(bh);
            auto m                        = getView<"charge">(bh);
            const size_type first         = bh.startIndex();
            const size_type n             = bh.getLocalNum();
            const std::uint64_t imageFlag = kImageFlag;
            SlotField allSources("BarnesHutAlgorithm::allSources", n);
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::keepAllSources", n, KOKKOS_LAMBDA(const size_type slot) {
                        const size_type j   = first + slot;
                        allSources(slot, 0) = ax(j);
                        allSources(slot, 1) = ay(j);
                        allSources(slot, 2) = az(j);
                        if ((id(j) & imageFlag) == 0) {
                            m(j) = 0.0;
                        }
                    });
            Kokkos::fence();
            return allSources;
        }

        /// field(i) = prefactor * a(first + i) for the targets of a direct sum.
        template <class FieldView>
        void unpackDirectSum(
                const FieldView& field, const Kokkos::View<double*>& ax,
                const Kokkos::View<double*>& ay, const Kokkos::View<double*>& az,
                std::size_t first, size_type nLocal, double prefactor) {
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::directUnpack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const std::size_t j = first + i;
                        field(i)            = Vector_t<double, 3>(
                                prefactor * ax(j), prefactor * ay(j), prefactor * az(j));
                    });
            Kokkos::fence();
        }

        /**
         * @brief Softened Coulomb field by ryoanji's direct sum, for bunches too small to split
         * into a balanced octree. Particles do not migrate.
         *
         * ryoanji::directSum is single-process: every rank gathers all (r, |q|, h) and evaluates
         * its own slice of targets against all sources with the same P2P kernel as the
         * Barnes-Hut traversal. The GPU variant has no prefactor, so G and the charge sign are
         * applied here, exactly as in the copy-back. Sources are placed with rest-frame z;
         * @p planeZ is the rest-frame plane. With @p imageField a second sum over the image
         * charges alone writes their field to B (scratch for the boost), as the copy-back does.
         */
        void directSumFields(
                BarnesHutAlgorithm::ParticleContainer& particles, double chargeScale,
                double softening, double prefactor, bool images, double planeZ, ZStretch stretch,
                bool imageField) {
            const size_type nLocal   = particles.getLocalNum();
            const size_type nSources = images ? 2 * nLocal : nLocal;
            const auto R             = particles.R.getView();
            const auto dtView        = particles.dt.getView();

            // Pack (x, y, z, |q|, h) per local particle, followed by the image charges (-|q|), and
            // gather them on every rank. Each rank's block starts with its real particles (the
            // targets).
            constexpr int kValues = 5;
            Kokkos::View<double*> local("BarnesHutAlgorithm::directLocal", kValues * nSources);
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::directPack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const double z         = stretch.apply(R(i)[2]);
                        local(kValues * i)     = R(i)[0];
                        local(kValues * i + 1) = R(i)[1];
                        local(kValues * i + 2) = z;
                        local(kValues * i + 3) = chargeScale * dtView(i);
                        local(kValues * i + 4) = softening;
                        if (images) {
                            const size_type k  = nLocal + i;
                            local(kValues * k)     = R(i)[0];
                            local(kValues * k + 1) = R(i)[1];
                            local(kValues * k + 2) = 2.0 * planeZ - z;
                            local(kValues * k + 3) = -chargeScale * dtView(i);
                            local(kValues * k + 4) = softening;
                        }
                    });
            const auto localHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), local);

            MPI_Comm comm   = ippl::Comm->getCommunicator();
            const int ranks = ippl::Comm->size();
            const int count = static_cast<int>(kValues * nSources);
            std::vector<int> counts(ranks), offsets(ranks, 0);
            MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            for (int r = 1; r < ranks; ++r) {
                offsets[r] = offsets[r - 1] + counts[r - 1];
            }
            std::vector<double> gathered(
                    static_cast<std::size_t>(offsets.back() + counts.back()));
            MPI_Allgatherv(
                    localHost.data(), count, MPI_DOUBLE, gathered.data(), counts.data(),
                    offsets.data(), MPI_DOUBLE, comm);

            // ryoanji wants one array per component; this rank's particles are a contiguous slice.
            const std::size_t nAll  = gathered.size() / kValues;
            const std::size_t first = static_cast<std::size_t>(offsets[ippl::Comm->rank()]) / kValues;
            std::array<std::vector<double>, kValues> columns;
            for (auto& column : columns) {
                column.resize(nAll);
            }
            for (std::size_t j = 0; j < nAll; ++j) {
                for (int c = 0; c < kValues; ++c) {
                    columns[c][j] = gathered[kValues * j + c];
                }
            }
            using HostUnmanaged = Kokkos::View<const double*, Kokkos::HostSpace,
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
            std::array<Kokkos::View<double*>, kValues> device;
            for (int c = 0; c < kValues; ++c) {
                device[c] = Kokkos::View<double*>("BarnesHutAlgorithm::directSource", nAll);
                Kokkos::deep_copy(device[c], HostUnmanaged(columns[c].data(), nAll));
            }
            // Outputs accumulate (+=); the views are zero-initialised.
            Kokkos::View<double*> potential("BarnesHutAlgorithm::directP", nAll);
            Kokkos::View<double*> ax("BarnesHutAlgorithm::directAx", nAll);
            Kokkos::View<double*> ay("BarnesHutAlgorithm::directAy", nAll);
            Kokkos::View<double*> az("BarnesHutAlgorithm::directAz", nAll);
            Kokkos::fence();

            // The collectives above include every rank; a rank without targets stops here. The
            // GPU kernel sizes its grid as (last - first - 1) / threads + 1, which wraps around
            // for an empty range and launches an enormous grid.
            if (nLocal == 0) {
                return;
            }
            using FieldView     = std::decay_t<decltype(particles.E.getView())>;
            const auto evaluate = [&](const FieldView& destination) {
                const ryoanji::Vec3<double> unusedBox{1.0, 1.0, 1.0};  // numShells = 0: open
#if defined(USE_CUDA)
                ryoanji::directSum<double>(
                        first, first + nLocal, nAll, unusedBox, 0, device[0].data(),
                        device[1].data(), device[2].data(), device[3].data(), device[4].data(),
                        potential.data(), ax.data(), ay.data(), az.data());
#else
                // The CPU variant evaluates every target; only this rank's slice is used below.
                ryoanji::directSum(
                        device[0].data(), device[1].data(), device[2].data(), device[4].data(),
                        device[3].data(), static_cast<cstone::LocalIndex>(nAll), 1.0f, unusedBox,
                        0, ax.data(), ay.data(), az.data(), potential.data());
#endif
                synchronizeDevice();
                unpackDirectSum(destination, ax, ay, az, first, nLocal, prefactor);
            };
            evaluate(particles.E.getView());
            if (!imageField) {
                return;
            }

            // Image charges alone: each rank's block is its real particles, then as many images.
            std::vector<double> imageCharge = columns[3];
            for (int r = 0; r < ranks; ++r) {
                const std::size_t begin = static_cast<std::size_t>(offsets[r]) / kValues;
                const std::size_t real  = static_cast<std::size_t>(counts[r]) / kValues / 2;
                std::fill_n(imageCharge.begin() + begin, real, 0.0);
            }
            Kokkos::deep_copy(device[3], HostUnmanaged(imageCharge.data(), nAll));
            for (auto* output : {&potential, &ax, &ay, &az}) {
                Kokkos::deep_copy(*output, 0.0);
            }
            Kokkos::fence();
            evaluate(particles.B.getView());
        }

        /**
         * @brief Tree-frame description of per-particle source velocities: z stretched by the mean
         * gamma of @p boost, sources described by their rapidity relative to it.
         */
        ryoanji::HeavisideFrame<double> heavisideFrameOf(const RestFrameBoost& boost) {
            ryoanji::HeavisideFrame<double> frame;
            frame.gammaRef = boost.gamma;
            frame.wRef     = referenceRapidity(boost);
            frame.betaRef  = std::tanh(frame.wRef);
            // Weighting the second channel by 1 - beta_ref beta keeps E_perp - c beta B accurate for
            // nearly co-moving particles at high gamma (no cancellation of two large sums).
            frame.netForceChannel = std::abs(frame.betaRef) > 0.5;
            return frame;
        }

        /**
         * @brief Field of every source moving uniformly along z with its own velocity, by direct
         * summation (BHRELATIVISTIC=PARTICLE). Particles do not migrate.
         *
         * Positions are in the tree frame (z stretched by @p frame.gammaRef, @p stretch); each
         * source carries its rapidity relative to the frame. Image charges of the plane @p planeZ
         * (tree frame) carry -q and move against their particle, -w. On return E holds the
         * tree-frame field E~ and B the second channel C (x, y; scratch for
         * boostParticleKinematics), both with @p prefactor applied as in directSumFields.
         *
         * @return the largest gamma of the sources (for the field bound of verifySolve)
         */
        double directSumHeavisideFields(
                BarnesHutAlgorithm::ParticleContainer& particles, double chargeScale,
                double softening, double prefactor, bool images, double planeZ, ZStretch stretch,
                const ryoanji::HeavisideFrame<double>& frame) {
            const size_type nLocal   = particles.getLocalNum();
            const size_type nSources = images ? 2 * nLocal : nLocal;
            const auto R             = particles.R.getView();
            const auto P             = particles.P.getView();
            const auto dtView        = particles.dt.getView();
            const double wRef        = frame.wRef;

            // (x, y, z, q, h, w - wRef) per local particle, then the image charges.
            constexpr int kValues = 6;
            Kokkos::View<double*> local("BarnesHutAlgorithm::heavisideLocal", kValues * nSources);
            double maxRapidity = 0.0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::heavisidePack", nLocal,
                    KOKKOS_LAMBDA(const size_type i, double& maxW) {
                        const double z         = stretch.apply(R(i)[2]);
                        const double w         = longitudinalRapidity(P(i));
                        local(kValues * i)     = R(i)[0];
                        local(kValues * i + 1) = R(i)[1];
                        local(kValues * i + 2) = z;
                        local(kValues * i + 3) = chargeScale * dtView(i);
                        local(kValues * i + 4) = softening;
                        local(kValues * i + 5) = w - wRef;
                        if (images) {
                            const size_type k      = nLocal + i;
                            local(kValues * k)     = R(i)[0];
                            local(kValues * k + 1) = R(i)[1];
                            local(kValues * k + 2) = 2.0 * planeZ - z;
                            local(kValues * k + 3) = -chargeScale * dtView(i);
                            local(kValues * k + 4) = softening;
                            local(kValues * k + 5) = -w - wRef;
                        }
                        maxW = Kokkos::fabs(w) > maxW ? Kokkos::fabs(w) : maxW;
                    },
                    Kokkos::Max<double>(maxRapidity));
            const auto localHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), local);
            double globalMaxRapidity = 0.0;
            ippl::Comm->allreduce(maxRapidity, globalMaxRapidity, 1, std::greater<double>());

            MPI_Comm comm   = ippl::Comm->getCommunicator();
            const int ranks = ippl::Comm->size();
            const int count = static_cast<int>(kValues * nSources);
            std::vector<int> counts(ranks), offsets(ranks, 0);
            MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            for (int r = 1; r < ranks; ++r) {
                offsets[r] = offsets[r - 1] + counts[r - 1];
            }
            std::vector<double> gathered(
                    static_cast<std::size_t>(offsets.back() + counts.back()));
            MPI_Allgatherv(
                    localHost.data(), count, MPI_DOUBLE, gathered.data(), counts.data(),
                    offsets.data(), MPI_DOUBLE, comm);

            const std::size_t nAll  = gathered.size() / kValues;
            const std::size_t first = static_cast<std::size_t>(offsets[ippl::Comm->rank()]) / kValues;
            std::array<std::vector<double>, kValues> columns;
            for (auto& column : columns) {
                column.resize(nAll);
            }
            for (std::size_t j = 0; j < nAll; ++j) {
                for (int c = 0; c < kValues; ++c) {
                    columns[c][j] = gathered[kValues * j + c];
                }
            }
            using HostUnmanaged = Kokkos::View<const double*, Kokkos::HostSpace,
                                               Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
            std::array<Kokkos::View<double*>, kValues> device;
            for (int c = 0; c < kValues; ++c) {
                device[c] = Kokkos::View<double*>("BarnesHutAlgorithm::heavisideSource", nAll);
                Kokkos::deep_copy(device[c], HostUnmanaged(columns[c].data(), nAll));
            }
            // Outputs accumulate (+=); the views are zero-initialised.
            std::array<Kokkos::View<double*>, 5> out;
            for (auto& view : out) {
                view = Kokkos::View<double*>("BarnesHutAlgorithm::heavisideField", nAll);
            }
            Kokkos::fence();

            const double maxGamma = std::cosh(globalMaxRapidity);
            if (nLocal == 0) {
                return maxGamma;
            }
#if defined(USE_CUDA)
            ryoanji::directSumHeaviside<double>(
                    first, first + nLocal, nAll, frame, device[0].data(), device[1].data(),
                    device[2].data(), device[3].data(), device[4].data(), device[5].data(),
                    out[0].data(), out[1].data(), out[2].data(), out[3].data(), out[4].data());
#else
            ryoanji::directSumHeavisideCpu<double>(
                    static_cast<cstone::LocalIndex>(first),
                    static_cast<cstone::LocalIndex>(first + nLocal),
                    static_cast<cstone::LocalIndex>(nAll), frame, device[0].data(),
                    device[1].data(), device[2].data(), device[3].data(), device[4].data(),
                    device[5].data(), out[0].data(), out[1].data(), out[2].data(), out[3].data(),
                    out[4].data());
#endif
            synchronizeDevice();

            auto E = particles.E.getView();
            auto B = particles.B.getView();
            const auto ex = out[0], ey = out[1], ez = out[2], cx = out[3], cy = out[4];
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::heavisideUnpack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const std::size_t j = first + i;
                        E(i) = Vector_t<double, 3>(
                                prefactor * ex(j), prefactor * ey(j), prefactor * ez(j));
                        B(i) = Vector_t<double, 3>(prefactor * cx(j), prefactor * cy(j), 0.0);
                    });
            Kokkos::fence();
            return maxGamma;
        }

        /**
         * @brief Tree-frame fields of per-particle source velocities to the lab (solve axes): E holds
         * E~, B the second channel C = sum chi_j E~_perp,j.
         *
         * E_perp = gammaRef E~_perp, E_z = E~_z, and B = (gammaRef / c) z x sum_j beta_j E~_perp,j,
         * where the beta-weighted sum is C, or (E~_perp - C) / betaRef for the net-force channel
         * chi = 1 - betaRef beta. With all sources at the frame velocity this is boostToLab.
         */
        void boostParticleKinematics(
                BarnesHutAlgorithm::ParticleContainer& particles,
                const ryoanji::HeavisideFrame<double>& frame) {
            auto E                = particles.E.getView();
            auto B                = particles.B.getView();
            const double gamma    = frame.gammaRef;
            const double betaRef  = frame.betaRef;
            const bool net        = frame.netForceChannel;
            const double bScale   = gamma / Physics::c;
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::boostParticleKinematics", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const size_type i) {
                        const Vector_t<double, 3> field   = E(i);
                        const Vector_t<double, 3> channel = B(i);
                        const double vx = net ? (field[0] - channel[0]) / betaRef : channel[0];
                        const double vy = net ? (field[1] - channel[1]) / betaRef : channel[1];
                        E(i) = Vector_t<double, 3>(gamma * field[0], gamma * field[1], field[2]);
                        B(i) = Vector_t<double, 3>(-bScale * vy, bScale * vx, 0.0);
                    });
            Kokkos::fence();
        }

        /**
         * @brief Mean particle spacing in the solve frame with rest-frame z (BHSOFTENINGREL).
         *
         * The bunch counts as a box with edges sqrt(2 pi) sigma per axis (the volume of a
         * Gaussian with these rms sizes) holding @p total particles. An axis shorter than the
         * spacing holds a single layer, so a sheet (emission from a cathode) gets its 2D spacing
         * and a line its 1D spacing. Moments are taken about the centre of @p bounds.
         */
        double meanSpacing(
                BarnesHutAlgorithm::ParticleContainer& particles, ZStretch stretch,
                const std::array<double, 6>& bounds, unsigned long long total) {
            const auto R           = particles.R.getView();
            const size_type nLocal = particles.getLocalNum();
            const double ox = 0.5 * (bounds[0] + bounds[1]), oy = 0.5 * (bounds[2] + bounds[3]),
                         oz = 0.5 * (bounds[4] + bounds[5]);
            ippl::Vector<double, 6> local(0.0);
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::spacingMoments", nLocal,
                    KOKKOS_LAMBDA(const size_type i, ippl::Vector<double, 6>& sum) {
                        const double d[3] = {
                                R(i)[0] - ox, R(i)[1] - oy, stretch.apply(R(i)[2]) - oz};
                        for (unsigned a = 0; a < 3; ++a) {
                            sum[a] += d[a];
                            sum[3 + a] += d[a] * d[a];
                        }
                    },
                    local);
            std::array<double, 6> m{};
            for (unsigned c = 0; c < 6; ++c) {
                m[c] = local[c];
            }
            ippl::Comm->allreduce(m.data(), 6, std::plus<double>());
            const double n = static_cast<double>(total);
            std::array<double, 3> edge{};
            for (unsigned a = 0; a < 3; ++a) {
                const double mean = m[a] / n;
                edge[a] = std::sqrt(2.0 * Physics::pi * std::max(0.0, m[3 + a] / n - mean * mean));
            }
            std::sort(edge.begin(), edge.end());
            double spacing = 0.0;
            for (unsigned thin = 0; thin < 3; ++thin) {
                double volume = 1.0;
                for (unsigned a = thin; a < 3; ++a) {
                    volume *= edge[a];
                }
                spacing = std::pow(volume / n, 1.0 / static_cast<double>(3 - thin));
                if (edge[thin] > 0.0 && edge[thin] >= spacing) {
                    break;
                }
            }
            return spacing;
        }

        /**
         * @brief True if the image charges of the plane @p planeZ are far enough from the sources
         * in @p bounds (both rest frame) to be one multipole: their extent seen from the nearest
         * source passes the multipole acceptance criterion with @p theta.
         */
        bool imagesAreFar(const std::array<double, 6>& bounds, double planeZ, double theta) {
            const double lx = bounds[1] - bounds[0], ly = bounds[3] - bounds[2],
                         lz = bounds[5] - bounds[4];
            const double radius   = 0.5 * std::sqrt(lx * lx + ly * ly + lz * lz);
            const double distance = 2.0 * std::abs(0.5 * (bounds[4] + bounds[5]) - planeZ);
            return distance > radius && 2.0 * radius < theta * (distance - radius);
        }

        /**
         * @brief Add the field of far image charges, as one quadrupole about their centre of
         * charge, to E (all sources) and B (image field, scratch for the boost).
         *
         * The image of particle i has charge -sign |q| dt_i / dt_step at the rest-frame position
         * mirrored at @p planeZ. Second moments are summed about the centre of @p bounds (rest
         * frame) and shifted to the centre of charge, where the dipole vanishes. With
         * S = sum q d d^T about the centre, the field at r from the centre is
         * k [Q r / r^3 - 3 S r / r^5 + 15/2 (r.S r) r / r^7 - 3/2 tr(S) r / r^5].
         */
        void addFarImageField(
                BarnesHutAlgorithm::ParticleContainer& particles, ZStretch stretch,
                const std::array<double, 6>& bounds, double planeZ, double chargeScale,
                double sign) {
            const auto R           = particles.R.getView();
            const auto dt          = particles.dt.getView();
            const size_type nLocal = particles.getLocalNum();
            const double ox = 0.5 * (bounds[0] + bounds[1]), oy = 0.5 * (bounds[2] + bounds[3]),
                         oz = 0.5 * (bounds[4] + bounds[5]);

            // w, w d, w d d^T (xx, yy, zz, xy, xz, yz) of the real particles, w = dt_i.
            ippl::Vector<double, 10> local(0.0);
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::imageMoments", nLocal,
                    KOKKOS_LAMBDA(const size_type i, ippl::Vector<double, 10>& sum) {
                        const double w  = dt(i);
                        const double dx = R(i)[0] - ox, dy = R(i)[1] - oy,
                                     dz = stretch.apply(R(i)[2]) - oz;
                        sum[0] += w;
                        sum[1] += w * dx;
                        sum[2] += w * dy;
                        sum[3] += w * dz;
                        sum[4] += w * dx * dx;
                        sum[5] += w * dy * dy;
                        sum[6] += w * dz * dz;
                        sum[7] += w * dx * dy;
                        sum[8] += w * dx * dz;
                        sum[9] += w * dy * dz;
                    },
                    local);
            std::array<double, 10> m{};
            for (unsigned c = 0; c < 10; ++c) {
                m[c] = local[c];
            }
            ippl::Comm->allreduce(m.data(), 10, std::plus<double>());
            if (!(m[0] > 0.0) || nLocal == 0) {
                return;
            }

            const double cx = m[1] / m[0], cy = m[2] / m[0], cz = m[3] / m[0];
            // Image charge per unit weight times the Coulomb constant; mirroring z negates the
            // xz and yz moments.
            const double kq = -sign * chargeScale / (4.0 * Physics::pi * Physics::epsilon_0);
            const double charge = kq * m[0];
            const double sxx = kq * (m[4] - m[0] * cx * cx), syy = kq * (m[5] - m[0] * cy * cy),
                         szz = kq * (m[6] - m[0] * cz * cz);
            const double sxy = kq * (m[7] - m[0] * cx * cy), sxz = -kq * (m[8] - m[0] * cx * cz),
                         syz = -kq * (m[9] - m[0] * cy * cz);
            const double trace = sxx + syy + szz;
            const double centreX = ox + cx, centreY = oy + cy, centreZ = 2.0 * planeZ - (oz + cz);

            auto E = particles.E.getView();
            auto B = particles.B.getView();
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::farImageField", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const double rx = R(i)[0] - centreX, ry = R(i)[1] - centreY,
                                     rz = stretch.apply(R(i)[2]) - centreZ;
                        const double r2   = rx * rx + ry * ry + rz * rz;
                        const double inv  = 1.0 / Kokkos::sqrt(r2);
                        const double inv2 = inv * inv;
                        const double inv3 = inv * inv2;
                        const double inv5 = inv3 * inv2;
                        const double srx  = sxx * rx + sxy * ry + sxz * rz;
                        const double sry  = sxy * rx + syy * ry + syz * rz;
                        const double srz  = sxz * rx + syz * ry + szz * rz;
                        const double rsr  = rx * srx + ry * sry + rz * srz;
                        const double radial =
                                charge * inv3 + 7.5 * rsr * inv5 * inv2 - 1.5 * trace * inv5;
                        const Vector_t<double, 3> field(
                                radial * rx - 3.0 * inv5 * srx, radial * ry - 3.0 * inv5 * sry,
                                radial * rz - 3.0 * inv5 * srz);
                        E(i) += field;
                        B(i) += field;
                    });
            Kokkos::fence();
        }

        /// Image charges of the plane as one multipole with per-particle velocities.
        struct FarImageMultipole {
            ryoanji::HeavisideQpole<double> multipole;  ///< Evaluation form.
            ryoanji::Vec3<double> centre{0.0, 0.0, 0.0};  ///< Tree frame.
            double macScale2 = 1.0;  ///< heavisideMacScale2; negative: too wide in rapidity.
            bool empty       = true;
        };

        /**
         * @brief Multipole (evaluation form, about the centre of charge) of the image charges of
         * the plane @p planeZ, every image moving against its particle (BHRELATIVISTIC=PARTICLE).
         *
         * The image of particle i sits at the tree-frame position mirrored at @p planeZ, carries
         * -|q_i| in the packed charge convention and the rapidity -w_i.
         */
        FarImageMultipole farImageMultipole(
                BarnesHutAlgorithm::ParticleContainer& particles, ZStretch stretch,
                const std::array<double, 6>& bounds, double planeZ, double chargeScale,
                const ryoanji::HeavisideFrame<double>& frame, double thetaV, bool velocityMoments) {
            using ryoanji::Hqi;
            const auto R           = particles.R.getView();
            const auto P           = particles.P.getView();
            const auto dt          = particles.dt.getView();
            const size_type nLocal = particles.getLocalNum();
            // Raw moments about the mirrored centre of the bounds.
            const double ox = 0.5 * (bounds[0] + bounds[1]), oy = 0.5 * (bounds[2] + bounds[3]),
                         oz = 2.0 * planeZ - 0.5 * (bounds[4] + bounds[5]);
            const double wRef = frame.wRef;

            // mass, absMass, P (3), S (6), absW, v1, X (3), v2: the summed entries of the raw form
            constexpr unsigned kSums = 17;
            ippl::Vector<double, kSums> local(0.0);
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::farImageMoments", nLocal,
                    KOKKOS_LAMBDA(const size_type i, ippl::Vector<double, kSums>& sum) {
                        const double q  = -chargeScale * dt(i);
                        const double dx = R(i)[0] - ox, dy = R(i)[1] - oy,
                                     dz = 2.0 * planeZ - stretch.apply(R(i)[2]) - oz;
                        const double w  = -longitudinalRapidity(P(i)) - wRef;
                        const double terms[kSums] = {q,          -q,          q * dx,     q * dy,
                                                     q * dz,     q * dx * dx, q * dx * dy, q * dx * dz,
                                                     q * dy * dy, q * dy * dz, q * dz * dz, -q * w,
                                                     q * w,      q * dx * w,  q * dy * w, q * dz * w,
                                                     q * w * w};
                        for (unsigned c = 0; c < kSums; ++c) {
                            sum[c] += terms[c];
                        }
                    },
                    local);
            using MinMax = Kokkos::MinMax<double>;
            typename MinMax::value_type range;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::farImageRapidities", nLocal,
                    KOKKOS_LAMBDA(const size_type i, typename MinMax::value_type& update) {
                        const double w = -longitudinalRapidity(P(i)) - wRef;
                        update.min_val = w < update.min_val ? w : update.min_val;
                        update.max_val = w > update.max_val ? w : update.max_val;
                    },
                    MinMax(range));
            std::array<double, kSums> sums{};
            for (unsigned c = 0; c < kSums; ++c) {
                sums[c] = local[c];
            }
            ippl::Comm->allreduce(sums.data(), kSums, std::plus<double>());
            std::array<double, 2> extremes{-range.min_val, range.max_val}, global{};
            ippl::Comm->allreduce(extremes.data(), global.data(), 2, std::greater<double>());

            FarImageMultipole result;
            if (!(sums[1] > 0.0)) {
                return result;
            }
            ryoanji::HeavisideQpole<double> raw = ryoanji::heavisideZero<double>();
            const int order[kSums] = {Hqi::mass, Hqi::absMass, Hqi::px,  Hqi::py,  Hqi::pz,  Hqi::sxx,
                                      Hqi::sxy,  Hqi::sxz,     Hqi::syy, Hqi::syz, Hqi::szz, Hqi::absW,
                                      Hqi::v1,   Hqi::xx,      Hqi::xy,  Hqi::xz,  Hqi::v2};
            for (unsigned c = 0; c < kSums; ++c) {
                raw[order[c]] = sums[c];
            }
            raw[Hqi::wMax]    = global[1];
            raw[Hqi::wNegMin] = global[0];
            // All images share one sign: the centre of charge is the |q|-weighted centre.
            const double cx = raw[Hqi::px] / raw[Hqi::mass], cy = raw[Hqi::py] / raw[Hqi::mass],
                         cz = raw[Hqi::pz] / raw[Hqi::mass];
            result.multipole = ryoanji::heavisideZero<double>();
            ryoanji::heavisideM2M(result.multipole, ryoanji::Vec3<double>{cx, cy, cz}, raw);
            ryoanji::heavisideFinalize(result.multipole, frame, velocityMoments);
            result.macScale2 =
                    ryoanji::heavisideMacScale2(result.multipole, frame, thetaV, velocityMoments);
            result.centre = {ox + cx, oy + cy, oz + cz};
            result.empty  = false;
            return result;
        }

        /**
         * @brief Add the field of far image charges (FarImageMultipole) to E (tree-frame E~) and B
         * (second channel, scratch for boostParticleKinematics), with @p prefactor as the direct
         * sum.
         */
        void addFarImageFieldHeaviside(
                BarnesHutAlgorithm::ParticleContainer& particles, ZStretch stretch,
                const FarImageMultipole& images, double prefactor) {
            if (images.empty) {
                return;
            }
            const auto R                              = particles.R.getView();
            auto E                                    = particles.E.getView();
            auto B                                    = particles.B.getView();
            const ryoanji::HeavisideQpole<double> M   = images.multipole;
            const ryoanji::Vec3<double> centre        = images.centre;
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::farImageFieldHeaviside", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const size_type i) {
                        const ryoanji::Vec3<double> target{
                                R(i)[0], R(i)[1], stretch.apply(R(i)[2])};
                        ::util::array<double, 5> acc{0.0, 0.0, 0.0, 0.0, 0.0};
                        acc = ryoanji::M2PHeaviside(acc, target, centre, M);
                        E(i) += Vector_t<double, 3>(
                                prefactor * acc[0], prefactor * acc[1], prefactor * acc[2]);
                        B(i) += Vector_t<double, 3>(prefactor * acc[3], prefactor * acc[4], 0.0);
                    });
            Kokkos::fence();
        }

        /** @brief Largest longitudinal gamma of the local particles over all ranks. */
        double maxLongitudinalGamma(BarnesHutAlgorithm::ParticleContainer& particles) {
            const auto P = particles.P.getView();
            double maxW  = 0.0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::maxRapidity", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const size_type i, double& m) {
                        const double w = Kokkos::fabs(longitudinalRapidity(P(i)));
                        m              = w > m ? w : m;
                    },
                    Kokkos::Max<double>(maxW));
            double global = 0.0;
            ippl::Comm->allreduce(maxW, global, 1, std::greater<double>());
            return std::cosh(global);
        }

        /**
         * @brief Rest-frame fields to the lab frame: E holds the field of all sources, B the
         * field of the image charges alone (zero without images). Image charges move against
         * the bunch, so their B has the opposite sign: B = (gamma/c^2) v x (E'_all - 2 E'_image).
         */
        void boostToLab(
                BarnesHutAlgorithm::ParticleContainer& particles, const RestFrameBoost& boost) {
            auto E = particles.E.getView();
            auto B = particles.B.getView();
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::boostToLab", particles.getLocalNum(),
                    KOKKOS_LAMBDA(const size_type i) {
                        const Vector_t<double, 3> all   = E(i);
                        const Vector_t<double, 3> image = B(i);
                        const Vector_t<double, 3> moving = all - 2.0 * image;
                        E(i) = boost.electric(all);
                        B(i) = boost.magnetic(moving, 1.0);
                    });
            Kokkos::fence();
        }

        /** @brief Local field and charge statistics of the primary after a solve. */
        struct FieldStatistics {
            double nonFinite = 0.0;  ///< Particles with a non-finite E or B component.
            double maxField  = 0.0;  ///< Largest |E| in V/m.
            double sumDt     = 0.0;  ///< Sum of dt_i; total source charge is |Q| sum(dt_i) / dt_step.
        };

        FieldStatistics globalFieldStatistics(BarnesHutAlgorithm::ParticleContainer& particles) {
            const auto E           = particles.E.getView();
            const auto B           = particles.B.getView();
            const auto dt          = particles.dt.getView();
            const size_type nLocal = particles.getLocalNum();
            double nonFinite = 0.0, maxE2 = 0.0, sumDt = 0.0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::fieldStatistics", nLocal,
                    KOKKOS_LAMBDA(const size_type i, double& bad, double& m, double& s) {
                        const double e2 = E(i)[0] * E(i)[0] + E(i)[1] * E(i)[1] + E(i)[2] * E(i)[2];
                        const double b2 = B(i)[0] * B(i)[0] + B(i)[1] * B(i)[1] + B(i)[2] * B(i)[2];
                        if (!Kokkos::isfinite(e2) || !Kokkos::isfinite(b2)) {
                            bad += 1.0;
                        } else if (e2 > m) {
                            m = e2;
                        }
                        s += dt(i);
                    },
                    nonFinite, Kokkos::Max<double>(maxE2), sumDt);
            FieldStatistics global;
            ippl::Comm->allreduce(nonFinite, global.nonFinite, 1, std::plus<double>());
            ippl::Comm->allreduce(maxE2, global.maxField, 1, std::greater<double>());
            ippl::Comm->allreduce(sumDt, global.sumDt, 1, std::plus<double>());
            global.maxField = std::sqrt(global.maxField);
            return global;
        }

        /// BHSTATS: the live Kokkos allocations in device memory of this rank, from the Kokkos
        /// Tools allocation callbacks. They are registered by the first solver with BHSTATS unless
        /// a Kokkos tool is loaded; allocations made before that are not seen.
        struct KokkosDeviceAllocations {
            std::unordered_map<const void*, std::pair<std::string, std::uint64_t>> live;
            std::uint64_t bytes = 0;
            bool active         = false;
        };

        KokkosDeviceAllocations& kokkosDeviceAllocations() {
            static KokkosDeviceAllocations allocations;
            return allocations;
        }

        bool isDeviceSpace(const Kokkos_Profiling_SpaceHandle& space) {
            return std::strcmp(space.name, "Host") != 0
                   && std::strstr(space.name, "HostPinned") == nullptr;
        }

        void onKokkosAllocate(const Kokkos_Profiling_SpaceHandle space, const char* label,
                              const void* ptr, const std::uint64_t size) {
            if (!isDeviceSpace(space)) {
                return;
            }
            KokkosDeviceAllocations& a = kokkosDeviceAllocations();
            a.live[ptr]                = {label, size};
            a.bytes += size;
        }

        void onKokkosDeallocate(const Kokkos_Profiling_SpaceHandle space, const char*,
                                const void* ptr, const std::uint64_t) {
            if (!isDeviceSpace(space)) {
                return;
            }
            KokkosDeviceAllocations& a = kokkosDeviceAllocations();
            const auto it              = a.live.find(ptr);
            if (it == a.live.end()) {
                return;  // allocated before the callbacks were registered
            }
            a.bytes -= it->second.second;
            a.live.erase(it);
        }

        void trackKokkosDeviceAllocations() {
            KokkosDeviceAllocations& a = kokkosDeviceAllocations();
            if (a.active || Kokkos::Tools::profileLibraryLoaded()) {
                return;
            }
            Kokkos::Tools::Experimental::set_allocate_data_callback(onKokkosAllocate);
            Kokkos::Tools::Experimental::set_deallocate_data_callback(onKokkosDeallocate);
            a.active = true;
        }

        /// Device memory in use on this rank's GPU (all allocations, CUDA context included), in GB.
        double deviceUsedGB() {
#if defined(KOKKOS_ENABLE_CUDA)
            std::size_t free = 0, total = 0;
            if (cudaMemGetInfo(&free, &total) == cudaSuccess) {
                return static_cast<double>(total - free) / 1073741824.0;
            }
#endif
            return 0.0;
        }

        /// Memory of the device's default stream-ordered pool (cudaMallocAsync, used by cstone for
        /// sort and tree-build temporaries), in GB: {reserved, used, release threshold}.
        /// OPALX_BH_TRIM_POOL=1 trims the pool to zero first (diagnostic).
        std::array<double, 3> asyncPoolGB() {
            std::array<double, 3> gb{};
#if defined(KOKKOS_ENABLE_CUDA)
            int device = 0;
            cudaMemPool_t pool{};
            if (cudaGetDevice(&device) != cudaSuccess
                || cudaDeviceGetDefaultMemPool(&pool, device) != cudaSuccess) {
                return gb;
            }
            static const bool trim = [] {
                const char* value = std::getenv("OPALX_BH_TRIM_POOL");
                return value != nullptr && std::strcmp(value, "1") == 0;
            }();
            if (trim) {
                cudaDeviceSynchronize();
                cudaMemPoolTrimTo(pool, 0);
            }
            const cudaMemPoolAttr attributes[3] = {cudaMemPoolAttrReservedMemCurrent,
                                                   cudaMemPoolAttrUsedMemCurrent,
                                                   cudaMemPoolAttrReleaseThreshold};
            for (int i = 0; i < 3; ++i) {
                std::uint64_t bytes = 0;
                if (cudaMemPoolGetAttribute(pool, attributes[i], &bytes) == cudaSuccess) {
                    gb[i] = static_cast<double>(bytes) / 1073741824.0;
                }
            }
#endif
            return gb;
        }

        /**
         * BHSTATS: where the device memory of every rank is, after one solve. One [bh.mem] line per
         * rank (up to 16 ranks; beyond, the ranks with the least and the most memory), and the
         * largest live Kokkos allocations of the rank with the most memory as [bh.mem.top].
         * opal: capacity of the OPAL-X particle attributes; bh: reserved buffers of the BH
         * container; kokkos: all live Kokkos device allocations (opal included when they were
         * allocated after the first BHSTATS solver); other: the rest (cstone domain, multipoles,
         * MPI, CUDA context), of which pool is held by the cudaMallocAsync pool. Collective.
         */
        void printMemoryReport(unsigned long long step, double devStart, std::size_t opalLocal,
                               std::size_t opalCapacity, std::size_t opalBytes,
                               std::size_t bhSlots, std::size_t bhBytes) {
            constexpr double GB              = 1073741824.0;
            const KokkosDeviceAllocations& a = kokkosDeviceAllocations();
            const std::array<double, 3> pool = asyncPoolGB();
            const double dev                 = deviceUsedGB();
            const double kokkos              = a.active ? static_cast<double>(a.bytes) / GB : 0.0;
            const double bh                  = static_cast<double>(bhBytes) / GB;
            constexpr int nValues            = 11;
            const double local[nValues]      = {dev,
                                                devStart,
                                                kokkos,
                                                static_cast<double>(opalBytes) / GB,
                                                static_cast<double>(opalLocal),
                                                static_cast<double>(opalCapacity),
                                                bh,
                                                static_cast<double>(bhSlots),
                                                pool[0],
                                                pool[1],
                                                pool[2]};
            const MPI_Comm comm = ippl::Comm->getCommunicator();
            const int rank      = ippl::Comm->rank();
            const int ranks     = ippl::Comm->size();
            std::vector<double> all(rank == 0 ? static_cast<std::size_t>(nValues * ranks) : 0);
            MPI_Gather(local, nValues, MPI_DOUBLE, all.data(), nValues, MPI_DOUBLE, 0, comm);

            struct {
                double value;
                int rank;
            } mine{dev, rank}, top{};
            MPI_Allreduce(&mine, &top, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

            if (rank == 0) {
                std::vector<int> shown;
                if (ranks <= 16) {
                    for (int r = 0; r < ranks; ++r) {
                        shown.push_back(r);
                    }
                } else {
                    int least = 0;
                    for (int r = 1; r < ranks; ++r) {
                        least = all[nValues * r] < all[nValues * least] ? r : least;
                    }
                    shown = {least, top.rank};
                }
                for (const int r : shown) {
                    const double* v = &all[static_cast<std::size_t>(nValues * r)];
                    std::printf("[bh.mem] step=%llu rank=%d dev=%.2f devStart=%.2f kokkos=%.2f "
                                "opal=%.2f (local %.3fM cap %.3fM) bh=%.2f (slots %.3fM) "
                                "other=%.2f pool=%.2f (used %.2f threshold %.3g) GB\n",
                                step, r, v[0], v[1], v[2], v[3], v[4] / 1e6, v[5] / 1e6, v[6],
                                v[7] / 1e6, v[0] - v[2] - v[6], v[8], v[9], v[10]);
                }
                std::fflush(stdout);
            }
            if (rank == top.rank && a.active) {
                std::unordered_map<std::string, std::uint64_t> byLabel;
                for (const auto& [ptr, allocation] : a.live) {
                    byLabel[allocation.first] += allocation.second;
                }
                std::vector<std::pair<std::string, std::uint64_t>> sorted(
                        byLabel.begin(), byLabel.end());
                std::sort(sorted.begin(), sorted.end(),
                          [](const auto& x, const auto& y) { return x.second > y.second; });
                std::ostringstream line;
                line << "[bh.mem.top] step=" << step << " rank=" << rank << " allocations="
                     << a.live.size();
                for (std::size_t i = 0; i < std::min<std::size_t>(sorted.size(), 6); ++i) {
                    line << " " << sorted[i].first << "="
                         << static_cast<double>(sorted[i].second) / GB;
                }
                std::printf("%s GB\n", line.str().c_str());
                std::fflush(stdout);
            }
        }

        /// Stage timers of one solve, in the order SolveTimeReport prints them.
        constexpr std::array<const char*, 15> kSolveTimers = {
                "bh.restFrame", "bh.domain",   "bh.syncGrav",  "bh.syncGrav.new", "bh.haloCharge",
                "bh.pack",      "bh.pack.new", "bh.groups",    "bh.upsweep",      "bh.compute",
                "bh.imagePass", "bh.copyBack", "bh.copyBack.new", "bh.farImages", "bh.directSum"};

        /**
         * BHSTATS: wall time of every stage of one solve, maximum over ranks, as one line on rank 0.
         * IpplTimings keeps every measurement, so the stages are the measurements added since
         * construction (the IPPL traversal stages included). Collective: all ranks must print().
         */
        class SolveTimeReport {
        public:
            SolveTimeReport() : start_m(MPI_Wtime()) {
                for (std::size_t i = 0; i < kSolveTimers.size(); ++i) {
                    timers_m[i] = IpplTimings::getTimer(kSolveTimers[i]);
                    counts_m[i] = IpplTimings::getMeasurementCount(timers_m[i]);
                }
            }

            void print(unsigned long long step, unsigned long long particles, double gamma,
                       double softening, bool images) const {
                std::array<double, kSolveTimers.size() + 1> local{}, global{};
                for (std::size_t i = 0; i < kSolveTimers.size(); ++i) {
                    const std::vector<double>& m = IpplTimings::getMeasurements(timers_m[i]);
                    for (std::size_t j = counts_m[i]; j < m.size(); ++j) {
                        local[i] += m[j];
                    }
                }
                local.back() = MPI_Wtime() - start_m;
                MPI_Allreduce(local.data(), global.data(), static_cast<int>(local.size()),
                              MPI_DOUBLE, MPI_MAX, ippl::Comm->getCommunicator());
                if (ippl::Comm->rank() != 0) {
                    return;
                }
                // pack.new / copyBack.new / syncGrav.new: this solve built the domain.
                const bool newDomain =
                        IpplTimings::getMeasurementCount(timers_m[6]) > counts_m[6];
                std::printf("[bh.time] step=%llu N=%llu gamma=%.5g h=%.4g images=%d new=%d "
                            "total=%.4f restFrame=%.4f domain=%.4f sync=%.4f halo=%.4f "
                            "pack=%.4f groups=%.4f upsweep=%.4f traverse=%.4f imagePass=%.4f "
                            "copyBack=%.4f farImages=%.4f direct=%.4f\n",
                            step, particles, gamma, softening, images ? 1 : 0, newDomain ? 1 : 0,
                            global.back(), global[0], global[1], global[2] + global[3], global[4],
                            global[5] + global[6], global[7], global[8], global[9], global[10],
                            global[11] + global[12], global[13], global[14]);
                std::fflush(stdout);
            }

        private:
            double start_m;
            std::array<IpplTimings::TimerRef, kSolveTimers.size()> timers_m{};
            std::array<std::size_t, kSolveTimers.size()> counts_m{};
        };

    }  // namespace

    /// Precision-independent interface of the implementation (BHPRECISION picks the policy).
    class BarnesHutAlgorithm::Impl {
    public:
        virtual ~Impl() = default;
        virtual SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) = 0;
        [[nodiscard]] virtual const Statistics& statistics() const    = 0;
        [[nodiscard]] virtual const BarnesHutConfig& config() const   = 0;
    };

    /// Implementation for one IPPL NBody precision policy (DoublePrecision, MixedPrecision).
    template <class Precision>
    class BarnesHutAlgorithm::PrecisionImpl final : public BarnesHutAlgorithm::Impl {
    public:
        using BHContainer  = ippl::nbody::NBodyParticleContainer<Precision, 3>;
        using BHSolver     = ippl::nbody::NBodySolver<Precision, 3>;
        using BHLocalIndex = typename BHContainer::LocalIndex;

        PrecisionImpl(
                BarnesHutConfig config, ParticleContainer& primary,
                std::shared_ptr<const BunchStateHandler> bunchState)
            : config_m(std::move(config)), primary_m(&primary), bunchState_m(std::move(bunchState)) {
            if (config_m.traversalStatistics) {
                trackKokkosDeviceAllocations();
            }
        }

        SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override;
        [[nodiscard]] const Statistics& statistics() const override { return statistics_m; }
        [[nodiscard]] const BarnesHutConfig& config() const override { return config_m; }

        BarnesHutConfig config_m;
        ParticleContainer* primary_m = nullptr;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<BHContainer> bh_m;
        std::unique_ptr<BHSolver> solver_m;
        ippl::nbody::MultipoleOrder solverOrder_m = ippl::nbody::MultipoleOrder::Quadrupole;
        unsigned bucketSize_m      = 0;
        unsigned bucketSizeFocus_m = 0;
        bool usingDirectSum_m      = false;
        std::string memoryTrace_m;  // BHSTATS: device memory at the stages of this solve
        bool newDomainSolve_m = false;  // this solve built a new cstone domain

        /// BHSTATS: record the device memory in use (after a device sync) at @p stage of this
        /// solve; printed as [bh.mem.trace] with the memory report.
        void traceMemory(const char* stage) {
            if (!config_m.traversalStatistics) {
                return;
            }
            Kokkos::fence("BarnesHutAlgorithm::traceMemory");
            char entry[64];
            std::snprintf(entry, sizeof(entry), " %s=%.2f", stage, deviceUsedGB());
            memoryTrace_m += entry;
        }
        bool imagesWereActive_m    = false;
        bool imagesWereFar_m       = false;
        bool warnedRemoteMac_m     = false;
        double softening_m         = 0.0;  ///< Uniform h of the current solve.
        Statistics statistics_m;

    private:
        void checkPreconditions(const SpaceChargeSolveContext& context) const;
        [[nodiscard]] bool imagesActive(std::size_t step) const;
        void ensureContainer(const std::array<double, 6>& bounds, unsigned long long total,
                             bool newDomain = false);
        void updateDomainBox(const std::array<double, 6>& bounds, unsigned long long total);
        void ensureSolver(bool images);
        BHLocalIndex prepareSlots(size_type nSlots);
        void pack(
                BHLocalIndex offset, size_type nLocal, double timeStep, bool images, double planeZ,
                ZStretch stretch, std::optional<double> referenceRapidity);
        void reportRemoteMacFailures(std::size_t step);
        void synchronize(bool newDomain);
        void exchangeChargeHalos();
        void traverse();
        SlotField imageOnlyTraversal();
        void copyBack(ZStretch stretch, const SlotField& allSources);
        void verifySolve(
                const SpaceChargeSolveContext& context, size_type totalBefore, bool images,
                double gamma);
    };

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::checkPreconditions(
            const SpaceChargeSolveContext& context) const {
        const char* where = "BarnesHutAlgorithm::solve";
        if (bunchState_m->fixedCartesianDomain().has_value()) {
            throw OpalException(
                    where, "TYPE=BH does not support a fixed Cartesian domain (BeamBeam).");
        }
        const double timeStep = context.stepState().timeStep;
        if (!std::isfinite(timeStep) || timeStep == 0.0) {
            throw OpalException(where, "Barnes-Hut charge weighting needs a nonzero time step.");
        }
        if (primary_m->getChargePerParticle() == 0.0) {
            throw OpalException(
                    where,
                    "Per-particle charge is zero but a space-charge solver is active (type=BH). "
                    "This almost always means the BEAM command is missing BCHARGE. Set BCHARGE on "
                    "the BEAM definition, or use TYPE=NONE when no space charge is intended.");
        }

        // Particle deletion runs after the space-charge step, so no particle is pending deletion
        // here. The copy-back resets InvalidMask, so a pending mark would be silently lost.
        const size_type localMarked = countMarkedForDeletion(*primary_m);
        size_type globalMarked      = 0;
        ippl::Comm->allreduce(localMarked, globalMarked, 1, std::plus<size_type>());
        if (globalMarked != 0) {
            throw OpalException(
                    where, "TYPE=BH requires that no particle is marked for deletion when the "
                           "space-charge step starts.");
        }
    }

    template <class Precision>
    bool BarnesHutAlgorithm::PrecisionImpl<Precision>::imagesActive(std::size_t step) const {
        // Same expiry rule as CartesianPIC3D (ZEROFACE_MAXSTEPS; zero never expires).
        const DirichletPlaneConfig& plane = config_m.dirichletPlane;
        return plane.enabled() && !(plane.maximumSteps != 0 && step >= plane.maximumSteps);
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::ensureContainer(
            const std::array<double, 6>& bounds, unsigned long long total, bool newDomain) {
        // total counts every Barnes-Hut source, including image charges.
        const int nRanks       = ippl::Comm->size();
        const auto ranks       = static_cast<unsigned long long>(nRanks);
        const auto perRank     = std::max<unsigned long long>(1, total / ranks);

        // ryoanji starts every traversal at node 1, the first child of the root, so the focus tree
        // must not consist of the root leaf alone: keep the focus bucket well below the per-rank
        // particle count.
        const auto configuredFocus = static_cast<unsigned>(config_m.bucketSizeFocus);
        const unsigned bucketSizeFocus = static_cast<unsigned>(std::min<unsigned long long>(
                configuredFocus, std::max<unsigned long long>(1, perRank / 8)));

        // cstone assigns whole global-tree leaves to ranks; with fewer leaves than ranks some ranks
        // own nothing and the first sync's focus-tree convergence loop never terminates. Small
        // bunches therefore get a global bucket of at most a quarter of the per-rank count.
        unsigned bucketSize =
                config_m.bucketSize != 0
                        ? static_cast<unsigned>(config_m.bucketSize)
                        : static_cast<unsigned>(std::min<unsigned long long>(
                                std::max<unsigned long long>(kMinimumBucketSize, total / (100ull * ranks)),
                                std::max<unsigned long long>(1, perRank / 4)));
        bucketSize = std::max(bucketSize, bucketSizeFocus);  // cstone: bucketSize >= focus

        // Bucket sizes are fixed per cstone domain. Emission can grow the bunch by orders of
        // magnitude and losses can shrink it; rebuild when the sizes are far off or the trees would
        // degenerate. All particles live in the OPALX container between solves, so a rebuild
        // loses no state. newDomain: updateDomainBox() resizes the box.
        const bool bucketsOff =
                bh_m != nullptr
                && ((config_m.bucketSize == 0
                     && (bucketSize > 4 * bucketSize_m || 4 * bucketSize < bucketSize_m))
                    || total <= bucketSizeFocus_m
                    || (bucketSizeFocus != bucketSizeFocus_m
                        && (bucketSizeFocus >= 2 * bucketSizeFocus_m
                            || bucketSizeFocus == configuredFocus)));
        if (bh_m != nullptr && !bucketsOff && !newDomain) {
            return;
        }

        static IpplTimings::TimerRef buildTimer = IpplTimings::getTimer("bh.newContainer");
        IpplTimings::startTimer(buildTimer);
        if (bucketsOff) {
            ++statistics_m.bucketRebuilds;
        }
        // Fixed boundaries: cstone keeps the box that updateDomainBox() sets before every sync
        // (an open box would be refitted without the guard against zero-extent axes).
        const cstone::Box<double> box = domainBox(bounds, config_m.boxMargin);
        if (bh_m != nullptr) {
            // A new cstone domain in the existing container: its buffers, and the solver with
            // its multipole buffers, are kept. Freeing them grew the device memory at every
            // rebuild (GPU-aware MPI retains freed buffers it mapped, see resetDomain()).
            bh_m->resetDomain(bucketSize, bucketSizeFocus, static_cast<float>(config_m.theta), box);
        } else {
            bh_m = std::make_unique<BHContainer>(
                    ippl::Comm->rank(), nRanks, bucketSize, bucketSizeFocus,
                    static_cast<float>(config_m.theta),
                    std::array<double, 6>{box.xmin(), box.xmax(), box.ymin(), box.ymax(),
                                          box.zmin(), box.zmax()},
                    std::array<cstone::BoundaryType, 3>{
                            box.boundaryX(), box.boundaryY(), box.boundaryZ()},
                    ippl::Comm->getCommunicator());
            bh_m->setUniformH(softening_m);
            bh_m->setLeafBasedH(config_m.leafBasedSoftening);
        }

        bucketSize_m      = bucketSize;
        bucketSizeFocus_m = bucketSizeFocus;
        ++statistics_m.containerBuilds;
        synchronizeDevice();
        IpplTimings::stopTimer(buildTimer);
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::updateDomainBox(
            const std::array<double, 6>& bounds, unsigned long long total) {
        // cstone updates its trees incrementally, which is only valid while the box keeps its
        // size: the box follows the sources by moving, and a resize rebuilds the domain (a resized
        // box in the old domain left ranks without particles and hung the traversal).
        constexpr unsigned maxLevel        = cstone::maxTreeLevel<BHKeyType>{};
        const cstone::Box<double> fitted   = domainBox(bounds, config_m.boxMargin);
        const cstone::Box<double> previous = bh_m->box();
        const double lower[3]    = {previous.xmin(), previous.ymin(), previous.zmin()};
        const double edge[3]     = {previous.lx(), previous.ly(), previous.lz()};
        const double refitted[3] = {fitted.lx(), fitted.ly(), fitted.lz()};
        bool outgrown = false, oversized = false;
        for (unsigned d = 0; d < 3; ++d) {
            const double extent = bounds[2 * d + 1] - bounds[2 * d];
            outgrown  = outgrown || extent > edge[d];
            oversized = oversized || edge[d] > (1.0 + config_m.boxSlack) * refitted[d];
        }

        if (outgrown || oversized) {
            traceMemory("rebuild");
            ensureContainer(bounds, total, /*newDomain=*/true);  // with domainBox(bounds, boxMargin)
            traceMemory("built");
            ++statistics_m.boxRebuilds;
            statistics_m.boxShrinkRebuilds += outgrown ? 0 : 1;
        } else {
            // Recentre at every solve so the sources keep their keys: cstone refines by one level
            // per sync, and a bunch drifting through a resting box left overfull leaves (several
            // times the P2P work at BHBOXMARGIN 0.5).
            // Skip shifts below the finest key cell (rounding of an unchanged centre).
            const auto            bitsBefore = previous.getBoxDimBits(maxLevel);
            std::array<double, 6> moved{};
            bool                  shifted = false;
            for (unsigned d = 0; d < 3; ++d) {
                const double centre = 0.5 * (bounds[2 * d] + bounds[2 * d + 1]);
                moved[2 * d]        = centre - 0.5 * edge[d];
                moved[2 * d + 1]    = centre + 0.5 * edge[d];
                const double finestCell = std::ldexp(edge[d], -static_cast<int>(bitsBefore[d]));
                shifted = shifted || std::abs(moved[2 * d] - lower[d]) > finestCell;
            }
            if (shifted) {
                using cstone::BoundaryType;
                bh_m->setBox(cstone::Box<double>(
                        moved[0], moved[1], moved[2], moved[3], moved[4], moved[5],
                        BoundaryType::fixed, BoundaryType::fixed, BoundaryType::fixed));
                ++statistics_m.boxMoves;
            }
        }

        const auto bits      = bh_m->box().getBoxDimBits(maxLevel);
        statistics_m.keyBits = {bits[0], bits[1], bits[2]};
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::ensureSolver(bool images) {
        // Cells that mix real particles and image charges have a dipole about their
        // |q|-weighted expansion center, which quadrupole-only multipoles drop (first-order
        // error). Image-free solves keep the quadrupole type, so their results do not change.
        // Per-particle velocities: HeavisideQpole, which includes the dipole.
        using ippl::nbody::MultipoleOrder;
        const MultipoleOrder wanted =
                config_m.kinematics == BarnesHutKinematics::Particle ? MultipoleOrder::Heaviside
                : images ? MultipoleOrder::DipoleQuadrupole
                         : MultipoleOrder::Quadrupole;
        // The solver refers to the container, which lives as long as this algorithm: a domain
        // rebuild keeps both.
        if (solver_m != nullptr && solverOrder_m == wanted) {
            return;
        }
        typename BHSolver::Params params;
        // ryoanji accumulates a_i = G sum_j m_j (r_j - r_i) / r^3. Real particles are packed with
        // m_j = |q_j| and the sign of the common charge is applied in the copy-back; image charges
        // carry -|q_j|. Mixed-sign sources need ryoanji's remote-leaf fallback to recognise any
        // leaf with sources, also one whose net charge is zero.
        params.G          = -1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
        params.theta      = static_cast<float>(config_m.theta);
        params.numShells  = 0;
        params.multipoles = wanted;
        // --info 3: per-solve consistency checks of the sources ([bh.check]); a failed traversal
        // reports them regardless.
        params.checks = ippl::Info->getOutputLevel() >= 3;
        static IpplTimings::TimerRef solverTimer = IpplTimings::getTimer("bh.newSolver");
        IpplTimings::startTimer(solverTimer);
        solver_m.reset();
        solver_m      = std::make_unique<BHSolver>(*bh_m, params);
        solverOrder_m = wanted;
        synchronizeDevice();
        IpplTimings::stopTimer(solverTimer);
    }

    template <class Precision>
    auto BarnesHutAlgorithm::PrecisionImpl<Precision>::prepareSlots(size_type nSlots)
            -> BHLocalIndex {
        if (nSlots > static_cast<size_type>(std::numeric_limits<BHLocalIndex>::max() / 2)) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve",
                    "Too many local particles for the 32-bit Barnes-Hut indices.");
        }
        const auto n = static_cast<BHLocalIndex>(nSlots);

        // After a copy-back the owned slots already belong to this rank's particle set; when the
        // count is unchanged they are overwritten in place and the next sync only moves particles
        // that crossed a space-filling-curve boundary.
        if (statistics_m.solves > 0 && bh_m->getLocalNum() == n && !config_m.repackEverySolve) {
            ++statistics_m.inPlacePacks;
            return bh_m->startIndex();
        }

        // Otherwise drop the previous owned set and append the current one.
        const BHLocalIndex owned = bh_m->getLocalNum();
        if (owned > 0) {
            Kokkos::View<bool*> all("BarnesHutAlgorithm::destroyMask", owned);
            Kokkos::deep_copy(all, true);
            Kokkos::fence();
            bh_m->destroy(all.data(), owned);
            synchronizeDevice();
        }
        return bh_m->create(n);
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::pack(
            BHLocalIndex offset, size_type nLocal, double timeStep, bool images, double planeZ,
            ZStretch stretch, std::optional<double> referenceRapidity) {
        // Deposition weight of CartesianPIC3D: q_i = Q dt_i / dt_step (emission fractions).
        const double chargeScale = std::abs(primary_m->getChargePerParticle()) / timeStep;
        packParticles(
                *bh_m, *primary_m, offset, nLocal, chargeScale, softening_m, images, planeZ,
                stretch, referenceRapidity);
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::reportRemoteMacFailures(std::size_t step) {
        // cstone builds the locally essential tree with the unscaled MAC; a remote leaf that fails
        // the rapidity-scaled MAC is applied as a multipole anyway (0 on one rank).
        const unsigned long long failures = solver_m->remoteMacFailures();
        statistics_m.remoteMacFailures += failures;
        if (failures > 0 && !warnedRemoteMac_m) {
            Inform m("BarnesHutAlgorithm");
            m << level1 << "WARNING: TYPE=BH with per-particle velocities on several ranks: "
              << failures << " remote octree leaves at step " << step << " were approximated "
              << "although their rapidity-scaled multipole acceptance criterion failed (the "
              << "locally essential tree is not yet velocity-aware). Counted in the run summary."
              << endl;
            warnedRemoteMac_m = true;
        }
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::synchronize(bool newDomain) {
        // A new domain converges its trees from scratch; its syncs are timed separately.
        static IpplTimings::TimerRef keptTimer = IpplTimings::getTimer("bh.syncGrav");
        static IpplTimings::TimerRef newTimer  = IpplTimings::getTimer("bh.syncGrav.new");
        static IpplTimings::TimerRef haloTimer = IpplTimings::getTimer("bh.haloCharge");
        const IpplTimings::TimerRef syncTimer  = newDomain ? newTimer : keptTimer;

        IpplTimings::startTimer(syncTimer);
        bh_m->setUniformH(softening_m);  // BHSOFTENINGREL changes it every solve
        if (config_m.kinematics == BarnesHutKinematics::Particle) {
            ippl::nbody::syncGravBH<Precision, RelConserved, RelDependent>(*bh_m);
        } else {
            ippl::nbody::syncGravBH<Precision, Conserved, Dependent>(*bh_m);
        }
        if (!config_m.leafBasedSoftening) {
            // Halo slots that cstone grew or reused must carry the uniform softening too.
            Kokkos::deep_copy(ippl::nbody::getView<"h">(*bh_m), softening_m);
        }
        synchronizeDevice();
        IpplTimings::stopTimer(syncTimer);

        IpplTimings::startTimer(haloTimer);
        exchangeChargeHalos();
        IpplTimings::stopTimer(haloTimer);
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::exchangeChargeHalos() {
        bh_m->exchangeHalos(
                std::tie(bh_m->template getDV<"charge">()), bh_m->haloSendBuf(),
                bh_m->haloRecvBuf());
        synchronizeDevice();
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::traverse() {
        try {
            // warmup = true skips the per-solve interaction counts and stage timers (BHSTATS).
            solver_m->runSolver(/*warmup=*/!config_m.traversalStatistics);
        } catch (const std::runtime_error& error) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve",
                    std::string(error.what())
                            + (newDomainSolve_m ? " This solve built a new domain."
                                                : " This solve kept the previous domain.")
                            + " Increase FIELDSOLVER BHTHETA to shorten the traversal.");
        }
        synchronizeDevice();
    }

    template <class Precision>
    SlotField BarnesHutAlgorithm::PrecisionImpl<Precision>::imageOnlyTraversal() {
        // Keep the field of all sources, then traverse with the real charges zeroed to get the
        // image field. The next solve repacks every slot, so nothing needs restoring.
        SlotField allSources = keepFieldAndClearRealCharges(*bh_m);
        exchangeChargeHalos();
        traverse();
        ++statistics_m.imagePasses;
        return allSources;
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::copyBack(
            ZStretch stretch, const SlotField& allSources) {
        // Image charges stay behind; only the real particles return to OPALX.
        const size_type nOwned = countOwnedReal(*bh_m);
        // Collective: every rank resizes, and the global count stays unchanged.
        primary_m->replaceLocalCount(nOwned);
        const double sign = primary_m->getChargePerParticle() < 0.0 ? -1.0 : 1.0;
        copyBackParticles(
                *bh_m, *primary_m, sign, stretch, allSources,
                config_m.kinematics == BarnesHutKinematics::Particle);
        primary_m->markMomentsDirty();
    }

    template <class Precision>
    SpaceChargeSolveResult BarnesHutAlgorithm::PrecisionImpl<Precision>::solve(
            const SpaceChargeSolveContext& context) {
        static IpplTimings::TimerRef domainTimer   = IpplTimings::getTimer("bh.domain");
        static IpplTimings::TimerRef packTimer     = IpplTimings::getTimer("bh.pack");
        static IpplTimings::TimerRef solveTimer    = IpplTimings::getTimer("bh.solve");
        static IpplTimings::TimerRef copyBackTimer = IpplTimings::getTimer("bh.copyBack");
        // Solves on a domain built in the same call (first solve, rebuilds) are timed separately.
        static IpplTimings::TimerRef packNewTimer     = IpplTimings::getTimer("bh.pack.new");
        static IpplTimings::TimerRef solveNewTimer    = IpplTimings::getTimer("bh.solve.new");
        static IpplTimings::TimerRef copyBackNewTimer = IpplTimings::getTimer("bh.copyBack.new");

        SpaceChargeSolveResult result;
        // Equivalent to clearSelfFields() without ParticleAttrib::operator=, whose extended-lambda
        // wrapper can be miscompiled across TUs (nvcc helper statics + GCC IPA-SRA clones).
        Kokkos::deep_copy(primary_m->E.getView(), Vector_t<double, 3>(0.0));
        Kokkos::deep_copy(primary_m->B.getView(), Vector_t<double, 3>(0.0));
        // The global count handles early-emission ranks that own no particles.
        if (primary_m->getTotalNum() <= 1) {
            return result;
        }
        checkPreconditions(context);
        const size_type totalBefore = primary_m->getTotalNum();
        std::optional<SolveTimeReport> timeReport;
        double devStart = 0.0;
        if (config_m.traversalStatistics) {
            timeReport.emplace();
            devStart = deviceUsedGB();
        }

        enterSolveFrame(context.stepState().frames, *primary_m);

        // Rest frame (BINS with MAXBINS=1): z stretched by the mean gamma about the centre of
        // the bunch, unchanged charges (the rho / gamma of CartesianPIC3D), boost afterwards.
        // Per-particle velocities use the same stretched frame as tree coordinates; every source
        // then carries its own gamma relative to it.
        static IpplTimings::TimerRef restFrameTimer = IpplTimings::getTimer("bh.restFrame");
        const std::array<double, 6> labBounds = globalBounds(*primary_m);
        const bool relativistic       = config_m.kinematics != BarnesHutKinematics::Static;
        const bool particleKinematics = config_m.kinematics == BarnesHutKinematics::Particle;
        RestFrameBoost boost;
        ZStretch stretch;
        if (relativistic) {
            IpplTimings::startTimer(restFrameTimer);
            boost   = restFrameOf(primary_m->P.getView(), primary_m->getLocalNum());
            stretch = ZStretch{0.5 * (labBounds[4] + labBounds[5]), boost.gamma};
            IpplTimings::stopTimer(restFrameTimer);
            statistics_m.gamma = boost.gamma;
        }
        std::array<double, 6> bounds = labBounds;
        bounds[4]                    = stretch.apply(labBounds[4]);
        bounds[5]                    = stretch.apply(labBounds[5]);
        const double planeZ          = stretch.apply(config_m.dirichletPlane.planeZ);

        const double extent = std::max({bounds[1] - bounds[0], bounds[3] - bounds[2],
                                        bounds[5] - bounds[4]});
        if (!std::isfinite(extent)) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve", "Particle positions are not finite.");
        }
        const auto total = static_cast<unsigned long long>(primary_m->getTotalNum());
        softening_m      = config_m.softening;
        if (config_m.relativeSoftening > 0.0 && extent > 0.0) {
            // The density changes by orders of magnitude along a beamline; h follows it.
            softening_m = config_m.relativeSoftening
                          * meanSpacing(*primary_m, stretch, bounds, total);
        }
        statistics_m.softening = softening_m;
        const bool images = imagesActive(context.stepState().step);
        if (config_m.dirichletPlane.enabled() && images != imagesWereActive_m) {
            Inform m("BarnesHutAlgorithm");
            m << level2 << "TYPE=BH image charges of the plane z = "
              << config_m.dirichletPlane.planeZ << " m are " << (images ? "active" : "expired")
              << " at step " << context.stepState().step << "." << endl;
            imagesWereActive_m = images;
        }
        // In the rest frame the images recede with 2 gamma times the distance to the plane; far
        // ones are one quadrupole instead of sources, which keeps the box around the bunch. The
        // per-particle direct sum always takes the images as sources.
        // Per-particle velocities: the images become one multipole with rapidity moments once
        // their extent passes the MAC scaled by kappa (heavisideMacScale2) and their rapidity
        // spread passes BHTHETAV.
        const ryoanji::HeavisideFrame<double> heavisideFrame =
                particleKinematics ? heavisideFrameOf(boost) : ryoanji::HeavisideFrame<double>{};
        const double chargeScale =
                std::abs(primary_m->getChargePerParticle()) / context.stepState().timeStep;
        FarImageMultipole farImageField;
        bool farImages = images && relativistic && extent > 0.0
                         && imagesAreFar(bounds, planeZ, config_m.theta);
        if (farImages && particleKinematics) {
            farImageField = farImageMultipole(
                    *primary_m, stretch, bounds, planeZ, chargeScale, heavisideFrame,
                    config_m.thetaV, config_m.velocityMoments);
            farImages = !farImageField.empty && farImageField.macScale2 > 0.0
                        && imagesAreFar(bounds, planeZ,
                                        config_m.theta / std::sqrt(farImageField.macScale2));
        }
        if (images && farImages != imagesWereFar_m) {
            Inform m("BarnesHutAlgorithm");
            m << level2 << "TYPE=BH evaluates the image charges "
              << (farImages ? "as one far-field quadrupole" : "as Barnes-Hut sources")
              << " from step " << context.stepState().step << " (gamma = " << boost.gamma << ")."
              << endl;
            imagesWereFar_m = farImages;
        }
        const bool sourceImages = images && !farImages;
        // Image charges move against the bunch: the boost needs their field separately. With
        // per-particle velocities every image carries its own (reversed) velocity instead.
        const bool imageField = sourceImages && boost.moving() && !particleKinematics;
        const double sign = primary_m->getChargePerParticle() < 0.0 ? -1.0 : 1.0;

        const auto directLimit =
                kDirectSumParticlesPerRank * static_cast<unsigned long long>(ippl::Comm->size());
        // BHDIRECT forces the direct sum (reference runs); particles then keep the OPALX layout.
        const bool forcedDirect = config_m.directSum;
        const bool useDirectSum = extent > 0.0 && (forcedDirect || total < directLimit);
        if (useDirectSum != usingDirectSum_m && extent > 0.0) {
            Inform m("BarnesHutAlgorithm");
            if (forcedDirect) {
                m << level2 << "TYPE=BH evaluates the field by direct O(N^2) summation "
                  << "(BHDIRECT=TRUE)." << endl;
            } else if (useDirectSum) {
                m << level1 << "WARNING: TYPE=BH uses a direct O(N^2) sum instead of Barnes-Hut: "
                  << total << " particles are fewer than " << kDirectSumParticlesPerRank
                  << " per rank (" << directLimit << " on " << ippl::Comm->size()
                  << " ranks). Switching to Barnes-Hut once the bunch is larger." << endl;
            } else {
                m << level2 << "TYPE=BH switches from the direct sum to Barnes-Hut at " << total
                  << " particles." << endl;
            }
            usingDirectSum_m = useDirectSum;
        }

        double fieldGamma = boost.gamma;  // bound factor of verifySolve
        const double prefactor = -sign / (4.0 * Physics::pi * Physics::epsilon_0);
        if (useDirectSum) {
            static IpplTimings::TimerRef directTimer = IpplTimings::getTimer("bh.directSum");
            IpplTimings::startTimer(directTimer);
            // Leaf-based softening has no tree here; a negligible floor only defuses the self pair.
            const double softening = softening_m > 0.0 ? softening_m : 1.0e-9 * extent;
            if (particleKinematics) {
                const double sourceGamma = directSumHeavisideFields(
                        *primary_m, chargeScale, softening, prefactor, sourceImages, planeZ,
                        stretch, heavisideFrame);
                fieldGamma = std::max(fieldGamma, sourceGamma);
                ++statistics_m.particleSolves;
                statistics_m.maxGamma = sourceGamma;
            } else {
                directSumFields(
                        *primary_m, chargeScale, softening, prefactor, sourceImages, planeZ,
                        stretch, imageField);
            }
            IpplTimings::stopTimer(directTimer);
            ++statistics_m.directSolves;
            statistics_m.imagePasses += imageField ? 1 : 0;
            result.backendSolves = 1;
        } else if (extent > 0.0) {
            IpplTimings::startTimer(domainTimer);
            const std::size_t buildsBefore     = statistics_m.containerBuilds;
            std::array<double, 6> sourceBounds = bounds;
            if (sourceImages) {
                // The mirrored bunch extends the source region beyond the plane.
                sourceBounds[4] = std::min(bounds[4], 2.0 * planeZ - bounds[5]);
                sourceBounds[5] = std::max(bounds[5], 2.0 * planeZ - bounds[4]);
            }
            const unsigned long long sources = sourceImages ? 2 * total : total;
            ensureContainer(sourceBounds, sources);
            updateDomainBox(sourceBounds, sources);
            ensureSolver(sourceImages);
            if (particleKinematics) {
                solver_m->setHeavisideParams(
                        {heavisideFrame, config_m.thetaV, config_m.velocityMoments});
            }
            IpplTimings::stopTimer(domainTimer);
            const bool newDomain = statistics_m.containerBuilds != buildsBefore;
            newDomainSolve_m     = newDomain;

            IpplTimings::startTimer(newDomain ? packNewTimer : packTimer);
            const size_type nLocal    = primary_m->getLocalNum();
            const BHLocalIndex offset = prepareSlots(sourceImages ? 2 * nLocal : nLocal);
            pack(offset, nLocal, context.stepState().timeStep, sourceImages, planeZ, stretch,
                 particleKinematics ? std::optional<double>(heavisideFrame.wRef) : std::nullopt);
            IpplTimings::stopTimer(newDomain ? packNewTimer : packTimer);

            traceMemory("packed");
            synchronize(newDomain);
            traceMemory("synced");

            IpplTimings::startTimer(newDomain ? solveNewTimer : solveTimer);
            traverse();
            IpplTimings::stopTimer(newDomain ? solveNewTimer : solveTimer);
            traceMemory("traversed");

            SlotField allSources;
            if (imageField) {
                static IpplTimings::TimerRef imageTimer = IpplTimings::getTimer("bh.imagePass");
                IpplTimings::startTimer(imageTimer);
                allSources = imageOnlyTraversal();
                IpplTimings::stopTimer(imageTimer);
            }

            IpplTimings::startTimer(newDomain ? copyBackNewTimer : copyBackTimer);
            copyBack(stretch, allSources);
            IpplTimings::stopTimer(newDomain ? copyBackNewTimer : copyBackTimer);
            traceMemory("copied");

            ++statistics_m.solves;
            if (solverOrder_m == ippl::nbody::MultipoleOrder::DipoleQuadrupole) {
                ++statistics_m.dipoleSolves;
            }
            if (particleKinematics) {
                const double sourceGamma = maxLongitudinalGamma(*primary_m);
                fieldGamma               = std::max(fieldGamma, sourceGamma);
                ++statistics_m.particleSolves;
                statistics_m.maxGamma = sourceGamma;
                reportRemoteMacFailures(context.stepState().step);
            }
            result.backendSolves = 1;
        }
        // All particles coincide otherwise: the softened field vanishes and E stays zero.

        if (result.backendSolves > 0) {
            if (farImages) {
                static IpplTimings::TimerRef farTimer = IpplTimings::getTimer("bh.farImages");
                IpplTimings::startTimer(farTimer);
                if (particleKinematics) {
                    addFarImageFieldHeaviside(*primary_m, stretch, farImageField, prefactor);
                } else {
                    addFarImageField(*primary_m, stretch, bounds, planeZ, chargeScale, sign);
                }
                IpplTimings::stopTimer(farTimer);
                ++statistics_m.farImageSolves;
            }
            if (particleKinematics) {
                IpplTimings::startTimer(restFrameTimer);
                boostParticleKinematics(*primary_m, heavisideFrame);
                IpplTimings::stopTimer(restFrameTimer);
            } else if (relativistic) {
                IpplTimings::startTimer(restFrameTimer);
                boostToLab(*primary_m, boost);
                IpplTimings::stopTimer(restFrameTimer);
                ++statistics_m.restFrameSolves;
            }
            verifySolve(context, totalBefore, images, fieldGamma);
        }

        leaveSolveFrame(context.stepState().frames, *primary_m);
        if (timeReport) {
            timeReport->print(static_cast<unsigned long long>(context.stepState().step), total,
                              boost.gamma, softening_m, images);
            std::size_t opalBytes = 0;
            primary_m->forAllAttributes([&opalBytes]<typename Attribute>(Attribute& attribute) {
                opalBytes += attribute->packedSize(attribute->size());
            });
            printMemoryReport(static_cast<unsigned long long>(context.stepState().step), devStart,
                              primary_m->getLocalNum(), primary_m->R.size(), opalBytes,
                              bh_m ? bh_m->nWithHalos() : 0, bh_m ? bh_m->reservedBytes() : 0);
            if (ippl::Comm->size() <= 16 && !memoryTrace_m.empty()) {
                std::printf("[bh.mem.trace] step=%llu rank=%d%s GB\n",
                            static_cast<unsigned long long>(context.stepState().step),
                            ippl::Comm->rank(), memoryTrace_m.c_str());
                std::fflush(stdout);
            }
        }
        memoryTrace_m.clear();
        return result;
    }

    template <class Precision>
    void BarnesHutAlgorithm::PrecisionImpl<Precision>::verifySolve(
            const SpaceChargeSolveContext& context, size_type totalBefore, bool images,
            double gamma) {
        // A failed halo/layout construction in cstone does not raise an error; it shows up as lost
        // or duplicated particles and absurd fields. Stop the run instead of tracking garbage.
        const char* where = "BarnesHutAlgorithm::solve";
        if (primary_m->getTotalNum() != totalBefore) {
            throw OpalException(
                    where, "The Barnes-Hut redistribution changed the particle count from "
                                   + std::to_string(totalBefore) + " to "
                                   + std::to_string(primary_m->getTotalNum())
                                   + ". This indicates an inconsistent domain decomposition.");
        }
        const FieldStatistics stats = globalFieldStatistics(*primary_m);
        if (stats.nonFinite > 0.0) {
            throw OpalException(
                    where, "Barnes-Hut produced non-finite fields for "
                                   + std::to_string(static_cast<long long>(stats.nonFinite))
                                   + " particles.");
        }
        // No pair can exceed k q_j / (2h)^2 with the softened kernel, so k sum|q| / (2h)^2 bounds
        // Bound for the exact sum (multipole errors are far smaller), times gamma for the boost;
        // only for a uniform softening.
        if (!config_m.leafBasedSoftening && softening_m > 0.0) {
            // Image charges double the source charge.
            const double totalCharge = (images ? 2.0 : 1.0)
                                       * std::abs(primary_m->getChargePerParticle()) * stats.sumDt
                                       / context.stepState().timeStep;
            const double coulomb     = 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
            const double bound = 2.0 * gamma * coulomb * totalCharge
                                 / (4.0 * softening_m * softening_m);
            if (stats.maxField > bound) {
                std::ostringstream message;
                message << "Barnes-Hut field |E| = " << stats.maxField
                        << " V/m exceeds the softened Coulomb bound " << bound
                        << " V/m. This indicates an inconsistent domain decomposition.";
                throw OpalException(where, message.str());
            }
        }
    }

    BarnesHutAlgorithm::BarnesHutAlgorithm(
            BarnesHutConfig config, std::span<ParticleContainer* const> particles,
            std::shared_ptr<const BunchStateHandler> bunchState) {
        const char* where = "BarnesHutAlgorithm::BarnesHutAlgorithm";
        validateSpaceChargeConfig(SpaceChargeConfig(config));
        if (bunchState == nullptr) {
            throw OpalException(where, "The bunch state handler is null.");
        }
        ParticleContainer& primary = requirePrimary(particles);
        if (primary.getQMStorageMode() != ParticleContainer::QMStorageMode::SingleValue) {
            throw OpalException(where, "TYPE=BH requires QM_MODE=SINGLE.");
        }
        if (primary.hasSpin()) {
            throw OpalException(where, "TYPE=BH does not yet carry spin (POLARIZATION).");
        }
        if (primary.getAttributeNum() != kSupportedAttributeCount) {
            throw OpalException(
                    where, "TYPE=BH migrates particles and must carry every attribute; the primary "
                           "container has attributes the Barnes-Hut copy-back does not know.");
        }
        static_assert(
                sizeof(typename ParticleContainer::index_type) == sizeof(std::uint64_t),
                "Particle IDs are carried through the 64-bit Barnes-Hut payload slot.");

        primary.setDecompositionOwnedExternally(true);
        switch (config.precision) {
            case BarnesHutPrecision::Mixed:
                impl_m = std::make_unique<PrecisionImpl<ippl::nbody::MixedPrecision>>(
                        std::move(config), primary, std::move(bunchState));
                break;
            case BarnesHutPrecision::Double:
            default:
                impl_m = std::make_unique<PrecisionImpl<ippl::nbody::DoublePrecision>>(
                        std::move(config), primary, std::move(bunchState));
                break;
        }
    }

    BarnesHutAlgorithm::~BarnesHutAlgorithm() = default;

    namespace {
        // An exception thrown on some ranks only (a failed GPU allocation, an exhausted traversal
        // stack) leaves the other ranks in the next collective, and unwinding then blocks in the
        // collective MPI_Win_free of the particle layout: the run hangs without a message, which
        // only rank 0 would print. With several ranks, report the error on the rank that saw it
        // and abort the job instead.
        [[noreturn]] void abortOnRank(const std::string& where, const std::string& what) {
            int rank = -1;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            std::fprintf(stderr, "*** TYPE=BH rank %d: %s: %s\n*** Aborting all ranks.\n", rank,
                         where.c_str(), what.c_str());
            std::fflush(stderr);
            MPI_Abort(MPI_COMM_WORLD, 1);
            std::abort();
        }
    }  // namespace

    SpaceChargeSolveResult BarnesHutAlgorithm::solve(const SpaceChargeSolveContext& context) {
        int ranks = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        if (ranks == 1) {
            return impl_m->solve(context);
        }
        try {
            return impl_m->solve(context);
        } catch (const OpalException& error) {
            abortOnRank(error.where(), error.what());
        } catch (const std::exception& error) {
            abortOnRank("BarnesHutAlgorithm::solve", error.what());
        }
    }

    const BarnesHutAlgorithm::Statistics& BarnesHutAlgorithm::statistics() const {
        return impl_m->statistics();
    }

    std::string BarnesHutAlgorithm::runSummary() const {
        const Statistics& s = impl_m->statistics();
        std::ostringstream out;
        out << "TYPE=BH: " << s.solves << " tree solves (" << s.dipoleSolves << " with dipoles), "
            << s.directSolves << " direct sums, " << s.particleSolves
            << " solves with per-particle velocities (max gamma " << s.maxGamma << ", "
            << s.remoteMacFailures << " remote leaves failing their rapidity MAC), "
            << s.restFrameSolves << " rest-frame boosts, "
            << s.imagePasses << " image-only passes, " << s.farImageSolves
            << " far-image quadrupoles; " << s.containerBuilds << " domain builds = 1 initial + "
            << s.boxRebuilds - s.boxShrinkRebuilds << " outgrown box (BHBOXMARGIN "
            << impl_m->config().boxMargin << ") + " << s.boxShrinkRebuilds
            << " oversized box (BHBOXSLACK " << impl_m->config().boxSlack << ") + "
            << s.bucketRebuilds << " bucket size; " << s.boxMoves << " box moves.";
        return out.str();
    }

}  // namespace opalx::spacecharge
