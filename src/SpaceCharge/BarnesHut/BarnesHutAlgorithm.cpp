/**
 * @file BarnesHutAlgorithm.cpp
 * @brief Implements the Barnes-Hut space-charge algorithm on top of IPPL's NBody module.
 *
 * This is the only OPALX translation unit that includes NBody/cstone/ryoanji headers.
 */

#include "SpaceCharge/BarnesHut/BarnesHutAlgorithm.h"

#include "PartBunch/BunchStateHandler.h"
#include "SpaceCharge/SpaceChargeFrames.h"
#include "Utilities/OpalException.h"

#include "NBody/NBodyParticleContainer.hpp"
#include "NBody/NBodySolver.hpp"
#include "NBody/core/BHFieldLists.hpp"
#include "NBody/wrappers/NBodyKokkosView.hpp"
#include "ryoanji/interface/multipole_holder.cuh"
#include "ryoanji/nbody/traversal_cpu.hpp"

#include "Physics/Physics.h"
#include "Utility/IpplTimings.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mpi.h>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace opalx::spacecharge {
    namespace {

        using Precision    = ippl::nbody::DoublePrecision;
        using BHContainer  = ippl::nbody::NBodyParticleContainer<Precision, 3>;
        using BHSolver     = ippl::nbody::NBodySolver<Precision, 3>;
        using size_type    = ippl::detail::size_type;
        using Conserved    = ippl::nbody::fields::AuxConserved;
        using Dependent    = ippl::nbody::fields::StdDependent;
        using BHLocalIndex = BHContainer::LocalIndex;

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

        /// Wait for raw CUDA (cstone/ryoanji, default stream) and Kokkos work alike.
        void synchronizeDevice() {
            Kokkos::fence();
#if defined(USE_CUDA)
            syncGpu();
#endif
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
         * @brief Write the local OPALX particles into BH slots [offset, offset + nLocal); with
         * images, their mirrors (z' = 2 planeZ - z, opposite charge) follow in
         * [offset + nLocal, offset + 2 nLocal).
         */
        void packParticles(
                BHContainer& bh, BarnesHutAlgorithm::ParticleContainer& particles,
                BHLocalIndex offset, size_type nLocal, double chargeScale, double softening,
                bool images, double planeZ) {
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

            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::pack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const size_type j = off + i;
                        x(j)              = R(i)[0];
                        y(j)              = R(i)[1];
                        z(j)              = R(i)[2];
                        px(j)             = P(i)[0];
                        py(j)             = P(i)[1];
                        pz(j)             = P(i)[2];
                        m(j)              = chargeScale * dtOut(i);
                        h(j)              = softening;
                        dt(j)             = dtOut(i);
                        id(j)             = static_cast<std::uint64_t>(ID(i));
                        if (images) {
                            const size_type k = j + nLocal;
                            x(k)              = R(i)[0];
                            y(k)              = R(i)[1];
                            z(k)              = 2.0 * planeZ - R(i)[2];
                            px(k)             = P(i)[0];
                            py(k)             = P(i)[1];
                            pz(k)             = -P(i)[2];
                            m(k)              = -chargeScale * dtOut(i);
                            h(k)              = softening;
                            dt(k)             = dtOut(i);
                            id(k)             = static_cast<std::uint64_t>(ID(i)) | imageFlag;
                        }
                    });
            Kokkos::fence();
        }

        /** @brief Number of real (non-image) particles among the BH-owned slots. */
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

        /**
         * @brief Overwrite the (already resized) OPALX particles with the real particles of the
         * BH-owned set, in slot order; image charges are skipped.
         */
        void copyBackParticles(
                BHContainer& bh, BarnesHutAlgorithm::ParticleContainer& particles, double sign) {
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
                        R(k)              = Vector_t<double, 3>(x(j), y(j), z(j));
                        P(k)              = Vector_t<double, 3>(px(j), py(j), pz(j));
                        dtOut(k)          = dt(j);
                        ID(k)             = static_cast<index_type>(id(j));
                        E(k)       = Vector_t<double, 3>(sign * ax(j), sign * ay(j), sign * az(j));
                        B(k)       = Vector_t<double, 3>(0.0);
                        Phi(k)     = 0.0;
                        Bin(k)     = 0;
                        invalid(k) = false;
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
         * applied here, exactly as in the copy-back.
         */
        void directSumFields(
                BarnesHutAlgorithm::ParticleContainer& particles, double chargeScale,
                double softening, double prefactor, bool images, double planeZ) {
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
                        local(kValues * i)     = R(i)[0];
                        local(kValues * i + 1) = R(i)[1];
                        local(kValues * i + 2) = R(i)[2];
                        local(kValues * i + 3) = chargeScale * dtView(i);
                        local(kValues * i + 4) = softening;
                        if (images) {
                            const size_type k  = nLocal + i;
                            local(kValues * k)     = R(i)[0];
                            local(kValues * k + 1) = R(i)[1];
                            local(kValues * k + 2) = 2.0 * planeZ - R(i)[2];
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
            const ryoanji::Vec3<double> unusedBox{1.0, 1.0, 1.0};  // numShells = 0: open
#if defined(USE_CUDA)
            ryoanji::directSum<double>(
                    first, first + nLocal, nAll, unusedBox, 0, device[0].data(), device[1].data(),
                    device[2].data(), device[3].data(), device[4].data(), potential.data(),
                    ax.data(), ay.data(), az.data());
#else
            // The CPU variant evaluates every target; only this rank's slice is used below.
            ryoanji::directSum(
                    device[0].data(), device[1].data(), device[2].data(), device[4].data(),
                    device[3].data(), static_cast<cstone::LocalIndex>(nAll), 1.0f, unusedBox, 0,
                    ax.data(), ay.data(), az.data(), potential.data());
#endif
            synchronizeDevice();

            auto E = particles.E.getView();
            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::directUnpack", nLocal, KOKKOS_LAMBDA(const size_type i) {
                        const size_type j = first + i;
                        E(i) = Vector_t<double, 3>(
                                prefactor * ax(j), prefactor * ay(j), prefactor * az(j));
                    });
            Kokkos::fence();
        }

        /** @brief Local field and charge statistics of the primary after a solve. */
        struct FieldStatistics {
            double nonFinite = 0.0;  ///< Particles with a non-finite field component.
            double maxField  = 0.0;  ///< Largest |E| in V/m.
            double sumDt     = 0.0;  ///< Sum of dt_i; total source charge is |Q| sum(dt_i) / dt_step.
        };

        FieldStatistics globalFieldStatistics(BarnesHutAlgorithm::ParticleContainer& particles) {
            const auto E           = particles.E.getView();
            const auto dt          = particles.dt.getView();
            const size_type nLocal = particles.getLocalNum();
            double nonFinite = 0.0, maxE2 = 0.0, sumDt = 0.0;
            Kokkos::parallel_reduce(
                    "BarnesHutAlgorithm::fieldStatistics", nLocal,
                    KOKKOS_LAMBDA(const size_type i, double& bad, double& m, double& s) {
                        const double e2 = E(i)[0] * E(i)[0] + E(i)[1] * E(i)[1] + E(i)[2] * E(i)[2];
                        if (!Kokkos::isfinite(e2)) {
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

    }  // namespace

    class BarnesHutAlgorithm::Impl {
    public:
        Impl(BarnesHutConfig config, ParticleContainer& primary,
             std::shared_ptr<const BunchStateHandler> bunchState)
            : config_m(std::move(config)), primary_m(&primary), bunchState_m(std::move(bunchState)) {}

        SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context);

        BarnesHutConfig config_m;
        ParticleContainer* primary_m = nullptr;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<BHContainer> bh_m;
        std::unique_ptr<BHSolver> solver_m;
        ippl::nbody::MultipoleOrder solverOrder_m = ippl::nbody::MultipoleOrder::Quadrupole;
        unsigned bucketSize_m      = 0;
        unsigned bucketSizeFocus_m = 0;
        bool usingDirectSum_m      = false;
        bool imagesWereActive_m    = false;
        Statistics statistics_m;

    private:
        void checkPreconditions(const SpaceChargeSolveContext& context) const;
        [[nodiscard]] bool imagesActive(std::size_t step) const;
        void ensureContainer(const std::array<double, 6>& bounds, unsigned long long total);
        void ensureSolver(bool images);
        BHLocalIndex prepareSlots(size_type nSlots);
        void pack(BHLocalIndex offset, size_type nLocal, double timeStep, bool images);
        void synchronize();
        void copyBack();
        void verifySolve(
                const SpaceChargeSolveContext& context, size_type totalBefore, bool images);
    };

    void BarnesHutAlgorithm::Impl::checkPreconditions(const SpaceChargeSolveContext& context) const {
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

    bool BarnesHutAlgorithm::Impl::imagesActive(std::size_t step) const {
        // Same expiry rule as CartesianPIC3D (ZEROFACE_MAXSTEPS; zero never expires).
        const DirichletPlaneConfig& plane = config_m.dirichletPlane;
        return plane.enabled() && !(plane.maximumSteps != 0 && step >= plane.maximumSteps);
    }

    void BarnesHutAlgorithm::Impl::ensureContainer(
            const std::array<double, 6>& bounds, unsigned long long total) {
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
        // loses no state.
        const bool rebuild =
                bh_m != nullptr
                && ((config_m.bucketSize == 0
                     && (bucketSize > 4 * bucketSize_m || 4 * bucketSize < bucketSize_m))
                    || total <= bucketSizeFocus_m
                    || (bucketSizeFocus != bucketSizeFocus_m
                        && (bucketSizeFocus >= 2 * bucketSizeFocus_m
                            || bucketSizeFocus == configuredFocus)));
        if (bh_m != nullptr && !rebuild) {
            return;
        }

        solver_m.reset();
        bh_m.reset();
        // Same cube cstone's cubic_open refit produces; avoids zero-extent axes (flat bunches).
        const double side = std::max({bounds[1] - bounds[0], bounds[3] - bounds[2],
                                      bounds[5] - bounds[4]});
        const std::array<double, 6> box{bounds[0], bounds[0] + side, bounds[2], bounds[2] + side,
                                        bounds[4], bounds[4] + side};
        using cstone::BoundaryType;
        // cstone refits open boxes on every sync, but on the first sync a rank without particles
        // contributes the constructor box to the global extent; seed it with the true bounds.
        bh_m = std::make_unique<BHContainer>(
                ippl::Comm->rank(), nRanks, bucketSize, bucketSizeFocus,
                static_cast<float>(config_m.theta),
                box,
                std::array<BoundaryType, 3>{
                        BoundaryType::cubic_open, BoundaryType::cubic_open,
                        BoundaryType::cubic_open},
                ippl::Comm->getCommunicator());
        bh_m->setUniformH(config_m.softening);
        bh_m->setLeafBasedH(config_m.leafBasedSoftening);

        bucketSize_m      = bucketSize;
        bucketSizeFocus_m = bucketSizeFocus;
        ++statistics_m.containerBuilds;
    }

    void BarnesHutAlgorithm::Impl::ensureSolver(bool images) {
        // Cells that mix real particles and image charges have a dipole about their
        // |q|-weighted expansion center, which quadrupole-only multipoles drop (first-order
        // error). Image-free solves keep the quadrupole type, so their results do not change.
        using ippl::nbody::MultipoleOrder;
        const MultipoleOrder wanted =
                images ? MultipoleOrder::DipoleQuadrupole : MultipoleOrder::Quadrupole;
        // The solver refers to the container: ensureContainer() drops both on a rebuild.
        if (solver_m != nullptr && solverOrder_m == wanted) {
            return;
        }
        BHSolver::Params params;
        // ryoanji accumulates a_i = G sum_j m_j (r_j - r_i) / r^3. Real particles are packed with
        // m_j = |q_j| and the sign of the common charge is applied in the copy-back; image charges
        // carry -|q_j|. Mixed-sign sources need ryoanji's remote-leaf fallback to recognise any
        // leaf with sources, also one whose net charge is zero.
        params.G          = -1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
        params.theta      = static_cast<float>(config_m.theta);
        params.numShells  = 0;
        params.multipoles = wanted;
        solver_m.reset();
        solver_m      = std::make_unique<BHSolver>(*bh_m, params);
        solverOrder_m = wanted;
    }

    BHLocalIndex BarnesHutAlgorithm::Impl::prepareSlots(size_type nSlots) {
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

    void BarnesHutAlgorithm::Impl::pack(
            BHLocalIndex offset, size_type nLocal, double timeStep, bool images) {
        // Deposition weight of CartesianPIC3D: q_i = Q dt_i / dt_step (emission fractions).
        const double chargeScale = std::abs(primary_m->getChargePerParticle()) / timeStep;
        packParticles(
                *bh_m, *primary_m, offset, nLocal, chargeScale, config_m.softening, images,
                config_m.dirichletPlane.planeZ);
    }

    void BarnesHutAlgorithm::Impl::synchronize() {
        static IpplTimings::TimerRef syncTimer = IpplTimings::getTimer("bh.syncGrav");
        static IpplTimings::TimerRef haloTimer = IpplTimings::getTimer("bh.haloCharge");

        IpplTimings::startTimer(syncTimer);
        ippl::nbody::syncGravBH<Precision, Conserved, Dependent>(*bh_m);
        if (!config_m.leafBasedSoftening) {
            // Halo slots that cstone grew or reused must carry the uniform softening too.
            Kokkos::deep_copy(ippl::nbody::getView<"h">(*bh_m), config_m.softening);
        }
        synchronizeDevice();
        IpplTimings::stopTimer(syncTimer);

        IpplTimings::startTimer(haloTimer);
        bh_m->exchangeHalos(
                std::tie(bh_m->template getDV<"charge">()), bh_m->haloSendBuf(),
                bh_m->haloRecvBuf());
        synchronizeDevice();
        IpplTimings::stopTimer(haloTimer);
    }

    void BarnesHutAlgorithm::Impl::copyBack() {
        // Image charges stay behind; only the real particles return to OPALX.
        const size_type nOwned = countOwnedReal(*bh_m);
        // Collective: every rank resizes, and the global count stays unchanged.
        primary_m->replaceLocalCount(nOwned);
        const double sign = primary_m->getChargePerParticle() < 0.0 ? -1.0 : 1.0;
        copyBackParticles(*bh_m, *primary_m, sign);
        primary_m->markMomentsDirty();
    }

    SpaceChargeSolveResult BarnesHutAlgorithm::Impl::solve(const SpaceChargeSolveContext& context) {
        static IpplTimings::TimerRef packTimer     = IpplTimings::getTimer("bh.pack");
        static IpplTimings::TimerRef solveTimer    = IpplTimings::getTimer("bh.solve");
        static IpplTimings::TimerRef copyBackTimer = IpplTimings::getTimer("bh.copyBack");

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

        enterSolveFrame(context.stepState().frames, *primary_m);

        const std::array<double, 6> bounds = globalBounds(*primary_m);
        const double extent = std::max({bounds[1] - bounds[0], bounds[3] - bounds[2],
                                        bounds[5] - bounds[4]});
        if (!std::isfinite(extent)) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve", "Particle positions are not finite.");
        }
        const auto total = static_cast<unsigned long long>(primary_m->getTotalNum());
        const bool images = imagesActive(context.stepState().step);
        if (config_m.dirichletPlane.enabled() && images != imagesWereActive_m) {
            Inform m("BarnesHutAlgorithm");
            m << level2 << "TYPE=BH image charges of the plane z = "
              << config_m.dirichletPlane.planeZ << " m are " << (images ? "active" : "expired")
              << " at step " << context.stepState().step << "." << endl;
            imagesWereActive_m = images;
        }
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

        if (useDirectSum) {
            static IpplTimings::TimerRef directTimer = IpplTimings::getTimer("bh.directSum");
            IpplTimings::startTimer(directTimer);
            // Leaf-based softening has no tree here; a negligible floor only defuses the self pair.
            const double softening =
                    config_m.softening > 0.0 ? config_m.softening : 1.0e-9 * extent;
            const double charge = primary_m->getChargePerParticle();
            const double prefactor =
                    (charge < 0.0 ? 1.0 : -1.0) / (4.0 * Physics::pi * Physics::epsilon_0);
            directSumFields(
                    *primary_m, std::abs(charge) / context.stepState().timeStep, softening,
                    prefactor, images, config_m.dirichletPlane.planeZ);
            IpplTimings::stopTimer(directTimer);
            ++statistics_m.directSolves;
            result.backendSolves = 1;
        } else if (extent > 0.0) {
            IpplTimings::startTimer(packTimer);
            std::array<double, 6> sourceBounds = bounds;
            if (images) {
                // The mirrored bunch extends the source region beyond the plane.
                const double planeZ = config_m.dirichletPlane.planeZ;
                sourceBounds[4]     = std::min(bounds[4], 2.0 * planeZ - bounds[5]);
                sourceBounds[5]     = std::max(bounds[5], 2.0 * planeZ - bounds[4]);
            }
            ensureContainer(sourceBounds, images ? 2 * total : total);
            ensureSolver(images);
            const size_type nLocal    = primary_m->getLocalNum();
            const BHLocalIndex offset = prepareSlots(images ? 2 * nLocal : nLocal);
            pack(offset, nLocal, context.stepState().timeStep, images);
            IpplTimings::stopTimer(packTimer);

            synchronize();

            IpplTimings::startTimer(solveTimer);
            try {
                solver_m->runSolver(/*warmup=*/true);
            } catch (const std::runtime_error& error) {
                throw OpalException(
                        "BarnesHutAlgorithm::solve",
                        std::string(error.what())
                                + " Increase FIELDSOLVER BHTHETA to shorten the traversal.");
            }
            synchronizeDevice();
            IpplTimings::stopTimer(solveTimer);

            IpplTimings::startTimer(copyBackTimer);
            copyBack();
            IpplTimings::stopTimer(copyBackTimer);

            ++statistics_m.solves;
            if (solverOrder_m == ippl::nbody::MultipoleOrder::DipoleQuadrupole) {
                ++statistics_m.dipoleSolves;
            }
            result.backendSolves = 1;
        }
        // All particles coincide otherwise: the softened field vanishes and E stays zero.

        if (result.backendSolves > 0) {
            verifySolve(context, totalBefore, images);
        }

        leaveSolveFrame(context.stepState().frames, *primary_m);
        return result;
    }

    void BarnesHutAlgorithm::Impl::verifySolve(
            const SpaceChargeSolveContext& context, size_type totalBefore, bool images) {
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
        // |E| for the exact sum; multipole errors are orders of magnitude smaller. Leaf-based
        // softening has no single h, so only the uniform case is bounded.
        if (!config_m.leafBasedSoftening && config_m.softening > 0.0) {
            // Image charges double the source charge.
            const double totalCharge = (images ? 2.0 : 1.0)
                                       * std::abs(primary_m->getChargePerParticle()) * stats.sumDt
                                       / context.stepState().timeStep;
            const double coulomb     = 1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
            const double bound =
                    2.0 * coulomb * totalCharge / (4.0 * config_m.softening * config_m.softening);
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
        impl_m = std::make_unique<Impl>(std::move(config), primary, std::move(bunchState));
    }

    BarnesHutAlgorithm::~BarnesHutAlgorithm() = default;

    SpaceChargeSolveResult BarnesHutAlgorithm::solve(const SpaceChargeSolveContext& context) {
        return impl_m->solve(context);
    }

    const BarnesHutAlgorithm::Statistics& BarnesHutAlgorithm::statistics() const {
        return impl_m->statistics_m;
    }

}  // namespace opalx::spacecharge
