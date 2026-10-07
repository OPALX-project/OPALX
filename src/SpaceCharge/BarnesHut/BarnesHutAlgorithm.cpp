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

#include "Physics/Physics.h"
#include "Utility/IpplTimings.h"

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>

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

        /** @brief Write the local OPALX particles into BH slots [offset, offset + nLocal). */
        void packParticles(
                BHContainer& bh, BarnesHutAlgorithm::ParticleContainer& particles,
                BHLocalIndex offset, size_type nLocal, double chargeScale, double softening) {
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
                    });
            Kokkos::fence();
        }

        /** @brief Overwrite the (already resized) OPALX particles with the BH-owned set. */
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

            Kokkos::parallel_for(
                    "BarnesHutAlgorithm::copyBack", n, KOKKOS_LAMBDA(const size_type k) {
                        const size_type j = first + k;
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
        unsigned bucketSize_m = 0;
        Statistics statistics_m;

    private:
        void checkPreconditions(const SpaceChargeSolveContext& context) const;
        void ensureContainer(const std::array<double, 6>& bounds);
        BHLocalIndex prepareSlots(size_type nLocal);
        void pack(BHLocalIndex offset, size_type nLocal, double timeStep);
        void synchronize();
        void copyBack();
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

    void BarnesHutAlgorithm::Impl::ensureContainer(const std::array<double, 6>& bounds) {
        const int nRanks = ippl::Comm->size();
        const auto total = static_cast<unsigned long long>(primary_m->getTotalNum());
        unsigned bucketSize =
                config_m.bucketSize != 0
                        ? static_cast<unsigned>(config_m.bucketSize)
                        : std::max<unsigned>(
                                kMinimumBucketSize,
                                static_cast<unsigned>(
                                        total / (100ull * static_cast<unsigned long long>(nRanks))));

        // The global-tree bucket size is fixed per cstone domain. Emission can grow the bunch by
        // orders of magnitude; rebuild once the automatic size is far off. All particles live in
        // the OPALX container between solves, so a rebuild loses no state.
        const bool rebuild =
                bh_m != nullptr && config_m.bucketSize == 0 && bucketSize > 4 * bucketSize_m;
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
                ippl::Comm->rank(), nRanks, bucketSize,
                static_cast<unsigned>(config_m.bucketSizeFocus), static_cast<float>(config_m.theta),
                box,
                std::array<BoundaryType, 3>{
                        BoundaryType::cubic_open, BoundaryType::cubic_open,
                        BoundaryType::cubic_open},
                ippl::Comm->getCommunicator());
        bh_m->setUniformH(config_m.softening);
        bh_m->setLeafBasedH(config_m.leafBasedSoftening);

        BHSolver::Params params;
        // ryoanji accumulates a_i = G sum_j m_j (r_j - r_i) / r^3. Packing m_j = |q_j| keeps every
        // multipole mass positive (ryoanji's remote-leaf fallback skips cells with mass <= 0),
        // and the sign of the common charge is applied in the copy-back.
        params.G         = -1.0 / (4.0 * Physics::pi * Physics::epsilon_0);
        params.theta     = static_cast<float>(config_m.theta);
        params.numShells = 0;
        solver_m         = std::make_unique<BHSolver>(*bh_m, params);
        bucketSize_m     = bucketSize;
        ++statistics_m.containerBuilds;
    }

    BHLocalIndex BarnesHutAlgorithm::Impl::prepareSlots(size_type nLocal) {
        if (nLocal > static_cast<size_type>(std::numeric_limits<BHLocalIndex>::max() / 2)) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve",
                    "Too many local particles for the 32-bit Barnes-Hut indices.");
        }
        const auto n = static_cast<BHLocalIndex>(nLocal);

        // After a copy-back the owned slots already belong to this rank's particle set; when the
        // count is unchanged they are overwritten in place and the next sync only moves particles
        // that crossed a space-filling-curve boundary.
        if (statistics_m.solves > 0 && bh_m->getLocalNum() == n) {
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

    void BarnesHutAlgorithm::Impl::pack(BHLocalIndex offset, size_type nLocal, double timeStep) {
        // Deposition weight of CartesianPIC3D: q_i = Q dt_i / dt_step (emission fractions).
        const double chargeScale = std::abs(primary_m->getChargePerParticle()) / timeStep;
        packParticles(*bh_m, *primary_m, offset, nLocal, chargeScale, config_m.softening);
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
        const size_type nOwned = bh_m->getLocalNum();
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

        enterSolveFrame(context.stepState().frames, *primary_m);

        const std::array<double, 6> bounds = globalBounds(*primary_m);
        const double extent = std::max({bounds[1] - bounds[0], bounds[3] - bounds[2],
                                        bounds[5] - bounds[4]});
        if (!std::isfinite(extent)) {
            throw OpalException(
                    "BarnesHutAlgorithm::solve", "Particle positions are not finite.");
        }
        if (extent > 0.0) {
            IpplTimings::startTimer(packTimer);
            ensureContainer(bounds);
            const size_type nLocal = primary_m->getLocalNum();
            const BHLocalIndex offset = prepareSlots(nLocal);
            pack(offset, nLocal, context.stepState().timeStep);
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
            result.backendSolves = 1;
        }
        // All particles coincide otherwise: the softened field vanishes and E stays zero.

        leaveSolveFrame(context.stepState().frames, *primary_m);
        return result;
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
