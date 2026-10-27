/**
 * @file BarnesHutAlgorithm.h
 * @brief Gridless Barnes-Hut space-charge algorithm (FIELDSOLVER TYPE=BH).
 */

#ifndef OPALX_SPACE_CHARGE_BARNES_HUT_ALGORITHM_H
#define OPALX_SPACE_CHARGE_BARNES_HUT_ALGORITHM_H

#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"

#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <string>

template <typename T, unsigned Dim>
class ParticleContainer;
class BunchStateHandler;

namespace opalx::spacecharge {

    /**
     * @brief Computes the electrostatic self-field of the primary container with IPPL's
     * Barnes-Hut module (cstone octree + ryoanji multipoles).
     *
     * Semantics match unbinned CartesianPIC3D with OPEN boundaries: the source charge of particle
     * i is Q dt_i / dt_step, the field is E = 1/(4 pi eps0) sum_j q_j (r_i - r_j) / max(|r_i -
     * r_j|^2, (h_i + h_j)^2)^{3/2} in V/m in the solve frame, gamma = 1 and B = 0. Other
     * containers keep their E/B.
     *
     * With BarnesHutKinematics::Frame (BINS with MAXBINS=1) it matches CartesianPIC3D with one
     * bin instead: the field is solved in the rest frame of the bunch (z stretched by the mean
     * gamma) and boosted to E and B, image charges with the opposite B sign. Image charges far
     * from the bunch in the rest frame are evaluated as one quadrupole.
     *
     * With BarnesHutKinematics::Particle (BHRELATIVISTIC=PARTICLE, or BINS with MAXBINS > 1)
     * every source contributes the field of a charge moving uniformly along the solve-frame z
     * with its own velocity, E_j = k q_j gamma_j R / (rho^2 + gamma_j^2 zeta^2)^{3/2} and
     * B_j = v_j x E_j / c^2, softened in the source's rest frame; the limit of energy binning for
     * infinitely many bins. Image charges carry -q and the reversed velocity. The octree cells
     * carry rapidity moments (ryoanji HeavisideQpole); far images are one such multipole.
     *
     * Barnes-Hut owns the primary's decomposition: every solve hands all primary particles to
     * cstone, which redistributes them along its space-filling curve, and copies the owned set
     * back (R, P, dt, ID and the new E). Particle order and per-rank counts change; the global
     * particle set is preserved. The IPPL spatial layout of the primary is stale afterwards and
     * ParticleContainer::update() is disabled.
     *
     * NBody/cstone headers are confined to the implementation file.
     */
    class BarnesHutAlgorithm final : public SpaceChargeAlgorithm {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;

        /** @brief Counters for diagnostics and tests. */
        struct Statistics {
            std::size_t solves          = 0;  ///< Completed Barnes-Hut solves.
            std::size_t inPlacePacks    = 0;  ///< Solves that reused the owned slots.
            std::size_t containerBuilds = 0;  ///< Constructions of the Barnes-Hut container.
            std::size_t directSolves    = 0;  ///< Small-bunch solves done by direct summation.
            std::size_t dipoleSolves    = 0;  ///< Barnes-Hut solves with dipole multipoles (images).
            /// Container rebuilds for a resized box (also counted in containerBuilds): the
            /// sources outgrew the box (BHBOXMARGIN) or the box became too large (BHBOXSLACK).
            std::size_t boxRebuilds = 0;
            std::size_t boxShrinkRebuilds = 0;  ///< Of boxRebuilds: the box was too large.
            /// Container rebuilds for automatic bucket sizes far off the particle count (emission,
            /// losses; also counted in containerBuilds).
            std::size_t bucketRebuilds = 0;
            /// Solves whose box was shifted, at unchanged size, to follow the sources.
            std::size_t boxMoves = 0;
            /// SFC key bits per axis of the last Barnes-Hut domain (fewer on short axes).
            std::array<unsigned, 3> keyBits{};
            std::size_t restFrameSolves = 0;  ///< Solves boosted from the bunch rest frame.
            /// Second traversals (or direct sums) for the field of the image charges alone.
            std::size_t imagePasses = 0;
            std::size_t farImageSolves = 0;  ///< Solves with the images as one quadrupole.
            /// Solves (tree or direct sum) with per-particle source velocities (PARTICLE).
            std::size_t particleSolves = 0;
            /// PARTICLE on several ranks: remote leaves applied as multipoles although they failed
            /// the rapidity-scaled MAC (cstone's locally essential tree is not velocity-aware).
            unsigned long long remoteMacFailures = 0;
            double maxGamma = 1.0;  ///< Largest source gamma of the last per-particle solve.
            double gamma = 1.0;  ///< Mean gamma of the last rest-frame solve.
            double softening = 0.0;  ///< Uniform softening h of the last solve (BHSOFTENINGREL).
        };

        BarnesHutAlgorithm(
                BarnesHutConfig config, std::span<ParticleContainer* const> particles,
                std::shared_ptr<const BunchStateHandler> bunchState);
        ~BarnesHutAlgorithm() override;

        BarnesHutAlgorithm(const BarnesHutAlgorithm&)            = delete;
        BarnesHutAlgorithm& operator=(const BarnesHutAlgorithm&) = delete;

        [[nodiscard]] SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override;

        [[nodiscard]] const Statistics& statistics() const;

        /** @brief Solve paths and domain events of the run (box moves, rebuilds by cause). */
        [[nodiscard]] std::string runSummary() const override;

    private:
        class Impl;
        /// Impl for one IPPL NBody precision policy (BarnesHutConfig::precision).
        template <class Precision>
        class PrecisionImpl;
        std::unique_ptr<Impl> impl_m;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_BARNES_HUT_ALGORITHM_H
