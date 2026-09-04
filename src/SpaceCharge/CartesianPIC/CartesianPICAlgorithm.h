/**
 * @file CartesianPICAlgorithm.h
 * @brief Declares the complete Cartesian 3D PIC space-charge algorithm.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H

#include "SpaceCharge/CartesianPIC/CartesianDomainUpdater.h"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/P3MShortRangeInteraction.h"
#include "SpaceCharge/CartesianPIC/ParticleBinTraversal.h"
#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"
#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeFrames.h"

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

class DataSink;
class BunchStateHandler;

namespace opalx::spacecharge {

    /**
     * @brief Owns all runtime orchestration for the existing Cartesian 3D PIC algorithm.
     *
     * The algorithm borrows stable particle containers and a data sink while owning its field
     * storage, Poisson solver, and all 3D orchestration components. It never retains PartBunch,
     * parser objects, per-call transforms, or native Kokkos views. The common solve() boundary is
     * always the tracker frame and explicit frame helpers restore R/E/B after successful solves.
     *
     * Each solve transforms primary positions to the beam frame, updates geometry/layouts, deposits
     * charge, invokes the configured backend passes, gathers and composes fields, then restores
     * primary R/E/B to the tracker frame. Normal mode migrates every container and next rebuilds a
     * reference-frame domain. Fixed-domain mode migrates only the primary, keeps the beam-frame
     * mesh and decomposition for BeamBeam reuse, and recomputes its restored-coordinate moments
     * without a second migration.
     */
    class CartesianPICAlgorithm final : public SpaceChargeAlgorithm {
    public:
        using ParticleContainer             = ::ParticleContainer<double, 3>;
        using FieldStorage                  = CartesianPICFieldStorage<double, 3>;
        using ParticleBinTraversalType      = ParticleBinTraversal;
        using ParticleBinType               = ParticleBin;
        using ParticleMeshTransfer          = ParticleMeshFieldTransfer;
        using RelativisticFieldComposerType = RelativisticFieldComposer;

        CartesianPICAlgorithm(
                CartesianPICConfig config, std::span<ParticleContainer* const> particles,
                std::unique_ptr<FieldStorage> fieldStorage, DataSink* dataSink,
                std::shared_ptr<const BunchStateHandler> bunchState);

        [[nodiscard]] SpaceChargeSolveResult solve(const SpaceChargeSolveContext& context) override;

    private:
        enum class PassKind { Primary, PrimaryAndImage, Image, ShiftedImage };

        struct SolvePlan {
            std::array<PassKind, 2> passes{};
            std::size_t passCount = 0;
            CorrectionConfig activeCorrection;
            bool correctionExpired = false;
        };

        struct PassProperties {
            ParticleMeshTransfer::DepositKind depositKind =
                    ParticleMeshTransfer::DepositKind::Primary;
            ParticleMeshTransfer::ImagePolicy imagePolicy;
            bool shiftedGreen            = false;
            bool suppressFieldDump       = false;
            double magneticSign          = 1.0;
            FieldSourceRule sourceRule   = FieldSourceRule::Direct;
            bool dumpDirichletPlaneAfter = false;
            const char* label            = "primary";
        };

        struct BinStatsRow final {
            long long binNumber;
            unsigned long long particleCount;
            double gamma;
        };

        [[nodiscard]] SolvePlan makeSolvePlan(std::size_t step) const;
        [[nodiscard]] static PassProperties passProperties(
                PassKind pass, double planeZ, bool binned);
        void solveInBeamFrame(
                const SpaceChargeSolveContext& context, const SolvePlan& plan,
                SpaceChargeSolveResult& result);
        void solveWholeBunch(
                const SpaceChargeSolveContext& context, const SolvePlan& plan,
                SpaceChargeSolveResult& result);
        void solveBinned(
                const SpaceChargeSolveContext& context, const SolvePlan& plan,
                SpaceChargeSolveResult& result);
        void solvePass(
                const SpaceChargeSolveContext& context, const SolvePlan& plan,
                const ParticleBinType* unit, PassKind pass, SpaceChargeSolveResult& result);
        void depositChargeForBin(
                const SpaceChargeSolveContext& context, const ParticleBinType& unit,
                const PassProperties& pass);
        void dumpBinSnapshot(
                const SpaceChargeSolveContext& context, const BinConfigurationSnapshot& snapshot,
                bool beforeMerge) const;
        void dumpDirichletPlaneDiagnosticsIfRequested(
                const SpaceChargeSolveContext& context, const std::string& solveTag, double planeZ);
        void printBinStatsTable() const;

        CartesianPICConfig config_m;
        ParticleContainer* primary_m = nullptr;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<FieldStorage> fieldStorage_m;
        std::unique_ptr<PoissonSolver> poissonSolver_m;
        std::optional<P3MShortRangeInteraction> shortRangeInteraction_m;
        PoissonSolverType poissonSolverType_m;
        DataSink* dataSink_m = nullptr;
        std::optional<BinningConfig> binningConfig_m;
        std::unique_ptr<ParticleBinTraversalType> particleBinTraversal_m;
        CartesianDomainUpdater domainUpdater_m;
        ParticleMeshTransfer particleMeshTransfer_m;
        RelativisticFieldComposerType relativisticFieldComposer_m;
        std::vector<BinStatsRow> binStats_m;
        bool warnedPlaneDumpParallelUnsupported_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H
