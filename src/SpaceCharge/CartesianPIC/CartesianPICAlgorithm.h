/**
 * @file CartesianPICAlgorithm.h
 * @brief Declares the complete Cartesian 3D PIC space-charge algorithm.
 */

#ifndef OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H
#define OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H

#include "SpaceCharge/CartesianPIC/CartesianDomainUpdater.h"
#include "SpaceCharge/CartesianPIC/CartesianPICFieldStorage.h"
#include "SpaceCharge/CartesianPIC/CorrectionPassSchedule.h"
#include "SpaceCharge/CartesianPIC/P3MShortRangeInteraction.h"
#include "SpaceCharge/CartesianPIC/ParticleBinTraversal.h"
#include "SpaceCharge/CartesianPIC/ParticleDomainOperations.h"
#include "SpaceCharge/CartesianPIC/ParticleMeshFieldTransfer.h"
#include "SpaceCharge/CartesianPIC/RelativisticFieldComposer.h"
#include "SpaceCharge/Poisson/PoissonSolver.h"
#include "SpaceCharge/SpaceChargeAlgorithm.h"
#include "SpaceCharge/SpaceChargeConfig.h"
#include "SpaceCharge/SpaceChargeFrameTransform.h"

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
     * always the tracker frame; SpaceChargeFrameTransform owns the 3D beam-frame conversion and
     * restoration contract.
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
        using ParticleBinTraversalType      = ParticleBinTraversal<double, 3>;
        using ParticleBinType               = ParticleBin<double, 3>;
        using PassProperties                = CartesianPICPassProperties<double, 3>;
        using ParticleMeshTransfer          = ParticleMeshFieldTransfer<double, 3>;
        using RelativisticFieldComposerType = RelativisticFieldComposer<double, 3>;

        CartesianPICAlgorithm(
                CartesianPICConfig config, std::span<const ParticleFieldBinding3D> particleBindings,
                std::unique_ptr<FieldStorage> fieldStorage, DataSink* dataSink,
                std::shared_ptr<const BunchStateHandler> bunchState);

        [[nodiscard]] SpaceChargeCapabilities capabilities() const override;
        void solve(SpaceChargeSolveContext& context, SpaceChargeDiagnostics& diagnostics) override;

    private:
        struct BinStatsRow final {
            long long binNumber;
            unsigned long long particleCount;
            double gamma;
        };

        void solveInBeamFrame(
                SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
                SpaceChargeDiagnostics& diagnostics);
        void solveWholeBunch(
                SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
                SpaceChargeDiagnostics& diagnostics);
        void solveBinned(
                SpaceChargeSolveContext& context, const CorrectionPassSequence& correction,
                SpaceChargeDiagnostics& diagnostics);
        void solvePass(
                SpaceChargeSolveContext& context, const ParticleBinType* unit,
                CartesianPICPass pass, bool& dumpedDirichletPlaneThisStep,
                SpaceChargeDiagnostics& diagnostics);
        void depositChargeForBin(
                const SpaceChargeSolveContext& context, const ParticleBinType& unit,
                const PassProperties& pass);
        void dumpBinSnapshot(
                const SpaceChargeSolveContext& context, const BinConfigurationSnapshot& snapshot,
                bool beforeMerge) const;
        void dumpDirichletPlaneDiagnosticsIfRequested(
                const SpaceChargeSolveContext& context, const std::string& solveTag, double planeZ,
                SpaceChargeDiagnostics& diagnostics);
        void printBinStatsTable() const;

        ParticleContainer* primary_m = nullptr;
        std::shared_ptr<const BunchStateHandler> bunchState_m;
        std::unique_ptr<FieldStorage> fieldStorage_m;
        std::unique_ptr<PoissonSolver> poissonSolver_m;
        std::optional<P3MShortRangeInteraction> shortRangeInteraction_m;
        PoissonSolverType poissonSolverType_m;
        DataSink* dataSink_m = nullptr;
        std::optional<BinningConfig> binningConfig_m;
        std::unique_ptr<ParticleBinTraversalType> particleBinTraversal_m;
        CorrectionPassSchedule correctionPassSchedule_m;
        ParticleDomainOperations particleDomainOperations_m;
        CartesianDomainUpdater domainUpdater_m;
        ParticleMeshTransfer particleMeshTransfer_m;
        RelativisticFieldComposerType relativisticFieldComposer_m;
        std::vector<BinStatsRow> binStats_m;
        bool warnedPlaneDumpParallelUnsupported_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_CARTESIAN_PIC_ALGORITHM_H
