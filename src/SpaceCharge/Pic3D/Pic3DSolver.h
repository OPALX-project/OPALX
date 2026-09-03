/**
 * @file Pic3DSolver.h
 * @brief Declares the complete Cartesian 3D PIC self-field algorithm.
 */

#ifndef OPALX_SPACE_CHARGE_PIC3D_PIC3D_SOLVER_H
#define OPALX_SPACE_CHARGE_PIC3D_PIC3D_SOLVER_H

#include "SpaceCharge/Ippl/IpplPoissonAdapter.h"
#include "SpaceCharge/ParticleFrameGuard.h"
#include "SpaceCharge/Pic3D/BinningPlan.h"
#include "SpaceCharge/Pic3D/CorrectionPlan.h"
#include "SpaceCharge/Pic3D/FieldComposer.h"
#include "SpaceCharge/Pic3D/P3MShortRangeInteraction.h"
#include "SpaceCharge/Pic3D/PicDomainManager.h"
#include "SpaceCharge/Pic3D/PicParticleDomainAdapter.h"
#include "SpaceCharge/Pic3D/PicScatterGather.h"
#include "SpaceCharge/Pic3D/PicWorkspace.h"
#include "SpaceCharge/SelfFieldAlgorithm.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

class DataSink;

namespace opalx::spacecharge {

    /**
     * @brief Owns all runtime orchestration for the existing Cartesian 3D PIC algorithm.
     *
     * The solver borrows stable particle containers and a data sink while owning its workspace,
     * Poisson backend, and all 3D orchestration components. It never retains PartBunch, parser
     * objects, per-call transforms, or native Kokkos views. The common execute() boundary is
     * always the tracker frame; ParticleFrameGuard owns the 3D beam-frame conversion and
     * restoration contract.
     */
    class Pic3DSolver final : public SelfFieldAlgorithm {
    public:
        using ParticleContainer = ::ParticleContainer<double, 3>;
        using Workspace         = PicWorkspace<double, 3>;
        using BinningPlan_t     = BinningPlan<double, 3>;
        using BinnedSolveUnit_t = BinnedSolveUnit<double, 3>;
        using PassPolicy        = SolvePassPolicy<double, 3>;
        using ScatterGather     = PicScatterGather<double, 3>;
        using FieldComposer_t   = FieldComposer<double, 3>;

        Pic3DSolver(
                Pic3DConfig config, std::span<const ParticleFieldBinding3d> particleBindings,
                std::unique_ptr<Workspace> workspace, DataSink* dataSink);

        [[nodiscard]] SolverCapabilities capabilities() const override;
        void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) override;

    private:
        struct BinStatsRow final {
            long long binNumber;
            unsigned long long particleCount;
            double gamma;
        };

        void solveInBeamFrame(
                SolveContext& context, const PreparedCorrection& correction,
                SelfFieldDiagnostics& diagnostics);
        void executeWholeBunch(
                SolveContext& context, const PreparedCorrection& correction,
                SelfFieldDiagnostics& diagnostics);
        void executeBinned(
                SolveContext& context, const PreparedCorrection& correction,
                SelfFieldDiagnostics& diagnostics);
        void executeSolvePass(
                SolveContext& context, const BinnedSolveUnit_t* unit, SolvePass pass,
                bool& dumpedDirichletPlaneThisStep, SelfFieldDiagnostics& diagnostics);
        void prepareRhoForUnit(
                const SolveContext& context, const BinnedSolveUnit_t& unit, const PassPolicy& pass);
        void dumpBinSnapshot(
                const SolveContext& context, const BinConfigurationSnapshot& snapshot,
                bool beforeMerge) const;
        void dumpDirichletPlaneDiagnosticsIfRequested(
                const SolveContext& context, const std::string& solveTag, double planeZ,
                SelfFieldDiagnostics& diagnostics);
        void printBinStatsTable() const;

        ParticleContainer* primary_m = nullptr;
        std::unique_ptr<Workspace> workspace_m;
        std::unique_ptr<IpplPoissonAdapter> backend_m;
        std::optional<P3MShortRangeInteraction> shortRangeInteraction_m;
        PoissonBackendKind backendKind_m;
        DataSink* dataSink_m = nullptr;
        std::optional<BinningConfig> binningConfig_m;
        std::unique_ptr<BinningPlan_t> binningPlan_m;
        CorrectionPlan correctionPlan_m;
        PicParticleDomainAdapter particleDomain_m;
        PicDomainManager domainManager_m;
        ScatterGather scatterGather_m;
        FieldComposer_t fieldComposer_m;
        std::vector<BinStatsRow> binStats_m;
        bool warnedPlaneDumpParallelUnsupported_m = false;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_PIC3D_PIC3D_SOLVER_H
