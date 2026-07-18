/**
 * @file LegacyPic3DAlgorithm.h
 * @brief Temporary adapter from SelfFieldAlgorithm to the current PartBunch solver path.
 */

#ifndef OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H
#define OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/Pic/CorrectionPlan.h"
#include "SpaceCharge/Pic/IterationPlan.h"
#include "SpaceCharge/Pic/PicDomainManager.h"
#include "SpaceCharge/Pic/PicParticleDomainAdapter.h"
#include "SpaceCharge/Pic/PicWorkspace.h"
#include "SpaceCharge/SelfFieldAlgorithm.h"
#include "SpaceCharge/SelfFieldConfig.h"

#include <memory>
#include <optional>

namespace opalx::spacecharge {

    /**
     * @brief Transitional 3D PIC implementation that borrows run-owned legacy state.
     *
     * No device kernel sees this adapter or its PartBunch pointer. Later stages replace this
     * class with a solver that owns its domain, workspace, and backend directly.
     */
    class LegacyPic3DAlgorithm final : public SelfFieldAlgorithm {
    public:
        using IterationPlan_t  = IterationPlan<double, 3>;
        using CorrectionPlan_t = CorrectionPlan<double, 3>;

        LegacyPic3DAlgorithm(
                PartBunch_t& bunch, Pic3DConfig config,
                std::shared_ptr<PicWorkspace<double, 3>> workspace);

        [[nodiscard]] SolverCapabilities capabilities() const override;
        void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) override;

    private:
        std::shared_ptr<PicWorkspace<double, 3>> workspace_m;
        std::optional<BinningConfig> binningConfig_m;
        std::unique_ptr<IterationPlan_t> iterationPlan_m;
        CorrectionPlan_t correctionPlan_m;
        PicParticleDomainAdapter particleDomain_m;
        PicDomainManager domainManager_m;
        PartBunch_t* bunch_m = nullptr;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H
