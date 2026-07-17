/**
 * @file LegacyPic3DAlgorithm.h
 * @brief Temporary adapter from SelfFieldAlgorithm to the current PartBunch solver path.
 */

#ifndef OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H
#define OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H

#include "PartBunch/PartBunchFwd.h"
#include "SpaceCharge/SelfFieldAlgorithm.h"

namespace opalx::spacecharge {

    /**
     * @brief Transitional 3D PIC implementation that borrows run-owned legacy state.
     *
     * No device kernel sees this adapter or its PartBunch pointer. Later stages replace this
     * class with a solver that owns its domain, workspace, and backend directly.
     */
    class LegacyPic3DAlgorithm final : public SelfFieldAlgorithm {
    public:
        explicit LegacyPic3DAlgorithm(PartBunch_t& bunch);

        [[nodiscard]] SolverCapabilities capabilities() const override;
        void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) override;

    private:
        PartBunch_t* bunch_m = nullptr;
    };

}  // namespace opalx::spacecharge

#endif  // OPALX_SPACE_CHARGE_LEGACY_PIC3D_ALGORITHM_H
