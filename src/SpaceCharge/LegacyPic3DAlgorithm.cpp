/**
 * @file LegacyPic3DAlgorithm.cpp
 * @brief Implements the temporary PartBunch self-field bridge.
 */

#include "SpaceCharge/LegacyPic3DAlgorithm.h"

#include "PartBunch/PartBunch.h"

#include <utility>

namespace opalx::spacecharge {

    LegacyPic3DAlgorithm::LegacyPic3DAlgorithm(
            PartBunch_t& bunch, std::shared_ptr<PicWorkspace<double, 3>> workspace)
        : workspace_m(std::move(workspace)), bunch_m(&bunch) {
        if (workspace_m == nullptr) {
            throw OpalException(
                    "LegacyPic3DAlgorithm::LegacyPic3DAlgorithm", "PIC workspace is null.");
        }
    }

    SolverCapabilities LegacyPic3DAlgorithm::capabilities() const {
        SolverCapabilities result;
        result.algorithm                  = SelfFieldAlgorithmKind::Pic3D;
        result.supportsBinning            = true;
        result.supportsImageCharge        = true;
        result.supportsShiftedGreen       = true;
        result.supportsRedistribution     = true;
        result.supportsMultipleContainers = false;
        result.supportsPotentialOutput    = true;

        result.requiredReadableAttributes =
                ParticleAttribute::Position | ParticleAttribute::Momentum
                | ParticleAttribute::Charge | ParticleAttribute::Mass | ParticleAttribute::TimeStep
                | ParticleAttribute::InvalidMask | ParticleAttribute::Bin;
        result.requiredWritableAttributes =
                ParticleAttribute::Position | ParticleAttribute::TimeStep
                | ParticleAttribute::ElectricField | ParticleAttribute::MagneticField
                | ParticleAttribute::Bin;
        return result;
    }

    void LegacyPic3DAlgorithm::execute(
            SolveContext& /*context*/, SelfFieldDiagnostics& diagnostics) {
        bunch_m->computeSelfFields(diagnostics);
    }

}  // namespace opalx::spacecharge
