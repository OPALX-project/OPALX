#include "SelfFieldClearTestOps.h"

namespace opalx::spacecharge::testing {

    void assignmentClearOtherTranslationUnit(SpaceChargeParticleContainer& particles) {
        assignmentClear(particles);
    }

    void fencedAssignmentClearOtherTranslationUnit(SpaceChargeParticleContainer& particles) {
        fencedAssignmentClear(particles);
    }

    void deepCopyClearOtherTranslationUnit(SpaceChargeParticleContainer& particles) {
        clearSelfFields(particles);
    }

}  // namespace opalx::spacecharge::testing
