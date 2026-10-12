#ifndef OPALX_SELF_FIELD_CLEAR_TEST_OPS_H
#define OPALX_SELF_FIELD_CLEAR_TEST_OPS_H

#include "SpaceCharge/SpaceChargeFrames.h"

namespace opalx::spacecharge::testing {

    using ClearOperation = void (*)(SpaceChargeParticleContainer&);

    // Deliberately instantiate the original IPPL assignment from two translation units.
    // CUDA extended-lambda launch-wrapper failures can depend on the linked binary.
    inline void assignmentClear(SpaceChargeParticleContainer& particles) {
        particles.E = Vector_t<double, 3>(0.0);
        particles.B = Vector_t<double, 3>(0.0);
    }

    inline void fencedAssignmentClear(SpaceChargeParticleContainer& particles) {
        Kokkos::fence("SelfFieldClearTest::beforeE");
        particles.E = Vector_t<double, 3>(0.0);
        Kokkos::fence("SelfFieldClearTest::afterE");
        Kokkos::fence("SelfFieldClearTest::beforeB");
        particles.B = Vector_t<double, 3>(0.0);
        Kokkos::fence("SelfFieldClearTest::afterB");
    }

    void assignmentClearOtherTranslationUnit(SpaceChargeParticleContainer& particles);
    void fencedAssignmentClearOtherTranslationUnit(SpaceChargeParticleContainer& particles);
    void deepCopyClearOtherTranslationUnit(SpaceChargeParticleContainer& particles);

}  // namespace opalx::spacecharge::testing

#endif
