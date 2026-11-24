#include "BeamlineCore/ConstantEFieldCavityRep.h"

ConstantEFieldCavityRep::ConstantEFieldCavityRep() : ConstantEFieldCavity() {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const ConstantEFieldCavityRep& right)
    : ConstantEFieldCavity(right) {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const std::string& name)
    : ConstantEFieldCavity(name) {}

ConstantEFieldCavityRep::~ConstantEFieldCavityRep() {}
