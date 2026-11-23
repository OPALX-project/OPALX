#include "BeamlineCore/ConstantEFieldCavityRep.h"

ConstantEFieldCavityRep::ConstantEFieldCavityRep() : ConstantEFieldCavity(), geometry() {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const ConstantEFieldCavityRep& right)
    : ConstantEFieldCavity(right), geometry(right.geometry) {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const std::string& name)
    : ConstantEFieldCavity(name), geometry() {}

ConstantEFieldCavityRep::~ConstantEFieldCavityRep() {}

ElementBase* ConstantEFieldCavityRep::clone() const { return new ConstantEFieldCavityRep(*this); }

Geometry& ConstantEFieldCavityRep::getGeometry() { return geometry; }

const Geometry& ConstantEFieldCavityRep::getGeometry() const { return geometry; }
