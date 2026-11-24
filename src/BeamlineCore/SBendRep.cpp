#include "BeamlineCore/SBendRep.h"

SBendRep::SBendRep() : SBend() {}

SBendRep::SBendRep(const SBendRep& right) : SBend(right) {}

SBendRep::SBendRep(const std::string& name) : SBend(name) {}

SBendRep::~SBendRep() = default;
