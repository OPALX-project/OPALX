#include "BeamlineCore/RBendRep.h"

RBendRep::RBendRep() : RBend() {}

RBendRep::RBendRep(const RBendRep& right) : RBend(right) {}

RBendRep::RBendRep(const std::string& name) : RBend(name) {}

RBendRep::~RBendRep() = default;
