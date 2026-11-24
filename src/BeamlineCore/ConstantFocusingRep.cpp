#include "BeamlineCore/ConstantFocusingRep.h"

ConstantFocusingRep::ConstantFocusingRep() : ConstantFocusing() {}

ConstantFocusingRep::ConstantFocusingRep(const ConstantFocusingRep& right)
    : ConstantFocusing(right) {}

ConstantFocusingRep::ConstantFocusingRep(const std::string& name) : ConstantFocusing(name) {}

ConstantFocusingRep::~ConstantFocusingRep() {}

ElementBase* ConstantFocusingRep::clone() const { return new ConstantFocusingRep(*this); }
