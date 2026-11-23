#include "BeamlineCore/ConstantFocusingRep.h"

ConstantFocusingRep::ConstantFocusingRep() : ConstantFocusing(), geometry() {}

ConstantFocusingRep::ConstantFocusingRep(const ConstantFocusingRep& right)
    : ConstantFocusing(right), geometry(right.geometry) {}

ConstantFocusingRep::ConstantFocusingRep(const std::string& name)
    : ConstantFocusing(name), geometry() {}

ConstantFocusingRep::~ConstantFocusingRep() {}

ElementBase* ConstantFocusingRep::clone() const { return new ConstantFocusingRep(*this); }

Geometry& ConstantFocusingRep::getGeometry() { return geometry; }

const Geometry& ConstantFocusingRep::getGeometry() const { return geometry; }
