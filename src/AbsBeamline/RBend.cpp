#include "AbsBeamline/RBend.h"

#include "AbsBeamline/BeamlineVisitor.h"

RBend::RBend() : RBend("") {}

RBend::RBend(const RBend& right) : BendBase(right) {}

RBend::RBend(const std::string& name) : BendBase(name) {}

RBend::~RBend() = default;

void RBend::accept(BeamlineVisitor& visitor) const { visitor.visitComponent(*this); }

ElementType RBend::getType() const { return ElementType::RBEND; }

double RBend::getExitAngle() const { return getBendAngle() - getEntranceAngle(); }

double RBend::getReferencePathLength() const { return getGeometry().getArcLength(); }
