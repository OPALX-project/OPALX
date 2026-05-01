#include "AbsBeamline/SBend.h"

#include "AbsBeamline/BeamlineVisitor.h"

SBend::SBend() : SBend("") {}

SBend::SBend(const SBend& right) : BendBase(right) {}

SBend::SBend(const std::string& name) : BendBase(name) {}

SBend::~SBend() = default;

void SBend::accept(BeamlineVisitor& visitor) const { visitor.visitSBend(*this); }

ElementType SBend::getType() const { return ElementType::SBEND; }

void SBend::setExitAngle(double exitAngle) { BendBase::setExitAngle(exitAngle); }

double SBend::getExitAngle() const { return getStoredExitAngle(); }
