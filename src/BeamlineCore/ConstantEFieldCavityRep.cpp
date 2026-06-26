#include "BeamlineCore/ConstantEFieldCavityRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (ConstantEFieldCavityRep::*get)() const;
        void (ConstantEFieldCavityRep::*set)(double);
    };

    const Entry entries[] = {
            {"L", &ConstantEFieldCavityRep::getElementLength,
             &ConstantEFieldCavityRep::setElementLength},
            {"EZ", &ConstantEFieldCavityRep::getEz, &ConstantEFieldCavityRep::setEz},
            {0, 0, 0}};
}  // namespace

ConstantEFieldCavityRep::ConstantEFieldCavityRep() : ConstantEFieldCavity(), geometry(), field() {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const ConstantEFieldCavityRep& right)
    : ConstantEFieldCavity(right), geometry(right.geometry), field(right.field) {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const std::string& name)
    : ConstantEFieldCavity(name), geometry(), field() {}

ConstantEFieldCavityRep::~ConstantEFieldCavityRep() {}

ElementBase* ConstantEFieldCavityRep::clone() const { return new ConstantEFieldCavityRep(*this); }

Channel* ConstantEFieldCavityRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<ConstantEFieldCavityRep>(*this, entry->get, entry->set);
        }
    }
    return ElementBase::getChannel(aKey, create);
}

ConstEzField& ConstantEFieldCavityRep::getField() { return field; }

const ConstEzField& ConstantEFieldCavityRep::getField() const { return field; }

StraightGeometry& ConstantEFieldCavityRep::getGeometry() { return geometry; }

const StraightGeometry& ConstantEFieldCavityRep::getGeometry() const { return geometry; }

void ConstantEFieldCavityRep::setElementLength(double length) { geometry.setElementLength(length); }

void ConstantEFieldCavityRep::setEz(double ez) {
    ConstantEFieldCavity::setEz(ez);
    field.setEz(ez);
}
