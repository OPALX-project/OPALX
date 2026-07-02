#include "BeamlineCore/ConstantEFieldCavityRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (ConstantEFieldCavityRep::*get)() const;
        void (ConstantEFieldCavityRep::*set)(double);
    };

    const Entry entries[] = {
                        {"EZ", &ConstantEFieldCavityRep::getEz, &ConstantEFieldCavityRep::setEz},
            {0, 0, 0}};
}  // namespace

ConstantEFieldCavityRep::ConstantEFieldCavityRep() : ConstantEFieldCavity(), geometry() {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const ConstantEFieldCavityRep& right)
    : ConstantEFieldCavity(right), geometry(right.geometry) {}

ConstantEFieldCavityRep::ConstantEFieldCavityRep(const std::string& name)
    : ConstantEFieldCavity(name), geometry() {}

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

Geometry& ConstantEFieldCavityRep::getGeometry() { return geometry; }

const Geometry& ConstantEFieldCavityRep::getGeometry() const { return geometry; }

