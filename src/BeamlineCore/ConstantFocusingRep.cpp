#include "BeamlineCore/ConstantFocusingRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (ConstantFocusingRep::*get)() const;
        void (ConstantFocusingRep::*set)(double);
    };

    const Entry entries[] = {
            {"STRENGTH", &ConstantFocusingRep::getStrength, &ConstantFocusingRep::setStrength},
            {"RADIUS", &ConstantFocusingRep::getRadius, &ConstantFocusingRep::setRadius},
            {nullptr, nullptr, nullptr}};
}  // namespace

ConstantFocusingRep::ConstantFocusingRep() : ConstantFocusing(), geometry() {}

ConstantFocusingRep::ConstantFocusingRep(const ConstantFocusingRep& right)
    : ConstantFocusing(right), geometry(right.geometry) {}

ConstantFocusingRep::ConstantFocusingRep(const std::string& name)
    : ConstantFocusing(name), geometry() {}

ConstantFocusingRep::~ConstantFocusingRep() {}

ElementBase* ConstantFocusingRep::clone() const { return new ConstantFocusingRep(*this); }

Channel* ConstantFocusingRep::getChannel(const std::string& key, bool create) {
    for (const Entry* entry = entries; entry->name != nullptr; ++entry) {
        if (key == entry->name) {
            return new IndirectChannel<ConstantFocusingRep>(*this, entry->get, entry->set);
        }
    }
    return ElementBase::getChannel(key, create);
}

Geometry& ConstantFocusingRep::getGeometry() { return geometry; }

const Geometry& ConstantFocusingRep::getGeometry() const { return geometry; }
