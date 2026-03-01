#include "BeamlineCore/ConstantEzRep.h"
#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (ConstantEzRep::*get)() const;
        void (ConstantEzRep::*set)(double);
    };

    const Entry entries[] = {
        {"L", &ConstantEzRep::getElementLength, &ConstantEzRep::setElementLength},
        {"EZ", &ConstantEzRep::getEz, &ConstantEzRep::setEz},
        {0, 0, 0}
    };
}

ConstantEzRep::ConstantEzRep() : ConstantEz(), geometry(), field() {}

ConstantEzRep::ConstantEzRep(const ConstantEzRep& right)
    : ConstantEz(right),
      geometry(right.geometry),
      field(right.field) {}

ConstantEzRep::ConstantEzRep(const std::string& name)
    : ConstantEz(name),
      geometry(),
      field() {}

ConstantEzRep::~ConstantEzRep() {}

ElementBase* ConstantEzRep::clone() const {
    return new ConstantEzRep(*this);
}

Channel* ConstantEzRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != 0; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<ConstantEzRep>(*this, entry->get, entry->set);
        }
    }
    return ElementBase::getChannel(aKey, create);
}

ConstEzField& ConstantEzRep::getField() {
    return field;
}

const ConstEzField& ConstantEzRep::getField() const {
    return field;
}

StraightGeometry& ConstantEzRep::getGeometry() {
    return geometry;
}

const StraightGeometry& ConstantEzRep::getGeometry() const {
    return geometry;
}

void ConstantEzRep::setElementLength(double length) {
    geometry.setElementLength(length);
}

void ConstantEzRep::setEz(double ez) {
    ConstantEz::setEz(ez);
    field.setEz(ez);
}
