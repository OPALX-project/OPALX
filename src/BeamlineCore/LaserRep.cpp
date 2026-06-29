#include "BeamlineCore/LaserRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (LaserRep::*get)() const;
        void (LaserRep::*set)(double);
    };

    const Entry entries[] = {
            {"L", &LaserRep::getElementLength, &LaserRep::setElementLength},
            {nullptr, nullptr, nullptr}};
}  // namespace

LaserRep::LaserRep() : Laser(), field_m(), geometry_m(0.0) {}

LaserRep::LaserRep(const LaserRep& right) : Laser(right), field_m(), geometry_m(right.geometry_m) {}

LaserRep::LaserRep(const std::string& name) : Laser(name), field_m(), geometry_m(0.0) {}

LaserRep::~LaserRep() {}

ElementBase* LaserRep::clone() const { return new LaserRep(*this); }

Channel* LaserRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != nullptr; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<LaserRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

NullField& LaserRep::getField() { return field_m; }

const NullField& LaserRep::getField() const { return field_m; }

StraightGeometry& LaserRep::getGeometry() { return geometry_m; }

const StraightGeometry& LaserRep::getGeometry() const { return geometry_m; }
