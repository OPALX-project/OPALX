#include "BeamlineCore/LaserRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (LaserRep::*get)() const;
        void (LaserRep::*set)(double);
    };

    const Entry entries[] = {{nullptr, nullptr, nullptr}};
}  // namespace

LaserRep::LaserRep() : Laser(), geometry_m() {}

LaserRep::LaserRep(const LaserRep& right) : Laser(right), geometry_m(right.geometry_m) {}

LaserRep::LaserRep(const std::string& name) : Laser(name), geometry_m() {}

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

Geometry& LaserRep::getGeometry() { return geometry_m; }

const Geometry& LaserRep::getGeometry() const { return geometry_m; }
