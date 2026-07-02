#include "BeamlineCore/RBendRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (RBendRep::*get)() const;
        void (RBendRep::*set)(double);
    };

    const Entry entries[] = {
            {"BY", &RBendRep::getB, &RBendRep::setB},
            {nullptr, nullptr, nullptr}};
}  // namespace

RBendRep::RBendRep() : RBend(), geometry_m(Geometry::makeRBend(0.0, 0.0)) {}

RBendRep::RBendRep(const RBendRep& right) : RBend(right), geometry_m(right.geometry_m) {}

RBendRep::RBendRep(const std::string& name)
    : RBend(name), geometry_m(Geometry::makeRBend(0.0, 0.0)) {}

RBendRep::~RBendRep() = default;

ElementBase* RBendRep::clone() const { return new RBendRep(*this); }

Channel* RBendRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != nullptr; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<RBendRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

Geometry& RBendRep::getGeometry() { return geometry_m; }

const Geometry& RBendRep::getGeometry() const { return geometry_m; }
