#include "BeamlineCore/RBendRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (RBendRep::*get)() const;
        void (RBendRep::*set)(double);
    };

    const Entry entries[] = {
            {"L", &RBendRep::getElementLength, &RBendRep::setElementLength},
            {"BY", &RBendRep::getB, &RBendRep::setB},
            {"E1", &RBendRep::getEntryFaceRotation, &RBendRep::setEntryFaceRotation},
            {"E2", &RBendRep::getExitFaceRotation, &RBendRep::setExitFaceRotation},
            {"H1", &RBendRep::getEntryFaceCurvature, &RBendRep::setEntryFaceCurvature},
            {"H2", &RBendRep::getExitFaceCurvature, &RBendRep::setExitFaceCurvature},
            {nullptr, nullptr, nullptr}};
}  // namespace

RBendRep::RBendRep() : RBend(), geometry_m(0.0, 0.0), field_m() {}

RBendRep::RBendRep(const RBendRep& right)
    : RBend(right), geometry_m(right.geometry_m), field_m(right.field_m) {}

RBendRep::RBendRep(const std::string& name) : RBend(name), geometry_m(0.0, 0.0), field_m() {}

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

BMultipoleField& RBendRep::getField() { return field_m; }

const BMultipoleField& RBendRep::getField() const { return field_m; }

RBendGeometry& RBendRep::getGeometry() { return geometry_m; }

const RBendGeometry& RBendRep::getGeometry() const { return geometry_m; }

void RBendRep::setField(const BMultipoleField& field) { field_m = field; }
