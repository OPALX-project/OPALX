#include "BeamlineCore/SBendRep.h"

#include "Channels/IndirectChannel.h"

namespace {
    struct Entry {
        const char* name;
        double (SBendRep::*get)() const;
        void (SBendRep::*set)(double);
    };

    const Entry entries[] = {
            {"L", &SBendRep::getElementLength, &SBendRep::setElementLength},
            {"BY", &SBendRep::getB, &SBendRep::setB},
            {"E1", &SBendRep::getEntryFaceRotation, &SBendRep::setEntryFaceRotation},
            {"E2", &SBendRep::getExitFaceRotation, &SBendRep::setExitFaceRotation},
            {"H1", &SBendRep::getEntryFaceCurvature, &SBendRep::setEntryFaceCurvature},
            {"H2", &SBendRep::getExitFaceCurvature, &SBendRep::setExitFaceCurvature},
            {nullptr, nullptr, nullptr}};
}  // namespace

SBendRep::SBendRep() : SBend(), geometry_m(Geometry::makeArc(0.0, 0.0)) {}

SBendRep::SBendRep(const SBendRep& right) : SBend(right), geometry_m(right.geometry_m) {}

SBendRep::SBendRep(const std::string& name) : SBend(name), geometry_m(Geometry::makeArc(0.0, 0.0)) {}

SBendRep::~SBendRep() = default;

ElementBase* SBendRep::clone() const { return new SBendRep(*this); }

Channel* SBendRep::getChannel(const std::string& aKey, bool create) {
    for (const Entry* entry = entries; entry->name != nullptr; ++entry) {
        if (aKey == entry->name) {
            return new IndirectChannel<SBendRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}

BGeometryBase& SBendRep::getGeometry() { return geometry_m; }

const BGeometryBase& SBendRep::getGeometry() const { return geometry_m; }
