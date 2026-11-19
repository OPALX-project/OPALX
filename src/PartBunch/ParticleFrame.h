#ifndef OPALX_PARTICLE_FRAME_H
#define OPALX_PARTICLE_FRAME_H

#include "Algorithms/CoordinateSystemTrafo.h"

class ElementBase;

/**
 * @class ParticleFrame
 * @brief A rigid frame for the particle data of one container, relative to its reference frame.
 *
 * The reference frame of a container has its origin at the reference particle and z along its
 * momentum; ParticleContainer::getToLabTrafo() maps it to the lab. Every other frame stores the
 * transform from the reference frame into it and the transform back. Keeping both directions
 * lets a change through the reference frame apply exactly the transforms the caller built,
 * without inverting a transform twice.
 *
 * A frame is only valid for the reference frame it was built from. Build it where it is
 * requested and do not keep it across a reference update.
 */
class ParticleFrame {
public:
    enum class Kind : unsigned char { Reference, Element, Beamline, SpaceChargeSolve };

    /// The container's reference frame.
    ParticleFrame() = default;

    /**
     * @param fromReference Maps reference-frame coordinates into this frame.
     * @param toReference   Inverse of @p fromReference.
     * @param element       Element whose local frame this is, if any.
     */
    ParticleFrame(
            Kind kind, const CoordinateSystemTrafo& fromReference,
            const CoordinateSystemTrafo& toReference, const ElementBase* element = nullptr)
        : kind_m(kind),
          element_m(element),
          fromReference_m(fromReference),
          toReference_m(toReference) {}

    ParticleFrame(
            Kind kind, const CoordinateSystemTrafo& fromReference,
            const ElementBase* element = nullptr)
        : ParticleFrame(kind, fromReference, fromReference.inverted(), element) {}

    Kind kind() const { return kind_m; }
    const char* name() const {
        switch (kind_m) {
            case Kind::Reference:
                return "reference";
            case Kind::Element:
                return "element";
            case Kind::Beamline:
                return "beamline";
            case Kind::SpaceChargeSolve:
                return "space-charge solve";
        }
        return "unknown";
    }
    bool isReference() const { return kind_m == Kind::Reference; }
    const ElementBase* element() const { return element_m; }
    const CoordinateSystemTrafo& fromReference() const { return fromReference_m; }
    const CoordinateSystemTrafo& toReference() const { return toReference_m; }

    /**
     * @brief Transform from this frame to @p target, composed on the host.
     *
     * A change from or to the reference frame returns the stored transform itself, so it is
     * bitwise the transform the caller built. Only a change between two other frames composes.
     */
    CoordinateSystemTrafo to(const ParticleFrame& target) const {
        if (isReference()) return target.fromReference_m;
        if (target.isReference()) return toReference_m;
        return target.fromReference_m * toReference_m;
    }

    /// Same kind, element and transform; then a change between the two is the identity.
    bool operator==(const ParticleFrame& other) const {
        if (kind_m != other.kind_m || element_m != other.element_m) return false;
        if (isReference()) return true;
        const auto a = fromReference_m.getOrigin(), b = other.fromReference_m.getOrigin();
        const auto p = fromReference_m.getRotation(), q = other.fromReference_m.getRotation();
        return a(0) == b(0) && a(1) == b(1) && a(2) == b(2) && p(0) == q(0) && p(1) == q(1)
               && p(2) == q(2) && p(3) == q(3);
    }

private:
    Kind kind_m                  = Kind::Reference;
    const ElementBase* element_m = nullptr;
    CoordinateSystemTrafo fromReference_m;
    CoordinateSystemTrafo toReference_m;
};

#endif
