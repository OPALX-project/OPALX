#include "AbsBeamline/RBend.h"

#include "AbsBeamline/BeamlineVisitor.h"

namespace {
    CoordinateSystemTrafo toCoordinateSystemTrafo(const Euclid3D& frame) {
        matrix3x3_t rotation;
        const Rotation3D& euclidRotation = frame.getRotation();
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                rotation(row, col) = euclidRotation(row, col);
            }
        }

        const Vector3D& displacement = frame.getVector();
        const Vector_t<double, 3> origin(
                displacement.getX(), displacement.getY(), displacement.getZ());

        return CoordinateSystemTrafo(origin, Quaternion(rotation).conjugate());
    }
}  // namespace

RBend::RBend() : RBend("") {}

RBend::RBend(const RBend& right) : BendBase(right) {}

RBend::RBend(const std::string& name) : BendBase(name) {}

RBend::~RBend() = default;

void RBend::accept(BeamlineVisitor& visitor) const { visitor.visitRBend(*this); }

ElementType RBend::getType() const { return ElementType::RBEND; }

CoordinateSystemTrafo RBend::getEdgeToBegin() const {
    const double halfAngle = 0.5 * getBendAngle();
    const double xMid      = -0.5 * getElementLength() * std::tan(0.5 * halfAngle);
    const Euclid3D bodyToReferenceEntry =
            Euclid3D::translation(xMid, 0.0, 0.5 * getElementLength()) * Euclid3D::YRotation(halfAngle);
    return toCoordinateSystemTrafo(bodyToReferenceEntry);
}

CoordinateSystemTrafo RBend::getEdgeToEnd() const {
    const double halfAngle = 0.5 * getBendAngle();
    const double xMid      = -0.5 * getElementLength() * std::tan(0.5 * halfAngle);
    const Euclid3D bodyToReferenceEntry =
            Euclid3D::translation(xMid, 0.0, 0.5 * getElementLength()) * Euclid3D::YRotation(halfAngle);
    const Euclid3D referenceEntryToExit =
            Euclid3D::translation(-getElementLength() * std::tan(halfAngle), 0.0, getElementLength())
            * Euclid3D::YRotation(-getBendAngle());
    const Euclid3D bodyToReferenceExit = bodyToReferenceEntry * referenceEntryToExit;
    return toCoordinateSystemTrafo(bodyToReferenceExit);
}

double RBend::getExitAngle() const { return getBendAngle() - getEntranceAngle(); }

double RBend::getReferencePathLength() const { return getGeometry().getArcLength(); }
