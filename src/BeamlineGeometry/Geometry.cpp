// ------------------------------------------------------------------------
// $RCSfile: BGeometryBase.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BGeometryBase
//    Pure virtual base class for all Beamline Geometries
//
// ------------------------------------------------------------------------
// Class category: BeamlineBGeometryBase
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/Euclid3D.h"

#include "Algorithms/Quaternion.hpp"

#include <algorithm>
#include <cmath>

// Class BGeometryBase.
// ------------------------------------------------------------------------

BGeometryBase::~BGeometryBase() {}

void BGeometryBase::setElementLength(double) {}

double BGeometryBase::getOrigin() const { return getArcLength() / 2.0; }

double BGeometryBase::getEntrance() const { return -getOrigin(); }

double BGeometryBase::getExit() const { return getArcLength() - getOrigin(); }

Euclid3D BGeometryBase::getTotalTransform() const { return getTransform(getExit(), getEntrance()); }

Euclid3D BGeometryBase::getTransform(double s) const { return getTransform(0.0, s); }

Euclid3D BGeometryBase::getEntranceFrame() const { return getTransform(0.0, getEntrance()); }

Euclid3D BGeometryBase::getExitFrame() const { return getTransform(0.0, getExit()); }

Euclid3D BGeometryBase::getEntrancePatch() const { return Euclid3D::identity(); }

Euclid3D BGeometryBase::getExitPatch() const { return Euclid3D::identity(); }

// Class Geometry.
// ------------------------------------------------------------------------

namespace {
    /// General transformation along an arc of length l and curvature h (XZ plane).
    /// Reproduces PlanarArcGeometry's ArcTransform.
    Euclid3D arcTransform(double l, double h) {
        Euclid3D t;
        if (h) {
            double phi = h * l;
            t          = Euclid3D::YRotation(-phi);
            t.setX((std::cos(phi) - 1.0) / h);
            t.setZ(std::sin(phi) / h);
        } else {
            t.setZ(l);
        }
        return t;
    }

    /// Convert a legacy Euclid3D body frame into the modern CoordinateSystemTrafo
    /// used for placement. (Formerly the file-static helper in BendBase.cpp.)
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

Geometry Geometry::makeNull() {
    Geometry g;
    g.kind_m = GeometryKind::Null;
    return g;
}

Geometry Geometry::makeStraight(double length) {
    Geometry g;
    g.kind_m = GeometryKind::Straight;
    g.len_m  = std::max(0.0, length);
    return g;
}

Geometry Geometry::makeArc(double length, double curvature) {
    Geometry g;
    g.kind_m   = GeometryKind::Arc;
    g.len_m    = length;
    g.h_m      = curvature;
    g.angle_m  = curvature * length;
    return g;
}

Geometry Geometry::makeRBend(double bodyLength, double angle) {
    Geometry g;
    g.kind_m  = GeometryKind::RBend;
    g.len_m   = std::max(0.0, bodyLength);
    g.angle_m = angle;
    return g;
}

double Geometry::getArcLength() const {
    switch (kind_m) {
        case GeometryKind::Null:
            return 0.0;
        case GeometryKind::RBend: {
            const double ha = halfAngle();
            return (ha == 0.0) ? len_m : len_m * ha / std::sin(ha);
        }
        default:  // Straight, Arc
            return len_m;
    }
}

void Geometry::setElementLength(double length) {
    if (kind_m == GeometryKind::Arc) {
        len_m = length;
        if (len_m != 0.0) {
            angle_m = h_m * len_m;
        }
    } else {  // Straight, RBend, Null
        len_m = std::max(0.0, length);
    }
}

double Geometry::getBendAngle() const {
    switch (kind_m) {
        case GeometryKind::Arc:
        case GeometryKind::RBend:
            return angle_m;
        default:
            return 0.0;
    }
}

void Geometry::setBendAngle(double angle) {
    angle_m = angle;
    if (kind_m == GeometryKind::Arc && len_m != 0.0) {
        h_m = angle_m / len_m;
    }
}

void Geometry::setCurvature(double curvature) {
    if (kind_m == GeometryKind::Arc && len_m != 0.0) {
        h_m     = curvature;
        angle_m = h_m * len_m;
    }
}

Euclid3D Geometry::getTransform(double s) const {
    if (kind_m == GeometryKind::Arc) {
        return arcTransform(s, h_m);
    }
    return Euclid3D::translation(0.0, 0.0, s);  // Straight, RBend body, Null
}

Euclid3D Geometry::getTransform(double fromS, double toS) const {
    if (kind_m == GeometryKind::Arc) {
        return arcTransform(toS - fromS, h_m);  // matches PlanarArcGeometry
    }
    return Euclid3D::translation(0.0, 0.0, fromS - toS);  // matches StraightGeometry
}

Euclid3D Geometry::getTotalTransform() const {
    switch (kind_m) {
        case GeometryKind::Arc:
            return arcTransform(len_m, h_m);
        case GeometryKind::RBend: {
            const Euclid3D patch = Euclid3D::YRotation(-halfAngle());
            const Euclid3D body  = Euclid3D::translation(0.0, 0.0, len_m);
            return patch * body * patch;  // matches RBendGeometry
        }
        case GeometryKind::Null:
            return Euclid3D::identity();
        default:  // Straight
            return Euclid3D::translation(0.0, 0.0, len_m);
    }
}

Euclid3D Geometry::getEntranceFrame() const {
    switch (kind_m) {
        case GeometryKind::Arc:
            return arcTransform(-len_m / 2.0, h_m);
        case GeometryKind::RBend:
            return Euclid3D::translation(0.0, 0.0, -len_m / 2.0) * Euclid3D::YRotation(halfAngle());
        case GeometryKind::Null:
            return Euclid3D::identity();
        default:  // Straight
            return Euclid3D::translation(0.0, 0.0, -len_m / 2.0);
    }
}

Euclid3D Geometry::getExitFrame() const {
    switch (kind_m) {
        case GeometryKind::Arc:
            return arcTransform(len_m / 2.0, h_m);
        case GeometryKind::RBend:
            return Euclid3D::translation(0.0, 0.0, len_m / 2.0) * Euclid3D::YRotation(-halfAngle());
        case GeometryKind::Null:
            return Euclid3D::identity();
        default:  // Straight
            return Euclid3D::translation(0.0, 0.0, len_m / 2.0);
    }
}

double Geometry::getChordLength() const {
    const Vector3D delta = getExitFrame().getVector() - getEntranceFrame().getVector();
    return std::sqrt(
            delta.getX() * delta.getX() + delta.getY() * delta.getY()
            + delta.getZ() * delta.getZ());
}

std::vector<Vector_t<double, 3>> Geometry::getDesignPath(std::size_t minSamples) const {
    const double sBegin = getEntrance();
    const double sEnd   = getExit();
    const double span   = std::abs(sEnd - sBegin);
    const std::size_t samples =
            std::max<std::size_t>(minSamples, static_cast<std::size_t>(std::ceil(span / 0.01)) + 1);

    std::vector<Vector_t<double, 3>> path;
    path.reserve(samples);
    for (std::size_t i = 0; i < samples; ++i) {
        const double alpha =
                (samples > 1) ? static_cast<double>(i) / static_cast<double>(samples - 1) : 0.0;
        const double s      = sBegin + alpha * (sEnd - sBegin);
        const Euclid3D pose = getTransform(s);
        path.emplace_back(pose.getX(), pose.getY(), pose.getZ());
    }

    return path;
}

CoordinateSystemTrafo Geometry::getEdgeToBegin() const {
    switch (kind_m) {
        case GeometryKind::Arc:
        case GeometryKind::RBend:
            return toCoordinateSystemTrafo(getEntranceFrame());
        default:  // Straight, Null: identity at the entrance edge
            return CoordinateSystemTrafo(
                    Vector_t<double, 3>({0.0, 0.0, 0.0}), Quaternion(1, 0, 0, 0));
    }
}

CoordinateSystemTrafo Geometry::getEdgeToEnd() const {
    switch (kind_m) {
        case GeometryKind::Arc:
        case GeometryKind::RBend:
            return toCoordinateSystemTrafo(getExitFrame());
        default:  // Straight, Null: +z shift to the exit edge
            return CoordinateSystemTrafo(
                    Vector_t<double, 3>({0.0, 0.0, len_m}), Quaternion(1, 0, 0, 0));
    }
}
