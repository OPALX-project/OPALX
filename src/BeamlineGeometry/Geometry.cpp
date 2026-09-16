//
// Class Geometry
//   Single geometry class for all beamline elements.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "BeamlineGeometry/Geometry.h"

#include "Algorithms/Quaternion.hpp"

#include <algorithm>
#include <cmath>

namespace {
    /// Rotation matrix for a rotation by `angle` about the local y-axis. Matches
    /// the convention of the former Rotation3D::YRotation used by the bend frames.
    matrix3x3_t yRotationMatrix(double angle) {
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        matrix3x3_t m;
        m(0, 0) = c;
        m(0, 1) = 0.0;
        m(0, 2) = s;
        m(1, 0) = 0.0;
        m(1, 1) = 1.0;
        m(1, 2) = 0.0;
        m(2, 0) = -s;
        m(2, 1) = 0.0;
        m(2, 2) = c;
        return m;
    }

    /// CoordinateSystemTrafo for a body frame located at `position` and rotated by
    /// `yAngle` about the y-axis. Equivalent to the former
    /// toCoordinateSystemTrafo(Euclid3D{position, YRotation(yAngle)}).
    CoordinateSystemTrafo frameTrafo(const Vector_t<double, 3>& position, double yAngle) {
        return CoordinateSystemTrafo(position, Quaternion(yRotationMatrix(yAngle)).conjugate());
    }

    CoordinateSystemTrafo identityTrafo(const Vector_t<double, 3>& position) {
        return CoordinateSystemTrafo(position, Quaternion(1, 0, 0, 0));
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

Geometry Geometry::makeSBend(double length, double curvature) {
    Geometry g;
    g.kind_m  = GeometryKind::SBend;
    g.len_m   = length;
    g.h_m     = curvature;
    g.angle_m = curvature * length;
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
        default:  // Straight, SBend
            return len_m;
    }
}

void Geometry::setElementLength(double length) {
    if (kind_m == GeometryKind::SBend) {
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
        case GeometryKind::SBend:
        case GeometryKind::RBend:
            return angle_m;
        default:
            return 0.0;
    }
}

void Geometry::setBendAngle(double angle) {
    angle_m = angle;
    if (kind_m == GeometryKind::SBend && len_m != 0.0) {
        h_m = angle_m / len_m;
    }
}

void Geometry::setCurvature(double curvature) {
    if (kind_m == GeometryKind::SBend && len_m != 0.0) {
        h_m     = curvature;
        angle_m = h_m * len_m;
    }
}

Vector_t<double, 3> Geometry::framePosition(double s) const {
    if (kind_m == GeometryKind::SBend && h_m != 0.0) {
        const double phi = h_m * s;
        return Vector_t<double, 3>({(std::cos(phi) - 1.0) / h_m, 0.0, std::sin(phi) / h_m});
    }
    return Vector_t<double, 3>({0.0, 0.0, s});  // straight / rbend body / null / sbend(h=0)
}

double Geometry::getChordLength() const {
    const Vector_t<double, 3> delta = framePosition(len_m) - framePosition(0.0);
    return std::sqrt(delta(0) * delta(0) + delta(1) * delta(1) + delta(2) * delta(2));
}

std::vector<Vector_t<double, 3>> Geometry::getDesignPath(std::size_t minSamples) const {
    const double sBegin = 0.0;    // entrance edge (chart anchored at the entrance)
    const double sEnd   = len_m;  // exit edge
    const double span   = std::abs(sEnd - sBegin);
    const std::size_t samples =
            std::max<std::size_t>(minSamples, static_cast<std::size_t>(std::ceil(span / 0.01)) + 1);

    std::vector<Vector_t<double, 3>> path;
    path.reserve(samples);
    for (std::size_t i = 0; i < samples; ++i) {
        const double alpha =
                (samples > 1) ? static_cast<double>(i) / static_cast<double>(samples - 1) : 0.0;
        const double s = sBegin + alpha * (sEnd - sBegin);
        path.emplace_back(framePosition(s));
    }

    return path;
}

CoordinateSystemTrafo Geometry::getEdgeToBegin() const {
    // The stored (global-to-local) frame IS the element's geometrical entrance face for every
    // kind (sector bend: entrance tangent; rectangular bend: box/chord face; straight body:
    // entrance edge), so the entrance-edge transform is always the identity.
    return identityTrafo(Vector_t<double, 3>({0.0, 0.0, 0.0}));
}

CoordinateSystemTrafo Geometry::getEdgeToEnd() const {
    // Transform from the entrance face (the stored frame) to the exit face of the body. A
    // straight body — straight element, rectangular-bend box, or null — has parallel faces, so
    // this is a pure +z shift by the body length. A sector bend's body curves, so its exit face
    // sits at the arc end (framePosition(len)) with the tangent turned by the full bend angle.
    // (The reference orbit meets an RBend's faces at half the bend angle, but that orbit tangent
    // is not part of the body edge geometry.)
    switch (kind_m) {
        case GeometryKind::SBend:
            return frameTrafo(framePosition(len_m), -h_m * len_m);
        default:  // Straight, RBend, Null: exit face parallel to the entrance face
            return identityTrafo(Vector_t<double, 3>({0.0, 0.0, len_m}));
    }
}
