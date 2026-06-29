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
        m(0, 0) = c;    m(0, 1) = 0.0; m(0, 2) = s;
        m(1, 0) = 0.0;  m(1, 1) = 1.0; m(1, 2) = 0.0;
        m(2, 0) = -s;   m(2, 1) = 0.0; m(2, 2) = c;
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

Geometry Geometry::makeArc(double length, double curvature) {
    Geometry g;
    g.kind_m  = GeometryKind::Arc;
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

Vector_t<double, 3> Geometry::framePosition(double s) const {
    if (kind_m == GeometryKind::Arc && h_m != 0.0) {
        const double phi = h_m * s;
        return Vector_t<double, 3>({(std::cos(phi) - 1.0) / h_m, 0.0, std::sin(phi) / h_m});
    }
    return Vector_t<double, 3>({0.0, 0.0, s});  // straight / rbend body / null / arc(h=0)
}

double Geometry::getChordLength() const {
    const Vector_t<double, 3> delta = framePosition(getExit()) - framePosition(getEntrance());
    return std::sqrt(delta(0) * delta(0) + delta(1) * delta(1) + delta(2) * delta(2));
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
        const double s = sBegin + alpha * (sEnd - sBegin);
        path.emplace_back(framePosition(s));
    }

    return path;
}

CoordinateSystemTrafo Geometry::getEdgeToBegin() const {
    switch (kind_m) {
        case GeometryKind::Arc:
            // Frame rotation is YRotation(-phi) with phi = h * (-len/2) = -angle/2.
            return frameTrafo(framePosition(getEntrance()), h_m * (len_m / 2.0));
        case GeometryKind::RBend:
            return frameTrafo(Vector_t<double, 3>({0.0, 0.0, -len_m / 2.0}), halfAngle());
        default:  // Straight, Null: identity at the entrance edge
            return identityTrafo(Vector_t<double, 3>({0.0, 0.0, 0.0}));
    }
}

CoordinateSystemTrafo Geometry::getEdgeToEnd() const {
    switch (kind_m) {
        case GeometryKind::Arc:
            return frameTrafo(framePosition(getExit()), -h_m * (len_m / 2.0));
        case GeometryKind::RBend:
            return frameTrafo(Vector_t<double, 3>({0.0, 0.0, len_m / 2.0}), -halfAngle());
        default:  // Straight, Null: +z shift to the exit edge
            return identityTrafo(Vector_t<double, 3>({0.0, 0.0, len_m}));
    }
}
