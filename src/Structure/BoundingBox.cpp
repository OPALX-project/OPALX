//
// Class BoundingBox
//
// This class provides functionality to compute bounding boxes, to compute if a position
// is inside the box and to compute the intersection point between a ray and the bounding box.
//
// Copyright (c) 201x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#include "BoundingBox.h"

#include <iomanip>
#include <limits>

BoundingBox::BoundingBox():
    lowerLeftCorner_m(std::numeric_limits<double>::max()),
    upperRightCorner_m(std::numeric_limits<double>::lowest())
{}

BoundingBox BoundingBox::getBoundingBox(const std::vector<Vector_t<double, 3>>& positions)
{
    BoundingBox boundingBox;
    for (const Vector_t<double, 3>& position : positions) {
        boundingBox.enlargeToContainPosition(position);
    }

    return boundingBox;
}

void BoundingBox::enlargeToContainPosition(const Vector_t<double, 3>& position)
{
    for (unsigned int d = 0; d < 3; ++ d) {
        lowerLeftCorner_m[d] = std::min(lowerLeftCorner_m[d], position[d]);
        upperRightCorner_m[d] = std::max(upperRightCorner_m[d], position[d]);
    }
}

void BoundingBox::enlargeToContainBoundingBox(const BoundingBox& boundingBox)
{
    for (unsigned int d = 0; d < 3; ++ d) {
        lowerLeftCorner_m[d] = std::min(lowerLeftCorner_m[d], boundingBox.lowerLeftCorner_m[d]);
        upperRightCorner_m[d] = std::max(upperRightCorner_m[d], boundingBox.upperRightCorner_m[d]);
    }
}

boost::optional<Vector_t<double, 3>> BoundingBox::getIntersectionPoint(const Vector_t<double, 3>& position,
                                                            const Vector_t<double, 3>& direction) const
{
    boost::optional<Vector_t<double, 3>> result = boost::none;
    double minDistance = std::numeric_limits<double>::max();
    const Vector_t<double, 3> dimensions = upperRightCorner_m - lowerLeftCorner_m;
    Vector_t<double, 3> normal(1, 0, 0);
    for (unsigned int d = 0; d < 3; ++ d) {
        double sign = -1;
        Vector_t<double, 3> upperCorner = lowerLeftCorner_m + dot(normal, upperRightCorner_m) * normal;
        for (const Vector_t<double, 3>& p0 : {lowerLeftCorner_m, upperCorner}) {
            double tau = dot(p0 - position, sign * normal) / dot(direction, sign * normal);
            if (tau < 0.0) {
                continue;
            }
            Vector_t<double, 3> pointOnPlane = position + tau * direction;
            Vector_t<double, 3> relativeP = pointOnPlane - p0;
            bool isOnFace = true;
            for (unsigned int i = 1; i < 3; ++ i) {
                unsigned int idx = (d + i) % 3;
                if (relativeP[idx] < 0 ||
                    relativeP[idx] > dimensions[idx]) {
                    isOnFace = false;
                    break;
                }
            }
            if (isOnFace) {
                double distance = euclidean_norm(pointOnPlane - position);
                if (distance < minDistance) {
                    minDistance = distance;
                    result = pointOnPlane;
                }
            }
            sign *= -1;
        }

        normal = Vector_t<double, 3>(normal[2], normal[0], normal[1]);
    }

    return result;
}

bool BoundingBox::isInside(const Vector_t<double, 3>& position) const
{
    Vector_t<double, 3> relPosition = position - lowerLeftCorner_m;
    Vector_t<double, 3> dimensions = upperRightCorner_m - lowerLeftCorner_m;
    for (unsigned int d = 0; d < 3; ++ d) {
        if (relPosition[d] < 0 || relPosition[d] > dimensions[d]) return false;
    }
    return true;
}

void BoundingBox::print(std::ostream& output) const
{
    int prePrecision = output.precision();
    output << std::setprecision(8);
    output << "Bounding box\n"
           << "lower left corner: " << lowerLeftCorner_m << "\n"
           << "upper right corner: " << upperRightCorner_m << std::endl;
    output.precision(prePrecision);
}
