#ifndef OPALX_PlacementResolver_HH
#define OPALX_PlacementResolver_HH

//
// Class PlacementResolver
//   The single PLACE stage of the input -> computeExternalFields() pipeline.
//
// Copyright (c) 2024, Paul Scherrer Institut, Villigen PSI, Switzerland
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

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Utilities/FieldList.h"

/**
 * @class PlacementResolver
 * @brief Resolves every beamline element's global-to-local transform in one place.
 *
 * There are two ways to place an element, and both are resolved here so placement has a single
 * owner:
 *   - 6D pose (Mode A): the element recorded an absolute global-to-local pose; it is composed
 *     with the lab frame and fixed.
 *   - ELEMEDGE (Mode B): the element is positioned by walking the reference path in s-order.
 *
 * @note The caller owns the "already resolved" guard; resolve() does the work unconditionally.
 */
class PlacementResolver {
public:
    /**
     * @brief Resolve and fix the global-to-local transform of every element.
     * @param elements the beamline element list, sorted by ascending field-start position
     * @param labFrame  the beamline lab-to-entry transform (OpalBeamline::coordTransformationTo_m)
     */
    static void resolve(FieldList& elements, const CoordinateSystemTrafo& labFrame);
};

#endif  // OPALX_PlacementResolver_HH
