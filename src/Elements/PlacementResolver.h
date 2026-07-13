#ifndef OPALX_PlacementResolver_HH
#define OPALX_PlacementResolver_HH

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
#include "Utilities/ElementList.h"

/**
 * @class PlacementResolver
 * @brief Resolves every beamline element's global-to-local transform in one place.
 *
 *  There are two ways to place elements in OPALX
 *   - 6D positioning (Mode A): The global X/Y/Z/PHI/PSI/THETA coordinates of the
 *     elements.
 *   - ELEMEDGE (Mode B): The elements are placed by their position along the
 *     reference path S.
 *
 *
 * @note These conventions cannot be mixed.
 * @note The the global-to-local transforms always transform to the elements'
 * geometrical entrance frame.
 */
class PlacementResolver {
public:
    /**
     * @brief Resolve and fix the global-to-local transform of every element.
     * @param elements the beamline element list, sorted by ascending field-start position
     * @param labFrame  the beamline lab-to-entry transform (OpalBeamline::coordTransformationTo_m)
     */
    static void resolve(ElementList& elements, const CoordinateSystemTrafo& labFrame);
};

#endif  // OPALX_PlacementResolver_HH
