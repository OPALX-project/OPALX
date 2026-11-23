//
// Function writeEMFieldDumps
//   Writes the external fields of the tracked line on the grids of the DUMPEMFIELDS requests.
//
// Copyright (c) 2017, Chris Rogers
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
#ifndef OPAL_EM_FIELD_DUMP_WRITER_H
#define OPAL_EM_FIELD_DUMP_WRITER_H

#include <vector>

#include "AbsBeamline/ElementBase.h"
#include "Track/RunDescription.h"

/**
 * @brief Write the summed E and B fields of @p elements on the grid of each request, one file
 *        per request; an existing file is overwritten.
 *
 * Each element is evaluated in its local frame and its fields are rotated back to the global
 * frame. Along the rows t varies fastest and x slowest. File format:
 *
 *     <number of grid points>
 *     1  x [m]
 *     2  y [m]
 *     3  z [m]
 *     4  t [ns]
 *     5  Bx [T]
 *     6  By [T]
 *     7  Bz [T]
 *     8  Ex [MV/m]
 *     9  Ey [MV/m]
 *     10 Ez [MV/m]
 *     0
 *     <x> <y> <z> <t> <Bx> <By> <Bz> <Ex> <Ey> <Ez>
 *     ...
 *
 * @throws OpalException if a file cannot be opened or written.
 */
void writeEMFieldDumps(
        const std::vector<opalx::run::EMFieldDumpDescription>& dumps, const ElementSet& elements);

#endif  // OPAL_EM_FIELD_DUMP_WRITER_H
