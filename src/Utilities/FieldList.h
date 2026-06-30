#ifndef OPALX_FieldList_HH
#define OPALX_FieldList_HH

//
// FieldList
//   The beamline's list of placed elements: a plain list of element pointers.
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

#include <list>
#include <memory>

class ElementBase;

/// The tracking-side list of placed beamline elements. Membership order is the s-sorted
/// order established by OpalBeamline::prepareSections(); the active s-range of each element
/// lives in the derived IndexMap, not here.
using FieldList = std::list<std::shared_ptr<ElementBase>>;

#endif  // OPALX_FieldList_HH
